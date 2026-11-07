// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Reverse-mode differentiation through a powerflow the outer loops solved: one transposed
// solve of the last Newton's Jacobian, every decision of the loops held -- on a single grid
// (LSGrid::solve_JT) and per row of a batch (keep_jacobian). Checked against what knows
// nothing about any of it: a central finite difference of the loss, the loops re-run, on
// steps small enough that no decision of theirs flips.

#include <cmath>
#include <complex>
#include <functional>
#include <memory>
#include <tuple>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "Solvers.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "powerflow_algorithm/outer_loop/DistributedSlackLoop.hpp"
#include "powerflow_algorithm/outer_loop/ReactiveLimitsLoop.hpp"
#include "powerflow_algorithm/outer_loop/TransformerVoltageControlLoop.hpp"

using Catch::Approx;
using ls2g::BatchAdjoint;
using ls2g::CplxVect;
using ls2g::DistributedSlackLoop;
using ls2g::InjectionSweep;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::ReactiveLimitsLoop;
using ls2g::RealVect;
using ls2g::TimeSeries;
using ls2g::TransformerVoltageControlLoop;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using AdjointMat = LSGrid::AdjointMat;
using RealMatRM = BatchAdjoint::RealMatRM;
const real_type SN_MVA = 100.;
const real_type TOL = 1e-12;

// The inputs a test perturbs, in MW / MVar / pu, one entry per load / generator.
struct Inputs
{
    std::vector<real_type> load_p, load_q, gen_p, gen_v;
};

// A loss depending on both the angle and the magnitude of every bus voltage; its cotangent
// in torch's convention (dL/dRe V + i dL/dIm V) is the constant wr + i wi.
real_type wr(int bus) { return 1. + 0.5 * bus; }
real_type wi(int bus) { return -0.75 + 0.25 * bus; }

real_type loss_of(const CplxVect & V)
{
    real_type res = 0.;
    for(Eigen::Index b = 0; b < V.size(); ++b) {
        res += wr(static_cast<int>(b)) * V(b).real() + wi(static_cast<int>(b)) * V(b).imag();
    }
    return res;
}

CplxVect flat(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
}

// A grid, ready to solve: built from the inputs, its algorithm and loops set.
using Builder = std::function<LSGrid(const Inputs &)>;

real_type run_loss(const Builder & make, const Inputs & in)
{
    LSGrid grid = make(in);
    const CplxVect V = grid.ac_pf(flat(grid), 30, TOL);
    REQUIRE(V.size() == static_cast<Eigen::Index>(grid.total_bus()));
    return loss_of(V);
}

struct Grad
{
    std::vector<real_type> load_p, load_q, gen_p, gen_v;
};

// The whole gradient of loss_of through the grid's last solve, as a caller assembles it:
// the cotangent onto the unknowns, one transposed solve, lambda read at each input's row,
// plus what the grid says beyond that (gen_v_indirect_grad, gen_p_indirect_grad) and the
// direct half of gen_v.
Grad grad_single(LSGrid & grid, const CplxVect & V,
                 const std::vector<int> & load_bus, const std::vector<int> & gen_bus)
{
    const int dim = grid.dim_J();
    const IntVect theta_col = grid.get_theta_col_of_bus();
    const IntVect vm_col = grid.get_vm_col_of_bus();
    AdjointMat xbar = AdjointMat::Zero(1, dim);
    for(Eigen::Index b = 0; b < V.size(); ++b) {
        const cplx_type gV(wr(static_cast<int>(b)), wi(static_cast<int>(b)));
        if(theta_col[b] >= 0) xbar(0, theta_col[b]) = (std::conj(V(b)) * gV).imag();
        if(vm_col[b] >= 0) xbar(0, vm_col[b]) = (std::conj(V(b) / std::abs(V(b))) * gV).real();
    }
    const RealVect lambda = grid.solve_JT(xbar).row(0).transpose();
    const IntVect p_row = grid.get_p_row_of_bus();
    const IntVect q_row = grid.get_q_row_of_bus();
    const auto at = [&lambda](const IntVect & rows, int bus) { return rows[bus] >= 0 ? lambda(rows[bus]) : 0.; };

    Grad res;
    for(int bus : load_bus) {
        res.load_p.push_back(-at(p_row, bus) / SN_MVA);
        res.load_q.push_back(-at(q_row, bus) / SN_MVA);
    }
    const RealVect gen_p_extra = grid.gen_p_indirect_grad(lambda);
    const RealVect gen_v_ind = grid.gen_v_indirect_grad(lambda);
    const IntVect target = grid.get_gen_v_target_bus();
    const RealVect share = grid.get_gen_v_share();
    for(std::size_t g = 0; g < gen_bus.size(); ++g) {
        const Eigen::Index k = static_cast<Eigen::Index>(g);
        res.gen_p.push_back(at(p_row, gen_bus[g]) / SN_MVA + gen_p_extra(k));
        real_type gv = gen_v_ind(k);
        if(target(k) >= 0) {
            const int b = target(k);
            const cplx_type gV(wr(b), wi(b));
            gv += (std::conj(V(b) / std::abs(V(b))) * gV).real() * share(k);
        }
        res.gen_v.push_back(gv);
    }
    return res;
}

// central differences of the loss, one input at a time
real_type fd(const Builder & make, const Inputs & in, std::vector<real_type> Inputs::* which,
             std::size_t k, real_type delta)
{
    Inputs plus = in, minus = in;
    (plus.*which)[k] += delta;
    (minus.*which)[k] -= delta;
    return (run_loss(make, plus) - run_loss(make, minus)) / (2. * delta);
}

// every input's gradient against its finite difference
void check_against_fd(const Builder & make, const Inputs & in, const Grad & grad,
                      real_type margin_pq = 1e-7, real_type margin_v = 2e-5)
{
    for(std::size_t k = 0; k < in.load_p.size(); ++k) {
        CHECK(grad.load_p[k] == Approx(fd(make, in, &Inputs::load_p, k, 1e-4)).margin(margin_pq));
        CHECK(grad.load_q[k] == Approx(fd(make, in, &Inputs::load_q, k, 1e-4)).margin(margin_pq));
    }
    for(std::size_t k = 0; k < in.gen_p.size(); ++k) {
        CHECK(grad.gen_p[k] == Approx(fd(make, in, &Inputs::gen_p, k, 1e-4)).margin(margin_pq));
        CHECK(grad.gen_v[k] == Approx(fd(make, in, &Inputs::gen_v, k, 1e-6)).margin(margin_v));
    }
}

// ---- a feeder whose reactive limits bite ------------------------------------------------
// buses 0-1-2-3 in a line, a load at bus 3; the slack generator at bus 0, a generator at
// bus 1 well inside its limits (it stays PV, its bus pinned) and one at bus 2 that cannot
// hold its set-point (frozen at its max)
const std::vector<int> RL_LOAD_BUS = {3};
const std::vector<int> RL_GEN_BUS = {0, 1, 2};

Inputs rl_inputs()
{
    Inputs in;
    in.load_p = {60.};
    in.load_q = {40.};
    in.gen_p = {0., 15., 20.};
    in.gen_v = {1.02, 1.015, 1.05};
    return in;
}

LSGrid rl_grid(const Inputs & in, const std::vector<real_type> & lim_q = {500., 500., 10.})
{
    LSGrid grid;
    grid.set_sn_mva(SN_MVA);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(4, 1, RealVect::Constant(4, 138.), 0, 0);
    Eigen::VectorXi from_id(3), to_id(3);
    from_id << 0, 1, 2;
    to_id << 1, 2, 3;
    grid.init_powerlines(RealVect::Constant(3, 0.01), RealVect::Constant(3, 0.08),
                         CplxVect::Zero(3), from_id, to_id);
    RealVect load_p(1), load_q(1);
    load_p << in.load_p[0];
    load_q << in.load_q[0];
    Eigen::VectorXi load_bus(1);
    load_bus << 3;
    grid.init_loads(load_p, load_q, load_bus);
    RealVect gen_p(3), gen_v(3), q(3), min_q(3), max_q(3);
    gen_p << in.gen_p[0], in.gen_p[1], in.gen_p[2];
    gen_v << in.gen_v[0], in.gen_v[1], in.gen_v[2];
    q << 0., 0., 0.;
    min_q << -lim_q[0], -lim_q[1], -lim_q[2];
    max_q << lim_q[0], lim_q[1], lim_q[2];
    Eigen::VectorXi gen_bus(3);
    gen_bus << 0, 1, 2;
    grid.init_generators_full(gen_p, gen_v, q, {true, true, true}, min_q, max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ReactiveLimitsLoop>());
    return grid;
}

real_type gen_q(const LSGrid & grid, int gen_id) { return std::get<1>(grid.get_gen_res())(gen_id); }
real_type gen_p(const LSGrid & grid, int gen_id) { return std::get<0>(grid.get_gen_res())(gen_id); }

}  // namespace


TEST_CASE("the gradient through the reactive limits is the one a finite difference measures", "[outer][adjoint]")
{
    const Inputs in = rl_inputs();
    const Builder make = [](const Inputs & x) { return rl_grid(x); };
    LSGrid grid = make(in);
    const CplxVect V = grid.ac_pf(flat(grid), 30, TOL);
    REQUIRE(V.size() == 4);
    // the state this is about: bus 2 frozen at its limit, bus 1 still holding its voltage
    REQUIRE(gen_q(grid, 2) == Approx(10.).margin(1e-6));
    REQUIRE(std::abs(V(1)) == Approx(in.gen_v[1]).margin(1e-10));
    REQUIRE(grid.get_vm_col_of_bus()[1] >= 0);     // its reserved magnitude unknown

    const Grad grad = grad_single(grid, V, RL_LOAD_BUS, RL_GEN_BUS);
    check_against_fd(make, in, grad);
    // the frozen unit's set-point reaches nothing; the pinned one's does
    CHECK(grad.gen_v[2] == 0.);
    CHECK(std::abs(grad.gen_v[1]) > 1e-3);
}


TEST_CASE("a unit the distributed slack holds at a bound has no gen_p gradient", "[outer][adjoint]")
{
    // three units sharing the slack inside the Newton (NROuterDist_), generator 1 capped:
    // its share would take it past 55 MW, so it is parked there, whatever its set-point
    Inputs in;
    in.load_p = {200.};
    in.load_q = {20.};
    in.gen_p = {50., 40., 60.};
    in.gen_v = {1.02, 1.01, 1.01};
    const Builder make = [](const Inputs & x) {
        LSGrid grid;
        grid.set_sn_mva(SN_MVA);
        grid.set_init_vm_pu(1.0);
        grid.init_bus(4, 1, RealVect::Constant(4, 138.), 0, 0);
        Eigen::VectorXi from_id(3), to_id(3);
        from_id << 0, 1, 2;
        to_id << 1, 2, 3;
        grid.init_powerlines(RealVect::Constant(3, 0.01), RealVect::Constant(3, 0.05),
                             CplxVect::Zero(3), from_id, to_id);
        RealVect load_p(1), load_q(1);
        load_p << x.load_p[0];
        load_q << x.load_q[0];
        Eigen::VectorXi load_bus(1);
        load_bus << 3;
        grid.init_loads(load_p, load_q, load_bus);
        RealVect gen_p(3), gen_v(3), min_q(3), max_q(3);
        gen_p << x.gen_p[0], x.gen_p[1], x.gen_p[2];
        gen_v << x.gen_v[0], x.gen_v[1], x.gen_v[2];
        min_q << -500., -500., -500.;
        max_q << 500., 500., 500.;
        Eigen::VectorXi gen_bus(3);
        gen_bus << 0, 1, 2;
        grid.init_generators(gen_p, gen_v, min_q, max_q, gen_bus);
        grid.add_gen_slackbus(0, 1.);
        grid.set_gen_can_participate_slack({true, true, true}, RealVect::Ones(3));
        RealVect p_min = RealVect::Zero(3), p_max(3);
        p_max << 1000., 55., 1000.;
        grid.set_gen_p_limits(p_min, p_max);
        grid.change_algorithm("NROuterDist_SparseLU");
        grid.clear_outer_loops();
        grid.add_outer_loop(std::make_shared<DistributedSlackLoop>());
        return grid;
    };
    LSGrid grid = make(in);
    const CplxVect V = grid.ac_pf(flat(grid), 30, TOL);
    REQUIRE(V.size() == 4);
    REQUIRE(gen_p(grid, 1) == Approx(55.).margin(1e-6));        // parked
    REQUIRE(gen_p(grid, 2) > in.gen_p[2] + 1.);                  // the others share

    const Grad grad = grad_single(grid, V, {3}, {0, 1, 2});
    check_against_fd(make, in, grad, 1e-6);
    CHECK(grad.gen_p[1] == Approx(0.).margin(1e-12));
    // the units share what the load takes: little moves, but something does
    CHECK(std::abs(grad.load_p[0]) > 1e-4);
}


TEST_CASE("NROuter_ differentiates its last single-slack Newton", "[outer][adjoint]")
{
    // The DistributedSlack loop around a single-slack Newton moves the units' targets
    // BETWEEN solves: the last Jacobian puts every change on the reference bus. That is
    // the documented gradient of NROuter_*, and it is the one of a plain single-slack solve
    // at the targets the loop settled on.
    Inputs in;
    in.load_p = {200.};
    in.load_q = {20.};
    in.gen_p = {50., 40., 60.};
    in.gen_v = {1.02, 1.01, 1.01};
    const auto make = [](const Inputs & x, bool with_loop) {
        LSGrid grid;
        grid.set_sn_mva(SN_MVA);
        grid.set_init_vm_pu(1.0);
        grid.init_bus(4, 1, RealVect::Constant(4, 138.), 0, 0);
        Eigen::VectorXi from_id(3), to_id(3);
        from_id << 0, 1, 2;
        to_id << 1, 2, 3;
        grid.init_powerlines(RealVect::Constant(3, 0.01), RealVect::Constant(3, 0.05),
                             CplxVect::Zero(3), from_id, to_id);
        RealVect load_p(1), load_q(1);
        load_p << x.load_p[0];
        load_q << x.load_q[0];
        Eigen::VectorXi load_bus(1);
        load_bus << 3;
        grid.init_loads(load_p, load_q, load_bus);
        RealVect gen_p(3), gen_v(3), min_q(3), max_q(3);
        gen_p << x.gen_p[0], x.gen_p[1], x.gen_p[2];
        gen_v << x.gen_v[0], x.gen_v[1], x.gen_v[2];
        min_q << -500., -500., -500.;
        max_q << 500., 500., 500.;
        Eigen::VectorXi gen_bus(3);
        gen_bus << 0, 1, 2;
        grid.init_generators(gen_p, gen_v, min_q, max_q, gen_bus);
        grid.add_gen_slackbus(0, 1.);
        if(with_loop) {
            grid.set_gen_can_participate_slack({true, true, true}, RealVect::Ones(3));
            grid.change_algorithm("NROuter_SparseLU");
            grid.clear_outer_loops();
            grid.add_outer_loop(std::make_shared<DistributedSlackLoop>());
        } else {
            grid.change_algorithm("NRSing_SparseLU");
        }
        return grid;
    };
    LSGrid grid = make(in, true);
    const CplxVect V = grid.ac_pf(flat(grid), 30, TOL);
    REQUIRE(V.size() == 4);
    const Grad grad = grad_single(grid, V, {3}, {0, 1, 2});

    Inputs settled = in;
    for(int g = 0; g < 3; ++g) settled.gen_p[static_cast<std::size_t>(g)] = gen_p(grid, g);
    REQUIRE(settled.gen_p[2] > in.gen_p[2] + 1.);   // the loop did move them
    LSGrid ref = make(settled, false);
    const CplxVect V_ref = ref.ac_pf(flat(ref), 30, TOL);
    REQUIRE((V - V_ref).cwiseAbs().maxCoeff() < 1e-8);
    const Grad grad_ref = grad_single(ref, V_ref, {3}, {0, 1, 2});
    CHECK(grad.load_p[0] == Approx(grad_ref.load_p[0]).margin(1e-8));
    CHECK(grad.load_q[0] == Approx(grad_ref.load_q[0]).margin(1e-8));
    for(std::size_t g = 0; g < 3; ++g) {
        CHECK(grad.gen_p[g] == Approx(grad_ref.gen_p[g]).margin(1e-8));
        CHECK(grad.gen_v[g] == Approx(grad_ref.gen_v[g]).margin(1e-8));
    }
}


TEST_CASE("the slack's set-point is differentiated on the taps the loops settled", "[outer][adjoint]")
{
    // bus 0 (the slack) feeds bus 1 through a transformer whose ratio tap the loop moves to
    // regulate bus 1; the slack's gen_v reaches the system through that transformer, at the
    // tap the loop left, not at the one the grid was given
    Inputs in;
    in.load_p = {40.};
    in.load_q = {20.};
    in.gen_p = {0.};
    in.gen_v = {1.0};
    const Builder make = [](const Inputs & x) {
        LSGrid grid;
        grid.set_sn_mva(SN_MVA);
        grid.set_init_vm_pu(1.0);
        RealVect vn(3);
        vn << 138., 20., 20.;
        grid.init_bus(3, 1, vn, 0, 0);
        Eigen::VectorXi line_from(1), line_to(1);
        line_from << 1;
        line_to << 2;
        grid.init_powerlines(RealVect::Constant(1, 0.02), RealVect::Constant(1, 0.1), CplxVect::Zero(1), line_from, line_to);
        RealVect r(1), xx(1), ratio(1), shift(1);
        r << 0.005;
        xx << 0.08;
        ratio << 1.;
        shift << 0.;
        Eigen::VectorXi b1(1), b2(1);
        b1 << 0;
        b2 << 1;
        grid.init_trafo(r, xx, CplxVect::Zero(1), ratio, shift, {false}, b1, b2, false);
        RealVect load_p(1), load_q(1);
        load_p << x.load_p[0];
        load_q << x.load_q[0];
        Eigen::VectorXi load_bus(1);
        load_bus << 2;
        grid.init_loads(load_p, load_q, load_bus);
        RealVect gen_p(1), gen_v(1), min_q(1), max_q(1);
        gen_p << x.gen_p[0];
        gen_v << x.gen_v[0];
        min_q << -1000.;
        max_q << 1000.;
        Eigen::VectorXi gen_bus(1);
        gen_bus << 0;
        grid.init_generators(gen_p, gen_v, min_q, max_q, gen_bus);
        grid.add_gen_slackbus(0, 1.);
        const int n_steps = 21;
        std::vector<real_type> rho, zero(static_cast<std::size_t>(n_steps), 0.);
        for(int k = 0; k < n_steps; ++k) rho.push_back(0.9 + 0.2 * k / (n_steps - 1));
        grid.set_trafo_ratio_tap_changer(0, 0, (n_steps - 1) / 2, rho, zero, zero, zero, zero);
        grid.set_trafo_ratio_tap_regulation(0, true, 1.03, 0.002, 1);
        grid.change_algorithm("NROuter_SparseLU");
        grid.clear_outer_loops();
        grid.add_outer_loop(std::make_shared<TransformerVoltageControlLoop>());
        return grid;
    };
    LSGrid grid = make(in);
    const CplxVect V = grid.ac_pf(flat(grid), 30, TOL);
    REQUIRE(V.size() == 3);
    REQUIRE(grid.get_trafos()[0].res_ratio_tap_position != 10);   // the tap moved

    const Grad grad = grad_single(grid, V, {2}, {0});
    check_against_fd(make, in, grad);
    CHECK(std::abs(grad.gen_v[0]) > 1e-3);
}


TEST_CASE("a remote regulator held at its limit has no gen_v gradient", "[outer][adjoint]")
{
    // generator 1 stands on bus 1 and regulates bus 2 remotely, but it cannot: held at its
    // max, its group's voltage row pins its output, and its set-point reaches nothing
    Inputs in = rl_inputs();
    in.gen_v = {1.02, 1.06, 1.05};
    const Builder make = [](const Inputs & x) {
        LSGrid grid = rl_grid(x, {500., 5., 500.});
        // generator 2 off the voltage problem: gen 1 regulates its bus remotely instead
        grid.deactivate_gen(2);
        grid.set_gen_regulated_bus(1, 2);
        return grid;
    };
    LSGrid grid = make(in);
    const CplxVect V = grid.ac_pf(flat(grid), 30, TOL);
    REQUIRE(V.size() == 4);
    REQUIRE(gen_q(grid, 1) == Approx(5.).margin(1e-6));   // held

    const Grad grad = grad_single(grid, V, RL_LOAD_BUS, {0, 1, 2});
    CHECK(grad.gen_v[1] == 0.);
    // the difference of two solves converged to TOL, over a 1e-6 step: zero to ~1e-8
    CHECK(fd(make, in, &Inputs::gen_v, 1, 1e-6) == Approx(0.).margin(1e-7));
    CHECK(grad.load_p[0] == Approx(fd(make, in, &Inputs::load_p, 0, 1e-4)).margin(1e-7));
    CHECK(grad.load_q[0] == Approx(fd(make, in, &Inputs::load_q, 0, 1e-4)).margin(1e-7));
    CHECK(grad.gen_v[0] == Approx(fd(make, in, &Inputs::gen_v, 0, 1e-6)).margin(2e-5));
}


// ---------------------------------------------------------------------------------------
// the batch
// ---------------------------------------------------------------------------------------

namespace {

const int NB_ROWS = 3;

// the rows of the reactive-limits feeder: the load and the units' set-points vary
void rl_rows(RealMatRM & load_p, RealMatRM & load_q, RealMatRM & gen_p_rows, RealMatRM & gen_v_rows)
{
    const Inputs base = rl_inputs();
    load_p.resize(NB_ROWS, 1);
    load_q.resize(NB_ROWS, 1);
    gen_p_rows.resize(NB_ROWS, 3);
    gen_v_rows.resize(NB_ROWS, 3);
    for(int i = 0; i < NB_ROWS; ++i) {
        load_p(i, 0) = base.load_p[0] + 5. * i;
        load_q(i, 0) = base.load_q[0] + 2. * i;
        for(int g = 0; g < 3; ++g) {
            gen_p_rows(i, g) = base.gen_p[static_cast<std::size_t>(g)] + (g > 0 ? 2. * i : 0.);
            gen_v_rows(i, g) = base.gen_v[static_cast<std::size_t>(g)];
        }
    }
}

Inputs row_inputs(const RealMatRM & load_p, const RealMatRM & load_q, const RealMatRM & gen_p_rows,
                  const RealMatRM & gen_v_rows, int i)
{
    Inputs in;
    in.load_p = {load_p(i, 0)};
    in.load_q = {load_q(i, 0)};
    for(int g = 0; g < 3; ++g) {
        in.gen_p.push_back(gen_p_rows(i, g));
        in.gen_v.push_back(gen_v_rows(i, g));
    }
    return in;
}

// the capability curve of generator 2: its max Q rises with its P
void add_curve(LSGrid & grid)
{
    Eigen::VectorXi ids(2);
    ids << 2, 2;
    RealVect p(2), qmin(2), qmax(2);
    p << 0., 100.;
    qmin << -10., -30.;
    qmax << 5., 25.;
    grid.set_gen_capability_curves(ids, p, qmin, qmax);
}

}  // namespace


TEST_CASE("a batch row's gradient is the one of the same solve on a single grid", "[outer][adjoint][batch]")
{
    RealMatRM load_p, load_q, gen_p_rows, gen_v_rows;
    rl_rows(load_p, load_q, gen_p_rows, gen_v_rows);
    LSGrid grid = rl_grid(rl_inputs());
    InjectionSweep sweep(grid);
    sweep.change_algorithm("NROuter_SparseLU");
    sweep.set_keep_jacobian(true);
    sweep.modify_load_p(load_p);
    sweep.modify_load_q(load_q);
    sweep.modify_gen_p(gen_p_rows);
    sweep.modify_gen_v(gen_v_rows);
    sweep.compute(flat(grid), 30, TOL);
    REQUIRE(sweep.get_status() == 1);

    const int dim = sweep.dim_J();
    const IntVect theta_col = sweep.get_theta_col_of_bus();
    const IntVect vm_col = sweep.get_vm_col_of_bus();
    const auto V = sweep.get_voltages();
    RealMatRM xbar = RealMatRM::Zero(NB_ROWS, dim);
    for(int i = 0; i < NB_ROWS; ++i) {
        for(int b = 0; b < 4; ++b) {
            const cplx_type gV(wr(b), wi(b));
            const cplx_type v = V(i, b);
            if(theta_col[b] >= 0) xbar(i, theta_col[b]) = (std::conj(v) * gV).imag();
            if(vm_col[b] >= 0) xbar(i, vm_col[b]) = (std::conj(v / std::abs(v)) * gV).real();
        }
    }
    const RealMatRM lambda = sweep.solve_JT(xbar);
    const RealMatRM gen_v_ind = sweep.gen_v_indirect_grad(lambda);
    const RealMatRM gen_p_ind = sweep.gen_p_indirect_grad(lambda);
    const RealMatRM share = sweep.get_gen_v_share_rows();

    for(int i = 0; i < NB_ROWS; ++i) {
        LSGrid one = rl_grid(row_inputs(load_p, load_q, gen_p_rows, gen_v_rows, i));
        const CplxVect V_one = one.ac_pf(flat(one), 30, TOL);
        REQUIRE(V_one.size() == 4);
        for(int b = 0; b < 4; ++b) REQUIRE(std::abs(V(i, b) - V_one(b)) < 1e-9);
        AdjointMat xbar_one = xbar.row(i);
        REQUIRE(one.dim_J() == dim);
        const RealVect lam_one = one.solve_JT(xbar_one).row(0).transpose();
        for(int k = 0; k < dim; ++k) CHECK(lambda(i, k) == Approx(lam_one(k)).margin(1e-8));
        const RealVect gv = one.gen_v_indirect_grad(lam_one);
        const RealVect gp = one.gen_p_indirect_grad(lam_one);
        const RealVect sh = one.get_gen_v_share();
        for(int g = 0; g < 3; ++g) {
            CHECK(gen_v_ind(i, g) == Approx(gv(g)).margin(1e-8));
            CHECK(gen_p_ind(i, g) == Approx(gp(g)).margin(1e-8));
            CHECK(share(i, g) == Approx(sh(g)).margin(1e-14));
        }
    }
}


TEST_CASE("a unit frozen on its capability curve owes gen_p the slope of its limit", "[outer][adjoint][batch]")
{
    // A batch row that gives generator 2 its target P reads its limits on its curve at that
    // P: frozen at its max, it injects max_q(P), and moving P moves that reactive power too.
    RealMatRM load_p, load_q, gen_p_rows, gen_v_rows;
    rl_rows(load_p, load_q, gen_p_rows, gen_v_rows);
    const auto run = [&](const RealMatRM & gp, bool keep, InjectionSweep * out) {
        LSGrid grid = rl_grid(rl_inputs());
        add_curve(grid);
        InjectionSweep local(grid);
        InjectionSweep & sweep = out != nullptr ? *out : local;
        sweep.change_algorithm("NROuter_SparseLU");
        sweep.set_keep_jacobian(keep);
        sweep.modify_load_p(load_p);
        sweep.modify_load_q(load_q);
        sweep.modify_gen_p(gp);
        sweep.modify_gen_v(gen_v_rows);
        sweep.compute(flat(grid), 30, TOL);
        REQUIRE(sweep.get_status() == 1);
        real_type res = 0.;
        const auto V = sweep.get_voltages();
        for(Eigen::Index i = 0; i < V.rows(); ++i) {
            for(Eigen::Index b = 0; b < V.cols(); ++b) {
                res += wr(static_cast<int>(b)) * V(i, b).real() + wi(static_cast<int>(b)) * V(i, b).imag();
            }
        }
        return res;
    };
    LSGrid grid = rl_grid(rl_inputs());
    add_curve(grid);
    InjectionSweep sweep(grid);
    run(gen_p_rows, true, &sweep);

    const int dim = sweep.dim_J();
    const IntVect theta_col = sweep.get_theta_col_of_bus();
    const IntVect vm_col = sweep.get_vm_col_of_bus();
    const IntVect p_row = sweep.get_p_row_of_bus();
    const auto V = sweep.get_voltages();
    RealMatRM xbar = RealMatRM::Zero(NB_ROWS, dim);
    for(int i = 0; i < NB_ROWS; ++i) {
        for(int b = 0; b < 4; ++b) {
            const cplx_type gV(wr(b), wi(b));
            const cplx_type v = V(i, b);
            if(theta_col[b] >= 0) xbar(i, theta_col[b]) = (std::conj(v) * gV).imag();
            if(vm_col[b] >= 0) xbar(i, vm_col[b]) = (std::conj(v / std::abs(v)) * gV).real();
        }
    }
    const RealMatRM lambda = sweep.solve_JT(xbar);
    const RealMatRM extra = sweep.gen_p_indirect_grad(lambda);
    for(int i = 0; i < NB_ROWS; ++i) {
        const real_type grad = lambda(i, p_row[2]) / SN_MVA + extra(i, 2);
        CHECK(std::abs(extra(i, 2)) > 1e-6);       // the slope does reach it
        const real_type delta = 1e-4;
        RealMatRM plus = gen_p_rows, minus = gen_p_rows;
        plus(i, 2) += delta;
        minus(i, 2) -= delta;
        const real_type fd_i = (run(plus, false, nullptr) - run(minus, false, nullptr)) / (2. * delta);
        CHECK(grad == Approx(fd_i).margin(1e-7));
    }
}


TEST_CASE("keeping the Jacobians does not change what the batch solves", "[outer][adjoint][batch]")
{
    // the forward with keep_jacobian on: the same voltages bit for bit -- a time series too,
    // whose rows resume from the one before, so a store that disturbed the algorithm's
    // state would show in the next row
    RealMatRM load_p, load_q, gen_p_rows, gen_v_rows;
    rl_rows(load_p, load_q, gen_p_rows, gen_v_rows);
    const auto run = [&](bool keep) {
        LSGrid grid = rl_grid(rl_inputs());
        TimeSeries sweep(grid);
        sweep.change_algorithm("NROuter_SparseLU");
        sweep.set_keep_jacobian(keep);
        sweep.modify_load_p(load_p);
        sweep.modify_load_q(load_q);
        sweep.modify_gen_p(gen_p_rows);
        sweep.modify_gen_v(gen_v_rows);
        sweep.compute(flat(grid), 30, TOL);
        REQUIRE(sweep.get_status() == 1);
        return Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic>(sweep.get_voltages());
    };
    const auto off = run(false);
    const auto on = run(true);
    REQUIRE(off.rows() == NB_ROWS);
    CHECK((off - on).cwiseAbs().maxCoeff() == 0.);
}


TEST_CASE("a single grid's gradient needs the grid as its last solve left it", "[outer][adjoint]")
{
    LSGrid grid = rl_grid(rl_inputs());
    REQUIRE_THROWS(grid.dim_J());                    // nothing solved yet
    REQUIRE(grid.ac_pf(flat(grid), 30, TOL).size() == 4);
    const int dim = grid.dim_J();
    REQUIRE(dim > 0);
    REQUIRE(grid.solve_JT(AdjointMat::Ones(1, dim)).cols() == dim);
    REQUIRE_THROWS(grid.solve_JT(AdjointMat::Ones(1, dim + 1)));   // not this system

    grid.change_p_load(0, 61.);                       // changed since: refused
    REQUIRE_THROWS(grid.solve_JT(AdjointMat::Ones(1, dim)));
    REQUIRE(grid.ac_pf(flat(grid), 30, TOL).size() == 4);
    REQUIRE(grid.solve_JT(AdjointMat::Ones(1, dim)).cols() == dim);

    grid.change_p_load(0, 1e7);                       // diverges
    grid.ac_pf(flat(grid), 30, TOL);
    REQUIRE_THROWS(grid.solve_JT(AdjointMat::Ones(1, dim)));

    LSGrid gs = rl_grid(rl_inputs());
    gs.change_algorithm("GaussSeidel");               // no Jacobian at all
    gs.ac_pf(flat(gs), 1000, 1e-8);
    REQUIRE_THROWS(gs.dim_J());
}
