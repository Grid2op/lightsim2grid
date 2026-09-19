// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// The gradient of a loss with respect to a generator's voltage setpoint (modify_gen_v).
//
// It is the one input of a batch whose gradient the adjoint alone cannot give. An
// injection enters the Newton-Raphson system as a right-hand side, so lambda IS its
// gradient; a voltage setpoint instead fixes a bus' magnitude, which the system does
// not solve for -- so it appears neither in x nor in Sbus, and its gradient needs the
// dS/d|V| column the Jacobian deliberately does not store for such a bus.
//
// Which makes it exactly the thing to check against something that knows nothing about
// any of it: a central finite difference of the loss, with the whole sweep re-run.

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "Solvers.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "case_exotic_elements.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::BatchAdjoint;
using ls2g::CplxVect;
using ls2g::InjectionSweep;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using RealMatRM = BatchAdjoint::RealMatRM;

const int NB_BUS = 4;
const int LOAD_BUS = NB_BUS - 1;
const int GEN_BUS = 1;
const real_type SN_MVA = 100.;
const int NB_STEPS = 4;

// The 4-bus radial feeder of test_batch_adjoint.cpp. `nb_gen_on_pv_bus` puts that many
// generators on bus 1, all regulating it: set_vm is last-writer-wins, so only the last
// of them reaches the solve, and only it can carry a gradient.
LSGrid make_grid(int nb_gen_on_pv_bus = 1, int remote_reg_bus = -1)
{
    const int nb_gen = 1 + nb_gen_on_pv_bus;
    LSGrid g;
    g.set_sn_mva(SN_MVA);
    g.set_init_vm_pu(1.0);
    g.init_bus(static_cast<unsigned int>(NB_BUS), 1, RealVect::Constant(NB_BUS, 138.), 0, 0);

    Eigen::VectorXi f(NB_BUS - 1), t(NB_BUS - 1);
    for(int i = 0; i < NB_BUS - 1; ++i){ f(i) = i; t(i) = i + 1; }
    g.init_powerlines(RealVect::Constant(NB_BUS - 1, 0.01), RealVect::Constant(NB_BUS - 1, 0.1),
                      CplxVect::Zero(NB_BUS - 1), f, t);

    RealVect lp(1), lq(1); lp << 50.; lq << 10.;
    Eigen::VectorXi lb(1); lb << LOAD_BUS;
    g.init_loads(lp, lq, lb);

    RealVect p(nb_gen), v(nb_gen), q(nb_gen), qmin(nb_gen), qmax(nb_gen);
    Eigen::VectorXi b(nb_gen);
    p(0) = 0.; v(0) = 1.02; b(0) = 0;                       // the slack
    for(int k = 1; k < nb_gen; ++k){ p(k) = 20. / nb_gen_on_pv_bus; v(k) = 1.01; b(k) = GEN_BUS; }
    q.setZero(); qmin.setConstant(-1000.); qmax.setConstant(1000.);
    g.init_generators_full(p, v, q, std::vector<bool>(static_cast<size_t>(nb_gen), true), qmin, qmax, b);
    g.add_gen_slackbus(0, 1.);
    // the generators of bus 1 regulating another bus instead: a voltage-control group
    if(remote_reg_bus >= 0){
        for(int k = 1; k < nb_gen; ++k) g.set_gen_regulated_bus(k, remote_reg_bus);
    }
    return g;
}

struct Inputs
{
    RealMatRM gen_p, sgen_p, load_p, load_q, gen_v;

    explicit Inputs(int nb_gen)
        : gen_p(NB_STEPS, nb_gen), sgen_p(NB_STEPS, 0),
          load_p(NB_STEPS, 1), load_q(NB_STEPS, 1), gen_v(NB_STEPS, nb_gen)
    {
        for(int i = 0; i < NB_STEPS; ++i){
            const real_type jitter = static_cast<real_type>((i * 37) % 11) / 10.;
            gen_p(i, 0) = 0.;
            for(int k = 1; k < nb_gen; ++k) gen_p(i, k) = (10. + 15. * jitter) / (nb_gen - 1);
            load_p(i, 0) = 30. + 40. * jitter;
            load_q(i, 0) = 5. + 10. * jitter;
            // a different setpoint per row, and a different one for the slack than for
            // the PV bus -- but every generator OF ONE BUS agrees, which a row must do:
            // a bus has one magnitude (see _row_gen_v_conflicts)
            gen_v(i, 0) = 1.02 + 0.004 * i;
            for(int k = 1; k < nb_gen; ++k) gen_v(i, k) = 1.005 + 0.003 * i;
        }
    }
};

// A loss depending on BOTH the angle and the magnitude of every bus voltage; its
// cotangent in torch's convention (dL/dRe V + i dL/dIm V) is the constant wr + i wi.
real_type wr(int bus) { return 1. + 0.5 * bus; }
real_type wi(int bus) { return -0.75 + 0.25 * bus; }

using CplxMatRM = Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

real_type loss_of(const Eigen::Ref<const CplxMatRM> & V)
{
    real_type res = 0.;
    for(Eigen::Index i = 0; i < V.rows(); ++i){
        for(Eigen::Index b = 0; b < V.cols(); ++b){
            res += wr(static_cast<int>(b)) * V(i, b).real() + wi(static_cast<int>(b)) * V(i, b).imag();
        }
    }
    return res;
}

void feed(InjectionSweep & sweep, const Inputs & in)
{
    sweep.modify_gen_p(in.gen_p);
    sweep.modify_sgen_p(in.sgen_p);
    sweep.modify_load_p(in.load_p);
    sweep.modify_load_q(in.load_q);
    sweep.modify_gen_v(in.gen_v);
}

// the whole sweep, from a fresh object, reduced to one number
real_type run_loss(const Inputs & in, int nb_gen_on_pv_bus, int remote_reg_bus = -1)
{
    LSGrid grid = make_grid(nb_gen_on_pv_bus, remote_reg_bus);
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    feed(sweep, in);
    sweep.compute(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-12);
    REQUIRE(sweep.get_status() == 1);
    return loss_of(sweep.get_voltages());
}

// The whole gen_v gradient: the indirect half from the sweep, the direct half here --
// which is the caller's job, and is the same arithmetic that built xbar (see
// BaseBatchSweep::get_gen_v_target_bus).
RealMatRM gen_v_gradient(InjectionSweep & sweep, int nb_gen)
{
    const int dim_J = sweep.dim_J();
    const IntVect theta_col = sweep.get_theta_col_of_bus();
    const IntVect vm_col = sweep.get_vm_col_of_bus();
    const auto V = sweep.get_voltages();

    RealMatRM xbar = RealMatRM::Zero(NB_STEPS, dim_J);
    for(int i = 0; i < NB_STEPS; ++i){
        for(int b = 0; b < NB_BUS; ++b){
            const cplx_type gV(wr(b), wi(b));
            const cplx_type v = V(i, b);
            if(theta_col[b] >= 0) xbar(i, theta_col[b]) = (std::conj(v) * gV).imag();
            if(vm_col[b] >= 0) xbar(i, vm_col[b]) = (std::conj(v / std::abs(v)) * gV).real();
        }
    }

    const RealMatRM lambda = sweep.solve_JT(xbar);
    RealMatRM grad = sweep.gen_v_indirect_grad(lambda);
    REQUIRE(grad.rows() == NB_STEPS);
    REQUIRE(grad.cols() == nb_gen);

    const IntVect target = sweep.get_gen_v_target_bus();
    const RealVect share = sweep.get_gen_v_share();
    for(int i = 0; i < NB_STEPS; ++i){
        for(int g = 0; g < nb_gen; ++g){
            const int b = target[g];
            if(b < 0) continue;
            const cplx_type gV(wr(b), wi(b));
            const cplx_type v = V(i, b);
            // the same share the indirect half already carries: tied set-points hold a
            // share of one derivative, not a partial each (see get_gen_v_share)
            grad(i, g) += (std::conj(v / std::abs(v)) * gV).real() * share[g];
        }
    }
    return grad;
}

}  // namespace


TEST_CASE("the gen_v gradient is the one a finite difference measures")
{
    const int NB_GEN = 2;   // one slack at bus 0, one PV generator at bus 1
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    sweep.set_keep_jacobian(true);

    const Inputs in(NB_GEN);
    feed(sweep, in);
    sweep.compute(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-12);
    REQUIRE(sweep.get_status() == 1);

    // BOTH generators fix a magnitude the solver never revisits: bus 1 is PV, and bus 0
    // is the slack -- neither owns a Vm unknown, which is exactly the condition
    const IntVect target = sweep.get_gen_v_target_bus();
    REQUIRE(target.size() == NB_GEN);
    REQUIRE(target[0] == 0);
    REQUIRE(target[1] == GEN_BUS);

    const RealMatRM grad = gen_v_gradient(sweep, NB_GEN);

    const real_type delta = 1e-6;   // pu
    for(int i = 0; i < NB_STEPS; ++i){
        for(int g = 0; g < NB_GEN; ++g){
            Inputs plus = in, minus = in;
            plus.gen_v(i, g) += delta;
            minus.gen_v(i, g) -= delta;
            const real_type fd = (run_loss(plus, 1) - run_loss(minus, 1)) / (2. * delta);
            REQUIRE(grad(i, g) == Approx(fd).margin(2e-5));
            // and it is not trivially zero, which would pass any margin
            REQUIRE(std::abs(fd) > 1e-3);
        }
    }
}


TEST_CASE("generators sharing a bus share the one derivative that exists")
{
    // Three generators regulate bus 1, so their set-points are TIED: a row must give
    // all three the same, because a bus has one magnitude. The loss is then a function
    // only on the diagonal v1 = v2 = v3, and off it there is nothing to compare against
    // -- such a row is refused, not solved differently. So no partial derivative of one
    // of them exists, and what these three numbers are is not a gradient in the usual
    // sense: what exists is the derivative along the tie, and they sum to it.
    //
    // Split equally, so that a caller driving the group from one parameter recovers it
    // through the chain rule, and one stepping on all three keeps them equal -- and so
    // that the answer does not depend on the order the generators sit in, which is what
    // "whichever set_vm writes last" would have made it depend on.
    const int NB_ON_PV = 3;
    const int NB_GEN = 1 + NB_ON_PV;
    LSGrid grid = make_grid(NB_ON_PV);
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    sweep.set_keep_jacobian(true);

    const Inputs in(NB_GEN);
    feed(sweep, in);
    sweep.compute(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-12);
    REQUIRE(sweep.get_status() == 1);

    const IntVect target = sweep.get_gen_v_target_bus();
    const RealVect share = sweep.get_gen_v_share();
    REQUIRE(target[0] == 0);                 // the slack, alone on its bus ...
    REQUIRE(share[0] == Approx(1.));         // ... so it owns its derivative outright
    for(int g = 1; g < NB_GEN; ++g){
        REQUIRE(target[g] == GEN_BUS);       // all three regulate it, none is privileged
        REQUIRE(share[g] == Approx(1. / static_cast<real_type>(NB_ON_PV)));
    }

    // A member of the group cannot be perturbed ALONE any more -- that is exactly the
    // contradiction the row now refuses -- so the group moves together, and what the
    // finite difference measures is the SUM of its gradients. That sum is the only part
    // of the split the constraint leaves meaningful, and it is what a chain rule through
    // a tied set of inputs consumes.
    const RealMatRM grad = gen_v_gradient(sweep, NB_GEN);
    const real_type delta = 1e-6;
    for(int i = 0; i < NB_STEPS; ++i){
        Inputs plus = in, minus = in;
        for(int g = 1; g < NB_GEN; ++g){      // every generator of bus 1, together
            plus.gen_v(i, g) += delta;
            minus.gen_v(i, g) -= delta;
        }
        const real_type fd = (run_loss(plus, NB_ON_PV) - run_loss(minus, NB_ON_PV)) / (2. * delta);
        real_type sum = 0.;
        for(int g = 1; g < NB_GEN; ++g) sum += grad(i, g);
        REQUIRE(sum == Approx(fd).margin(2e-5));
        REQUIRE(std::abs(fd) > 1e-3);
        // ... split equally among them, rather than heaped on whichever set_vm writes
        // last: nothing distinguishes them, and only the sum is a derivative at all
        for(int g = 1; g < NB_GEN; ++g){
            REQUIRE(grad(i, g) == Approx(fd / static_cast<real_type>(NB_ON_PV)).margin(2e-5));
        }

        // the slack, alone on its bus, is still perturbed on its own
        Inputs sp = in, sm = in;
        sp.gen_v(i, 0) += delta;
        sm.gen_v(i, 0) -= delta;
        const real_type fd0 = (run_loss(sp, NB_ON_PV) - run_loss(sm, NB_ON_PV)) / (2. * delta);
        REQUIRE(grad(i, 0) == Approx(fd0).margin(2e-5));
    }
}


TEST_CASE("a remote regulator's gen_v is its voltage-control group's set-point")
{
    // Generator 1 stands on bus 1 but regulates bus 2: the bordered voltage control.
    // Bus 2 keeps its magnitude unknown, held by  |V_2| - v_set = 0, so re-seeding |V_2|
    // only moves the starting point -- the set-point that row reads has to be the
    // row's. Each row must therefore land where a one-off solve given that target does.
    const int NB_GEN = 2;
    const int REG_BUS = 2;
    LSGrid grid = make_grid(1, REG_BUS);
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    sweep.set_keep_jacobian(true);

    const Inputs in(NB_GEN);
    feed(sweep, in);
    sweep.compute(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-12);
    REQUIRE(sweep.get_status() == 1);

    const auto V = sweep.get_voltages();
    for(int i = 0; i < NB_STEPS; ++i){
        REQUIRE(std::abs(V(i, REG_BUS)) == Approx(in.gen_v(i, 1)).margin(1e-9));

        LSGrid one_off = make_grid(1, REG_BUS);
        one_off.change_v_gen(0, in.gen_v(i, 0));
        one_off.change_v_gen(1, in.gen_v(i, 1));
        one_off.change_p_gen(1, in.gen_p(i, 1));
        one_off.change_p_load(0, in.load_p(i, 0));
        one_off.change_q_load(0, in.load_q(i, 0));
        one_off.change_algorithm(AlgorithmType::NR_SparseLU);
        const CplxVect ref = one_off.ac_pf(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-12);
        REQUIRE(ref.size() == NB_BUS);
        for(int b = 0; b < NB_BUS; ++b){
            REQUIRE(std::abs(V(i, b) - ref(b)) < 1e-9);
        }
    }

    // no |V| is fixed by it -- its gradient is read at the group's voltage row instead
    REQUIRE(sweep.get_gen_v_target_bus()[1] == -1);
    REQUIRE(sweep.get_gen_v_vc_row()[1] >= 0);
    REQUIRE(sweep.get_gen_v_share()[1] == Approx(1.));

    const RealMatRM grad = gen_v_gradient(sweep, NB_GEN);
    const real_type delta = 1e-6;
    for(int i = 0; i < NB_STEPS; ++i){
        for(int g = 0; g < NB_GEN; ++g){
            Inputs plus = in, minus = in;
            plus.gen_v(i, g) += delta;
            minus.gen_v(i, g) -= delta;
            const real_type fd = (run_loss(plus, 1, REG_BUS) - run_loss(minus, 1, REG_BUS)) / (2. * delta);
            REQUIRE(grad(i, g) == Approx(fd).margin(2e-5));
            REQUIRE(std::abs(fd) > 1e-3);
        }
    }
}


TEST_CASE("a row asking one bus for two magnitudes is not solved")
{
    // |V| at a bus is unique: two generators regulating it with two different targets
    // state two constraints that cannot both hold. set_vm resolves that by applying
    // whichever generator it visits last -- silently, and with no reason to prefer
    // either. Such a row has an unsatisfiable input, so it is not solved at all.
    const int NB_ON_PV = 2;
    const int NB_GEN = 1 + NB_ON_PV;
    LSGrid grid = make_grid(NB_ON_PV);
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);

    Inputs in(NB_GEN);
    // row 1 alone disagrees with itself: generators 1 and 2 both regulate bus 1
    in.gen_v(1, 2) = in.gen_v(1, 1) + 0.01;
    feed(sweep, in);
    sweep.compute(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-12);

    REQUIRE(sweep.converged_mask()[1] == 0);
    // and only that row: its neighbours, which agree, are solved as usual
    for(int i = 0; i < NB_STEPS; ++i){
        if(i == 1) continue;
        REQUIRE(sweep.converged_mask()[static_cast<size_t>(i)] == 1);
    }
    // a row that is not solved reads back as exact zero, like any skipped row
    const auto V = sweep.get_voltages();
    for(int b = 0; b < NB_BUS; ++b) REQUIRE(V(1, b) == cplx_type(0., 0.));
}


TEST_CASE("agreeing set-points on one bus are not a conflict")
{
    // the ordinary case: several machines on a busbar, all asked for the same
    // magnitude. Nothing to refuse -- and the tolerance is the float one, so a value
    // that differs only in the last bits still agrees.
    const int NB_ON_PV = 3;
    const int NB_GEN = 1 + NB_ON_PV;
    LSGrid grid = make_grid(NB_ON_PV);
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);

    Inputs in(NB_GEN);
    in.gen_v(2, 2) = in.gen_v(2, 1) * (1. + 1e-15);
    feed(sweep, in);
    sweep.compute(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-12);
    REQUIRE(sweep.get_status() == 1);
    for(int i = 0; i < NB_STEPS; ++i) REQUIRE(sweep.converged_mask()[static_cast<size_t>(i)] == 1);
}


TEST_CASE("a row that did not converge gets no gen_v gradient")
{
    const int NB_GEN = 2;
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    sweep.set_keep_jacobian(true);

    Inputs in(NB_GEN);
    in.load_p(2, 0) = 1e7;   // far past the nose of the PV curve: this row diverges
    feed(sweep, in);
    sweep.compute(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-12);
    REQUIRE(sweep.converged_mask()[2] == 0);

    const int dim_J = sweep.dim_J();
    const RealMatRM lambda = sweep.solve_JT(RealMatRM::Constant(NB_STEPS, dim_J, 1.));
    const RealMatRM grad = sweep.gen_v_indirect_grad(lambda);
    for(int g = 0; g < NB_GEN; ++g) REQUIRE(grad(2, g) == 0.);
}


TEST_CASE("the gen_v gradient is an AC quantity and says so in DC")
{
    const int NB_GEN = 2;
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::DC_SparseLU);

    const Inputs in(NB_GEN);
    feed(sweep, in);
    sweep.compute(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-12);
    REQUIRE_THROWS(sweep.gen_v_indirect_grad(RealMatRM::Zero(NB_STEPS, 1)));
}


// ===================== set-points modify_gen_v cannot move ====================

TEST_CASE("a row cannot move a generator sharing its bus with an hvdc station")
{
    // ONLY GENERATOR set-points vary per row: there is no modify_svc_v and no
    // modify_hvdc_v. So a generator whose regulated bus also carries a converter
    // station (or a voltage-mode SVC) cannot actually be moved -- that element keeps
    // asking the bus for its own, fixed magnitude, and a row that gives the generator
    // a different one asks one bus for two.
    //
    // Left undetected this is the quiet kind of wrong: _apply_step_gen_v writes only
    // the generators, so the station's target -- already applied when the starting
    // voltage was built -- is simply overwritten, and the row solves at a set-point
    // nobody asked for. See the TODO at the top of the changelog.
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    grid.deactivate_svc(0);
    grid.deactivate_storage(0);
    grid.change_bus_gen_python(1, 1);     // generator 1 onto the station's bus
    grid.change_v_gen(1, 1.0);            // ... agreeing with it, so the GRID is valid

    const int nb_gen = 5;                 // the fixture's generator count
    const Eigen::Index nb_bus = static_cast<Eigen::Index>(grid.total_bus());
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);

    RealMatRM gen_v = RealMatRM::Zero(3, nb_gen);
    for(int i = 0; i < 3; ++i){
        for(int g = 0; g < nb_gen; ++g) gen_v(i, g) = 1.0;
    }
    gen_v(1, 1) = 1.04;                   // row 1 alone tries to move generator 1
    sweep.modify_gen_v(gen_v);
    sweep.compute(CplxVect::Constant(nb_bus, cplx_type(1., 0.)), 30, 1e-10);

    REQUIRE(sweep.converged_mask()[1] == 0);   // refused, not silently applied
    REQUIRE(sweep.converged_mask()[0] == 1);
    REQUIRE(sweep.converged_mask()[2] == 1);
}


TEST_CASE("a generator alone on its bus is free to move")
{
    // the same grid without the co-location: nothing else asks that bus for anything,
    // so every row may set whatever it likes
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    grid.deactivate_svc(0);
    grid.deactivate_storage(0);

    const int nb_gen = 5;
    const Eigen::Index nb_bus = static_cast<Eigen::Index>(grid.total_bus());
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);

    RealMatRM gen_v = RealMatRM::Zero(3, nb_gen);
    for(int i = 0; i < 3; ++i){
        for(int g = 0; g < nb_gen; ++g) gen_v(i, g) = 1.0 + 0.01 * i;
    }
    sweep.modify_gen_v(gen_v);
    sweep.compute(CplxVect::Constant(nb_bus, cplx_type(1., 0.)), 30, 1e-10);
    for(int i = 0; i < 3; ++i) REQUIRE(sweep.converged_mask()[static_cast<size_t>(i)] == 1);
}
