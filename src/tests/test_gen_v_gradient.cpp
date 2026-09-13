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
LSGrid make_grid(int nb_gen_on_pv_bus = 1)
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
            // a different setpoint per row and per generator, so no two are confused
            gen_v(i, 0) = 1.02 + 0.004 * i;
            for(int k = 1; k < nb_gen; ++k) gen_v(i, k) = 1.005 + 0.003 * i + 0.002 * k;
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
real_type run_loss(const Inputs & in, int nb_gen_on_pv_bus)
{
    LSGrid grid = make_grid(nb_gen_on_pv_bus);
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
    for(int i = 0; i < NB_STEPS; ++i){
        for(int g = 0; g < nb_gen; ++g){
            const int b = target[g];
            if(b < 0) continue;
            const cplx_type gV(wr(b), wi(b));
            const cplx_type v = V(i, b);
            grad(i, g) += (std::conj(v / std::abs(v)) * gV).real();
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


TEST_CASE("only the generator whose setpoint reaches the solve carries its gradient")
{
    // Three generators regulate bus 1. set_vm walks them in order and the last one
    // wins, so the first two never reach the solve at all: their setpoints are
    // genuinely dead inputs, and a gradient spread over all three -- or given to the
    // wrong one -- would be silently wrong on any grid with more than one machine per
    // busbar, which is most of them.
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
    REQUIRE(target[0] == 0);                 // the slack, alone on its bus
    REQUIRE(target[1] == -1);                // overwritten by generator 2
    REQUIRE(target[2] == -1);                // overwritten by generator 3
    REQUIRE(target[3] == GEN_BUS);           // the last writer

    const RealMatRM grad = gen_v_gradient(sweep, NB_GEN);
    const real_type delta = 1e-6;
    for(int i = 0; i < NB_STEPS; ++i){
        for(int g = 0; g < NB_GEN; ++g){
            Inputs plus = in, minus = in;
            plus.gen_v(i, g) += delta;
            minus.gen_v(i, g) -= delta;
            const real_type fd = (run_loss(plus, NB_ON_PV) - run_loss(minus, NB_ON_PV)) / (2. * delta);
            REQUIRE(grad(i, g) == Approx(fd).margin(2e-5));
            if(target[g] < 0) REQUIRE(fd == Approx(0.).margin(1e-9));   // a dead input
            else REQUIRE(std::abs(fd) > 1e-3);
        }
    }
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
