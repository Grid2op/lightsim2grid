// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Tests of the Fast-Decoupled powerflow (src/core/powerflow_algorithm/BaseFDPFAlgo.*).
//
// The family had no answer-level coverage at all: it appeared in the suite only
// through test_cache_reuse.cpp (which checks caching, not numbers) and
// test_plugin_registration.cpp (which checks names). What is pinned here is that
// both flavours, XB and BX, land on the SAME solution as Newton-Raphson on the
// same grid -- the property that makes FDPF an alternative solver rather than a
// different model -- plus the (Vm, Va) <-> V bookkeeping BaseFDPFAlgo::
// has_converged performs on every iteration.
//
// That bookkeeping is the reason this file exists. has_converged rebuilds V from
// (Vm, Va) and then has to put the pair back in canonical form -- magnitude
// >= 0, angle wrapped -- because the Q iteration (`Vm_(pq) -= q_`) can drive a
// magnitude negative and the P iteration accumulates into Va without wrapping.
// It used to do that by going through V (`Vm_ = V_.abs(); Va_ = V_.arg()`), a
// hypot and an atan2 per bus; it now does it directly. The two agree on every
// grid that converges, which is what the sections below check, and the identity
// itself is checked head-on in "canonical form ... " at the bottom -- including
// the negative-magnitude and out-of-range-angle cases a converging grid never
// reaches (they only occur on trajectories that end up diverging).
//
// Grids are built by hand through the init_* API, like test_lsgrid.cpp: no
// python, no grid2op, no external file. C++14 only (project policy).

#include <cmath>
#include <complex>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "Solvers.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::BaseConstants;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

// an `nx` x `ny` mesh of identical lines, a load on every bus but bus 0, a
// slack generator at bus 0 and a PV generator every 25 buses. `load_scale`
// stretches the loading, so the same topology can be asked for an easy solve
// or a heavily loaded one.
LSGrid make_mesh(int nx, int ny, real_type load_scale)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    const int n_bus = nx * ny;
    grid.init_bus(static_cast<unsigned int>(n_bus), 1, RealVect::Constant(n_bus, 138.), 0, 0);

    std::vector<int> from, to;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int here = j * nx + i;
            if (i + 1 < nx) { from.push_back(here); to.push_back(here + 1); }
            if (j + 1 < ny) { from.push_back(here); to.push_back(here + nx); }
        }
    }
    const int n_line = static_cast<int>(from.size());
    Eigen::VectorXi from_id(n_line), to_id(n_line);
    for (int k = 0; k < n_line; ++k) { from_id(k) = from[k]; to_id(k) = to[k]; }
    grid.init_powerlines(RealVect::Constant(n_line, 0.002),
                         RealVect::Constant(n_line, 0.02),
                         CplxVect::Constant(n_line, cplx_type(0., 0.02)),
                         from_id, to_id);

    const int n_load = n_bus - 1;
    RealVect load_p(n_load), load_q(n_load);
    Eigen::VectorXi load_bus(n_load);
    for (int k = 0; k < n_load; ++k) {
        load_bus(k) = k + 1;
        load_p(k) = load_scale * (1. + 0.01 * (k % 17));
        load_q(k) = load_scale * (0.2 + 0.005 * (k % 7));
    }
    grid.init_loads(load_p, load_q, load_bus);

    std::vector<int> gen_buses;
    gen_buses.push_back(0);
    for (int b = 25; b < n_bus; b += 25) gen_buses.push_back(b);
    const int n_gen = static_cast<int>(gen_buses.size());
    RealVect gen_p(n_gen), gen_v(n_gen), gen_qmin(n_gen), gen_qmax(n_gen);
    Eigen::VectorXi gen_bus(n_gen);
    for (int k = 0; k < n_gen; ++k) {
        gen_bus(k) = gen_buses[k];
        gen_p(k) = (k == 0) ? 0. : 20. * load_scale;
        gen_v(k) = 1.02;
        gen_qmin(k) = -1e5;
        gen_qmax(k) = 1e5;
    }
    grid.init_generators(gen_p, gen_v, gen_qmin, gen_qmax, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    return grid;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
}

}  // namespace

TEST_CASE("FDPF converges to the Newton-Raphson solution", "[fdpf]")
{
    // FDPF is an approximation of the JACOBIAN, not of the equations: it takes
    // more iterations than Newton-Raphson but must land on the same voltages.
    // Anything that quietly changed what FDPF solves -- rather than how fast it
    // gets there -- shows up here.
    const int sizes[][2] = {{4, 4}, {7, 7}, {11, 11}};
    const real_type loads[] = {0.2, 1.0, 3.0};

    for (const auto & size : sizes) {
        for (real_type load : loads) {
            const std::string what = std::to_string(size[0]) + "x" + std::to_string(size[1])
                                   + " load " + std::to_string(load);
            LSGrid nr_grid = make_mesh(size[0], size[1], load);
            nr_grid.change_algorithm(AlgorithmType::NR_SparseLU);
            const CplxVect v_ref = nr_grid.ac_pf(flat_start(nr_grid), 50, 1e-11);
            INFO("newton-raphson reference: " << what);
            REQUIRE(v_ref.size() > 0);

            for (AlgorithmType algo : {AlgorithmType::FDPF_XB_SparseLU,
                                       AlgorithmType::FDPF_BX_SparseLU}) {
                LSGrid grid = make_mesh(size[0], size[1], load);
                grid.init_fdpf_coeffs();
                grid.change_algorithm(algo);
                const CplxVect v = grid.ac_pf(flat_start(grid), 100, 1e-9);
                INFO("fdpf: " << what);
                REQUIRE(v.size() == v_ref.size());
                for (Eigen::Index i = 0; i < v.size(); ++i) {
                    INFO("bus " << i);
                    // both are converged solutions of the same equations, to
                    // 1e-11 and 1e-9 respectively: they agree far tighter than
                    // the looser of the two tolerances
                    CHECK(std::abs(v(i) - v_ref(i)) < 1e-8);
                }
            }
        }
    }
}

TEST_CASE("FDPF reports magnitudes and angles consistent with its own voltages", "[fdpf]")
{
    // get_Vm() / get_Va() are not free-standing outputs: they are the state
    // has_converged canonicalises every iteration, and get_V() is rebuilt from
    // them. If the canonicalisation drifted from the voltages, this catches it.
    for (AlgorithmType algo : {AlgorithmType::FDPF_XB_SparseLU,
                               AlgorithmType::FDPF_BX_SparseLU}) {
        LSGrid grid = make_mesh(7, 7, 1.0);
        grid.init_fdpf_coeffs();
        grid.change_algorithm(algo);
        const CplxVect v = grid.ac_pf(flat_start(grid), 100, 1e-9);
        REQUIRE(v.size() > 0);
        const RealVect vm = grid.get_algo().get_Vm();
        const RealVect va = grid.get_algo().get_Va();
        REQUIRE(vm.size() == v.size());
        REQUIRE(va.size() == v.size());
        for (Eigen::Index i = 0; i < v.size(); ++i) {
            INFO("bus " << i);
            CHECK(vm(i) >= 0.);                                  // canonical magnitude
            CHECK(std::abs(va(i)) <= BaseConstants::my_pi + 1e-12);  // canonical angle
            CHECK(vm(i) == Approx(std::abs(v(i))).epsilon(1e-12));
            CHECK(std::cos(va(i)) == Approx(std::cos(std::arg(v(i)))).margin(1e-12));
            CHECK(std::sin(va(i)) == Approx(std::sin(std::arg(v(i)))).margin(1e-12));
        }
    }
}

TEST_CASE("canonical form of (Vm, Va) matches the abs/arg it replaces", "[fdpf]")
{
    // BaseFDPFAlgo::canonicalise_vm_va replaced `Vm_ = V_.abs(); Va_ = V_.arg();`
    // in has_converged. Here the real function is called (not a copy of its
    // formula) and checked against the pair it replaced, over a domain a
    // converging powerflow never visits: a negative magnitude (the Q iteration
    // overshooting) and an angle several turns out of range (the P iteration
    // accumulating without wrapping). Those states only occur on trajectories
    // that end in divergence, and a diverged solve clears Vm_ / Va_ / V_ -- so
    // this is the only level at which the branches can be observed at all.
    using FDPF = ls2g::BaseFDPFAlgo<ls2g::SparseLULinearSolver, ls2g::FDPFMethod::XB>;
    const real_type pi = BaseConstants::my_pi;

    // one batch, so the function is exercised on a real vector rather than
    // coefficient by coefficient
    std::vector<real_type> vm_in, va_in;
    for (int im = -60; im <= 60; ++im) {
        for (int ia = -400; ia <= 400; ++ia) {
            vm_in.push_back(0.05 * im);   // negative, zero and positive
            va_in.push_back(0.025 * ia);  // |Va| up to 10 rad: several wraps
        }
    }
    const Eigen::Index n = static_cast<Eigen::Index>(vm_in.size());
    RealVect vm(n), va(n);
    for (Eigen::Index i = 0; i < n; ++i) { vm(i) = vm_in[i]; va(i) = va_in[i]; }

    // what the abs()/arg() pair produced, from the same starting point
    CplxVect v_ref(n);
    RealVect vm_ref(n), va_ref(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        const cplx_type tmp_va(std::cos(va(i)), std::sin(va(i)));  // as compute_pf builds it
        v_ref(i) = vm(i) * tmp_va;
        vm_ref(i) = std::abs(v_ref(i));
        va_ref(i) = std::arg(v_ref(i));
    }

    FDPF::canonicalise_vm_va(vm, va);

    for (Eigen::Index i = 0; i < n; ++i) {
        INFO("Vm = " << vm_in[i] << ", Va = " << va_in[i]);
        CHECK(vm(i) >= 0.);                    // canonical magnitude
        CHECK(std::abs(va(i)) <= pi + 1e-12);  // canonical angle
        CHECK(vm(i) == Approx(vm_ref(i)).margin(1e-15));
        if (vm_in[i] != 0.) {
            // angles compared as phasors: +pi and -pi are the same direction,
            // and which one atan2 returns at the branch cut is not part of the
            // contract
            CHECK(std::cos(va(i)) == Approx(std::cos(va_ref(i))).margin(1e-12));
            CHECK(std::sin(va(i)) == Approx(std::sin(va_ref(i))).margin(1e-12));
        }
        // At Vm == 0 the phase of a zero phasor is undefined and the two
        // deliberately disagree (see canonicalise_vm_va). What still has to
        // hold everywhere is the voltage itself -- the only thing the solve
        // goes on to use.
        CHECK(std::abs(vm(i) * cplx_type(std::cos(va(i)), std::sin(va(i))) - v_ref(i)) < 1e-14);
    }
}
