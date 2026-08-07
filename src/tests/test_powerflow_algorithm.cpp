// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Dedicated tests of the Newton-Raphson hot path (src/core/powerflow_algorithm)
// and of its extensions -- hvdc angle droop, remote voltage control / SVC,
// distributed slack -- driven entirely from C++.
//
// Why a separate file from test_lsgrid.cpp: those tests check that a grid gives
// the right *answer*. These check that the NRSystem machinery underneath stays
// *internally consistent*, which is a different property and one that only
// shows up as memory corruption when it breaks. Release wheels are built
// -O3 -DNDEBUG, so every index the extensions use -- `ledger.p_row(bus)` is a
// `std::vector::operator[]`, `Va(bus)` an Eigen `operator()` -- is unchecked;
// the guarantee that they are in range comes from validation, not from the
// language. Run this file under ASan/UBSan (see .github/workflows/sanitizers.yml)
// and an out-of-bounds access becomes a failure instead of silent corruption.
//
// Two families here:
//   * end-to-end solves through LSGrid with the extensions engaged, alone and
//     together, plus long mutate-and-resolve loops -- the sanitizer bait;
//   * direct drives of NRSystem's 3-phase protocol, which is the only way to
//     reach the states LSGrid itself refuses to build (a topology rebuild that
//     never happened, a Jacobian smaller than the grid).

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "LSGrid.hpp"
#include "powerflow_algorithm/NRSystem.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using ls2g::CplxVect;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::MultiSlackNRSystem;
using ls2g::RealVect;
using ls2g::SingleSlackNRSystem;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

// ---------------------------------------------------------------------------
// A six-bus mesh, roomy enough to host every extension at once without the
// configurations LSGrid::fill_voltage_control_solver_data rejects (a controller
// on a PV bus, a regulated bus that is already voltage-pinned, an SVC sharing
// its regulated bus). Bus roles:
//   0 : slack generator, locally voltage-regulating (the reference)
//   1 : generator, remote voltage controller (regulates bus 3) -> stays PQ
//   2 : load, hvdc side 1
//   3 : load, regulated by the generator at bus 1
//   4 : load, SVC host
//   5 : load, hvdc side 2
// Lines: 0-1, 1-2, 2-3, 3-4, 4-5, 5-0 (a ring, so no radial dead end).
// ---------------------------------------------------------------------------
LSGrid make_six_bus_skeleton()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    RealVect bus_vn_kv(6);
    bus_vn_kv << 138., 138., 138., 138., 138., 138.;
    grid.init_bus(6, 1, bus_vn_kv, 0, 0);

    const int nb_line = 6;
    RealVect branch_r(nb_line), branch_x(nb_line);
    branch_r << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01;
    branch_x << 0.10, 0.10, 0.10, 0.10, 0.10, 0.10;
    const CplxVect branch_h = CplxVect::Zero(nb_line);
    Eigen::VectorXi from_id(nb_line), to_id(nb_line);
    from_id << 0, 1, 2, 3, 4, 5;
    to_id   << 1, 2, 3, 4, 5, 0;
    grid.init_powerlines(branch_r, branch_x, branch_h, from_id, to_id);

    // one load per non-generator bus
    RealVect load_p(4), load_q(4);
    load_p << 20., 25., 15., 30.;
    load_q << 5., 6., 4., 8.;
    Eigen::VectorXi load_bus(4);
    load_bus << 2, 3, 4, 5;
    grid.init_loads(load_p, load_q, load_bus);

    // gen 0 at bus 0 (slack, local PV), gen 1 at bus 1
    RealVect gen_p(2), gen_v(2), gen_min_q(2), gen_max_q(2);
    gen_p << 0., 30.;
    gen_v << 1.02, 1.03;
    gen_min_q << -1000., -1000.;
    gen_max_q << 1000., 1000.;
    Eigen::VectorXi gen_bus(2);
    gen_bus << 0, 1;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    return grid;
}

// gen 1 regulates bus 3 instead of its own bus -> VoltageControl group of one
void add_remote_voltage_control(LSGrid & grid)
{
    grid.set_gen_regulated_bus(1, 3);
    grid.tell_solver_need_reset();
}

// a voltage-mode SVC at bus 4 regulating its own bus
void add_svc(LSGrid & grid, real_type target_vm = 1.01, real_type slope_pu = 0.)
{
    RealVect target(1), q_set(1), slope(1), b_min(1), b_max(1);
    target << target_vm;
    q_set << 0.;
    slope << slope_pu;
    b_min << -100.;
    b_max << 100.;
    Eigen::VectorXi reg_bus(1), svc_bus(1);
    reg_bus << 4;
    svc_bus << 4;
    // SvcContainer::RegulationMode: VOLTAGE = 1
    grid.init_svcs({1}, target, q_set, slope, b_min, b_max, reg_bus, svc_bus);
    grid.tell_solver_need_reset();
}

// a VSC-VSC hvdc line bus 2 -> bus 5, angle-droop ("AC emulation") enabled.
// Neither station regulates voltage, so both ends stay PQ.
// `pmax_mw` is what a SATURATED regime (status_droop = +/-1) transfers, so it
// must stay physical for this grid (~90 MW of load): a 300 MW forced injection
// on a 6-bus ring simply has no solution and would make the test measure
// divergence rather than memory safety.
void add_droop_hvdc(LSGrid & grid, real_type p0_mw = 10., real_type k_mw_per_deg = 2.,
                    real_type pmax_mw = 40.)
{
    Eigen::VectorXi bus1(1), bus2(1);
    bus1 << 2;
    bus2 << 5;
    const std::vector<int> type_vsc{0};          // ConverterType: VSC = 0
    const std::vector<int> mode{0};              // ConvertersMode: SIDE_1_RECTIFIER = 0
    const std::vector<bool> vreg{false}, droop_on{true};
    const RealVect zero = RealVect::Zero(1);
    RealVect vm_pu(1), q_min(1), q_max(1), pf(1), p_set(1), p_max(1), p0(1), k_deg(1);
    vm_pu << 1.0;
    q_min << -1e6;
    q_max << 1e6;
    pf << 1.;
    p_set << 0.;
    p_max << pmax_mw;
    p0 << p0_mw;
    k_deg << k_mw_per_deg;
    grid.init_hvdc_lines(bus1, bus2, type_vsc, type_vsc,
                         zero, zero,             // no station losses
                         vreg, vreg, vm_pu, vm_pu,
                         zero, zero,             // q setpoints
                         q_min, q_max, q_min, q_max,
                         pf, pf, mode, p_set,
                         zero, zero,             // r_ohm, nominal_v
                         droop_on, p0, k_deg, p_max, p_max);
    grid.tell_solver_need_reset();
}

// An AlgoControl with every flag raised: "assume everything changed". The
// extensions cache their per-solve grid data and only re-pull it when the
// control says something moved (see nr_extension_data_is_stale), so a test
// driving update_state by hand must say so explicitly -- otherwise it would be
// measuring the cache, not the pull.
const ls2g::AlgoControl & all_changed()
{
    static ls2g::AlgoControl c = []{ ls2g::AlgoControl a; a.tell_all_changed(); return a; }();
    return c;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
}

CplxVect solve_ac(LSGrid & grid, int max_iter = 30)
{
    return grid.ac_pf(flat_start(grid), max_iter, 1e-9);
}

bool all_finite(const CplxVect & v)
{
    for (Eigen::Index i = 0; i < v.size(); ++i)
        if (!std::isfinite(v(i).real()) || !std::isfinite(v(i).imag())) return false;
    return true;
}

// Every entry of a bus -> J index map is either "this bus owns none" (-1) or a
// real index of a dim x dim Jacobian. Anything else would be used raw.
// Takes the vector by const reference: these maps are returned BY VALUE by the
// *_python accessors, so the caller must bind the result to a named variable
// before indexing it.
void check_index_map(const IntVect & map, int dim)
{
    for (Eigen::Index i = 0; i < map.size(); ++i) {
        REQUIRE(map(i) >= -1);
        REQUIRE(map(i) < dim);
    }
}

}  // anonymous namespace


// ===========================================================================
// Family 1: end-to-end through LSGrid, extensions engaged
// ===========================================================================

TEST_CASE("NR hot path: the six-bus reference grid solves with no extension", "[nr][baseline]")
{
    LSGrid grid = make_six_bus_skeleton();
    const CplxVect V = solve_ac(grid);
    REQUIRE(V.size() == 6);
    CHECK(all_finite(V));
    // bus 0 is the slack (1.02 pu), bus 1 an ordinary PV bus (1.03 pu)
    CHECK(std::abs(V(0)) == Approx(1.02).epsilon(1e-8));
    CHECK(std::abs(V(1)) == Approx(1.03).epsilon(1e-8));
}

TEST_CASE("NR hot path: each extension alone keeps the solve finite", "[nr][extensions]")
{
    SECTION("remote voltage control") {
        LSGrid grid = make_six_bus_skeleton();
        add_remote_voltage_control(grid);
        const CplxVect V = solve_ac(grid);
        REQUIRE(V.size() == 6);
        REQUIRE(all_finite(V));
        // the bordered voltage row pins the REGULATED bus, not the controller's
        CHECK(std::abs(V(3)) == Approx(1.03).epsilon(1e-7));
    }
    SECTION("SVC in voltage mode") {
        LSGrid grid = make_six_bus_skeleton();
        add_svc(grid, 1.01);
        const CplxVect V = solve_ac(grid);
        REQUIRE(V.size() == 6);
        REQUIRE(all_finite(V));
        CHECK(std::abs(V(4)) == Approx(1.01).epsilon(1e-7));
    }
    SECTION("hvdc angle droop") {
        LSGrid grid = make_six_bus_skeleton();
        add_droop_hvdc(grid);
        const CplxVect V = solve_ac(grid);
        REQUIRE(V.size() == 6);
        REQUIRE(all_finite(V));
        // linear regime: the transferred power follows p0 + k.(theta1 - theta2)
        const real_type k_per_rad = 2. * 180. / std::acos(real_type(-1.));
        const real_type expected = 10. + k_per_rad * (std::arg(V(2)) - std::arg(V(5)));
        CHECK(std::get<0>(grid.get_dcline_res2())(0) == Approx(expected).epsilon(1e-6));
    }
    SECTION("distributed slack") {
        LSGrid grid = make_six_bus_skeleton();
        grid.add_gen_slackbus(1, 0.4);
        grid.tell_solver_need_reset();
        const CplxVect V = solve_ac(grid);
        REQUIRE(V.size() == 6);
        REQUIRE(all_finite(V));
    }
}

TEST_CASE("NR hot path: all extensions engaged at once", "[nr][extensions]")
{
    // hvdc droop (buses 2-5) + remote voltage control (gen 1 -> bus 3) + SVC
    // (bus 4) + distributed slack, all in the same Jacobian. The bordered
    // VoltageControl block, MultiSlack's extra column and Hvdc's dP/dtheta
    // entries all land in one J here, which is the layout no single-feature
    // test exercises.
    LSGrid grid = make_six_bus_skeleton();
    add_droop_hvdc(grid);
    add_svc(grid, 1.01);
    add_remote_voltage_control(grid);
    grid.add_gen_slackbus(1, 0.4);
    grid.tell_solver_need_reset();

    const CplxVect V = solve_ac(grid);
    REQUIRE(V.size() == 6);
    REQUIRE(all_finite(V));
    CHECK(std::abs(V(3)) == Approx(1.03).epsilon(1e-6));  // remote target
    CHECK(std::abs(V(4)) == Approx(1.01).epsilon(1e-6));  // svc target

    // the converged reactive injections come back one per controller
    const RealVect q = grid.get_algo().get_controller_q();
    const IntVect kind = grid.get_algo().get_controller_kind();
    REQUIRE(q.size() == kind.size());
    REQUIRE(q.size() == 2);  // the remote generator and the SVC
    for (Eigen::Index i = 0; i < q.size(); ++i) CHECK(std::isfinite(q(i)));
}

TEST_CASE("NR hot path: the Jacobian layout stays self-consistent", "[nr][jacobian]")
{
    LSGrid grid = make_six_bus_skeleton();
    add_droop_hvdc(grid);
    add_svc(grid, 1.01);
    add_remote_voltage_control(grid);
    grid.add_gen_slackbus(1, 0.4);
    grid.tell_solver_need_reset();
    REQUIRE(solve_ac(grid).size() == 6);

    const auto & algo = grid.get_algo();
    const Eigen::SparseMatrix<real_type> J = algo.get_J_python();
    const int dim = static_cast<int>(J.rows());

    // square, and the augmented dimension is what the ledger handed out
    CHECK(J.cols() == J.rows());
    CHECK(dim > 0);

    // every bus -> row/column map has one entry per bus and points inside J
    // (or says "none")
    const IntVect theta_col = algo.get_theta_to_J_col_python();
    const IntVect vm_col    = algo.get_vm_to_J_col_python();
    const IntVect q_col     = algo.get_q_to_J_col_python();
    const IntVect p_row     = algo.get_p_to_J_row_python();
    const IntVect q_row     = algo.get_q_to_J_row_python();
    REQUIRE(theta_col.size() == 6);
    REQUIRE(vm_col.size() == 6);
    REQUIRE(q_col.size() == 6);
    REQUIRE(p_row.size() == 6);
    REQUIRE(q_row.size() == 6);
    check_index_map(theta_col, dim);
    check_index_map(vm_col, dim);
    check_index_map(q_col, dim);
    check_index_map(p_row, dim);
    check_index_map(q_row, dim);

    // the compact registration pair lists agree in length and stay in range
    const IntVect p_buses = algo.get_p_buses_python(), p_rows = algo.get_p_rows_python();
    const IntVect q_buses = algo.get_q_buses_python(), q_rows = algo.get_q_rows_python();
    REQUIRE(p_buses.size() == p_rows.size());
    REQUIRE(q_buses.size() == q_rows.size());
    for (Eigen::Index i = 0; i < p_rows.size(); ++i) {
        CHECK(p_buses(i) >= 0); CHECK(p_buses(i) < 6);
        CHECK(p_rows(i)  >= 0); CHECK(p_rows(i)  < dim);
    }
    for (Eigen::Index i = 0; i < q_rows.size(); ++i) {
        CHECK(q_buses(i) >= 0); CHECK(q_buses(i) < 6);
        CHECK(q_rows(i)  >= 0); CHECK(q_rows(i)  < dim);
    }

    // MultiSlack's slack_absorbed column, and VoltageControl's per-controller
    // Q columns, are dereferenced unconditionally on every NR iteration
    CHECK(algo.get_slack_col() >= 0);
    CHECK(algo.get_slack_col() < dim);
    const IntVect q_cols = algo.get_controller_q_col();
    for (Eigen::Index i = 0; i < q_cols.size(); ++i) {
        CHECK(q_cols(i) >= 0);
        CHECK(q_cols(i) < dim);
    }
    // J is never left with a non-finite coefficient after a converged solve
    for (int k = 0; k < J.nonZeros(); ++k) CHECK(std::isfinite(J.valuePtr()[k]));
}

TEST_CASE("NR hot path: repeated solves under grid mutation", "[nr][hotpath]")
{
    // The scenario the extensions' cached row/column arrays are exposed to in
    // production: the same solver object re-solving while the grid keeps
    // changing underneath. Each mutation may or may not trigger a topology
    // rebuild; either way the solve must end in a converged result or a clean
    // exception -- never in an out-of-bounds access. Worth running under ASan.
    LSGrid grid = make_six_bus_skeleton();
    add_droop_hvdc(grid);
    add_svc(grid, 1.01);
    add_remote_voltage_control(grid);
    grid.add_gen_slackbus(1, 0.4);
    grid.tell_solver_need_reset();
    REQUIRE(solve_ac(grid).size() == 6);

    int nb_solved = 0;
    for (int step = 0; step < 40; ++step) {
        // a deterministic but varied mix of the mutations grid2op drives
        switch (step % 8) {
            case 0: grid.change_p_load(0, 20. + step); break;
            case 1: grid.change_q_load(1, 5. + 0.5 * step); break;
            case 2: grid.change_p_gen(1, 10. + step); break;
            case 3: grid.change_v_gen(0, 1.01 + 0.001 * (step % 5)); break;
            // droop regime flip: an INPUT of the solver, and the one change the
            // Hvdc extension declares its entries unconditionally to absorb
            // without altering the J sparsity pattern
            case 4: grid.set_status_droop_hvdc(0, (step / 8) % 3 - 1); break;
            case 5: grid.deactivate_load(2); break;
            case 6: grid.reactivate_load(2); break;
            // Re-point the remote controller between two legal PQ buses. Both
            // must be electrically reachable from the controller at bus 1: bus 5
            // sits next to the slack, so pinning it to 1.03 against the slack's
            // own 1.02 across one 0.1 pu reactance has no solution and would
            // make this loop measure divergence instead of memory safety.
            case 7: grid.set_gen_regulated_bus(1, (step % 16 == 7) ? 3 : 2); break;
            default: break;
        }
        const CplxVect V = solve_ac(grid);
        if (V.size() != 0) {
            CHECK(all_finite(V));
            ++nb_solved;
        } else {
            WARN("step " << step << " (case " << (step % 8) << ") did not converge");
        }
    }
    // The sequence above is physically reasonable throughout, so every step
    // should converge. Asserting the full count (rather than "it did not
    // crash") is what keeps this a test of the hot path: a mutation that
    // quietly stopped converging would stop exercising the extensions at all.
    CHECK(nb_solved == 40);
}

TEST_CASE("NR hot path: connecting a droop line forces a topology rebuild", "[nr][hvdc]")
{
    // The set the droop equations apply to is {droop_enabled AND connected}, and
    // the Hvdc extension sizes its cached row / column arrays from it at
    // register_in time (topology rebuild) while re-pulling the data itself on
    // every solve. Disconnecting and reconnecting a droop line changes that set
    // in both directions; the GROWING one is what the caches cannot survive
    // without a rebuild (arrays sized for 0 lines, indexed for 1).
    //
    // unset_changes() after every solve, as the grid2op backend does, so the
    // solver really is on its cached path between steps rather than being saved
    // by a stale "everything changed" flag.
    //
    // Measured honestly: this scenario already rebuilds today without
    // AlgoControl::has_hvdc_droop_changed() (reconnecting an hvdc line moves
    // Sbus and trips one of the other rebuild flags), so this test does NOT fail
    // if that flag is removed. It is here as an end-to-end guard on the
    // invariant -- that the droop set and the caches indexing it stay in step
    // across a connect / disconnect cycle -- not as proof of the flag.
    LSGrid grid = make_six_bus_skeleton();
    add_droop_hvdc(grid);
    REQUIRE(solve_ac(grid).size() == 6);
    grid.unset_changes();

    grid.deactivate_dcline(0);          // droop set: 1 -> 0
    REQUIRE(solve_ac(grid).size() == 6);
    grid.unset_changes();

    grid.reactivate_dcline(0);          // droop set: 0 -> 1, arrays still sized 0
    const CplxVect V = solve_ac(grid);
    REQUIRE(V.size() == 6);
    CHECK(all_finite(V));

    // and the droop law is being applied again, i.e. the extension really did
    // re-register rather than silently stay empty
    const real_type k_per_rad = 2. * 180. / std::acos(real_type(-1.));
    const real_type expected = 10. + k_per_rad * (std::arg(V(2)) - std::arg(V(5)));
    CHECK(std::get<0>(grid.get_dcline_res2())(0) == Approx(expected).epsilon(1e-6));
}

TEST_CASE("NR hot path: droop status flips reuse the sparsity pattern", "[nr][hvdc]")
{
    // Hvdc::declare_feature_entries claims its 4 dP/dtheta entries whatever the
    // regime precisely so a `status_droop` flip is a pure value change. If that
    // ever regressed, the pattern would change while the symbolic factorization
    // is reused -- values written at positions resolved for a different matrix.
    LSGrid grid = make_six_bus_skeleton();
    add_droop_hvdc(grid);
    REQUIRE(solve_ac(grid).size() == 6);
    const Eigen::SparseMatrix<real_type> J_linear = grid.get_algo().get_J_python();

    grid.set_status_droop_hvdc(0, 1);  // saturated 1 -> 2
    REQUIRE(solve_ac(grid).size() == 6);
    const Eigen::SparseMatrix<real_type> J_saturated = grid.get_algo().get_J_python();

    REQUIRE(J_saturated.rows() == J_linear.rows());
    REQUIRE(J_saturated.nonZeros() == J_linear.nonZeros());
    for (int k = 0; k <= J_linear.rows(); ++k)
        CHECK(J_saturated.outerIndexPtr()[k] == J_linear.outerIndexPtr()[k]);
    for (int k = 0; k < J_linear.nonZeros(); ++k)
        CHECK(J_saturated.innerIndexPtr()[k] == J_linear.innerIndexPtr()[k]);
}

TEST_CASE("NR hot path: an hvdc droop end on a slack bus drops its entries", "[nr][hvdc]")
{
    // A slack end owns no theta unknown (single slack) and its P balance is not
    // an equation of the base block, so Hvdc::register_in gets -1 back from the
    // ledger and declare_feature_entries must skip those entries rather than
    // hand FeatureSink a negative row/column.
    LSGrid grid = make_six_bus_skeleton();
    Eigen::VectorXi bus1(1), bus2(1);
    bus1 << 0;  // the reference slack bus
    bus2 << 3;
    const std::vector<int> type_vsc{0};
    const std::vector<int> mode{0};
    const std::vector<bool> vreg{false}, droop_on{true};
    const RealVect zero = RealVect::Zero(1);
    RealVect vm_pu(1), q_min(1), q_max(1), pf(1), p_set(1), p_max(1), p0(1), k_deg(1);
    vm_pu << 1.0; q_min << -1e6; q_max << 1e6; pf << 1.;
    p_set << 0.; p_max << 300.; p0 << 10.; k_deg << 2.;
    grid.init_hvdc_lines(bus1, bus2, type_vsc, type_vsc, zero, zero,
                         vreg, vreg, vm_pu, vm_pu, zero, zero,
                         q_min, q_max, q_min, q_max, pf, pf, mode, p_set,
                         zero, zero, droop_on, p0, k_deg, p_max, p_max);
    grid.tell_solver_need_reset();

    const CplxVect V = solve_ac(grid);
    REQUIRE(V.size() == 6);
    CHECK(all_finite(V));
}

TEST_CASE("NR hot path: a sloped SVC borders the voltage row", "[nr][svc]")
{
    // A non-zero slope adds the (v_row, q_col) coupling. The entry is declared
    // for every SVC whatever its slope, so this must not change the pattern
    // either -- only the value, and the converged voltage.
    LSGrid grid = make_six_bus_skeleton();
    add_svc(grid, 1.01, /*slope_pu=*/0.02);
    const CplxVect V = solve_ac(grid);
    REQUIRE(V.size() == 6);
    REQUIRE(all_finite(V));
    // Vm(reg) = v_set - slope . Q  (generator convention): injecting reactive
    // power pulls the regulated bus below its setpoint
    const RealVect q = grid.get_algo().get_controller_q();
    REQUIRE(q.size() == 1);
    CHECK(std::abs(V(4)) == Approx(1.01 - 0.02 * q(0)).epsilon(1e-6));
}

TEST_CASE("NR hot path: two generators sharing one regulated bus", "[nr][vctrl]")
{
    // A control group of two: one bordered voltage row plus one reactive-sharing
    // row, and two Q columns. This is the `cnt - 1` sharing-row arithmetic in
    // VoltageControl::register_in, which no single-controller test reaches.
    // re-declare the generators: a third one at bus 5, also regulating bus 3
    // (init_generators replaces the skeleton's two, so the slack is re-declared
    // below as well)
    LSGrid grid = make_six_bus_skeleton();
    RealVect gen_p(3), gen_v(3), gen_min_q(3), gen_max_q(3);
    gen_p << 0., 30., 20.;
    gen_v << 1.02, 1.03, 1.03;
    gen_min_q << -1000., -1000., -500.;
    gen_max_q << 1000., 1000., 500.;
    Eigen::VectorXi gen_bus(3);
    gen_bus << 0, 1, 5;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    grid.set_gen_regulated_bus(1, 3);
    grid.set_gen_regulated_bus(2, 3);
    grid.tell_solver_need_reset();

    const CplxVect V = solve_ac(grid);
    REQUIRE(V.size() == 6);
    REQUIRE(all_finite(V));
    CHECK(std::abs(V(3)) == Approx(1.03).epsilon(1e-6));

    // proportional-to-range sharing: Q_i / (qmax_i - qmin_i) is the same for both
    const RealVect q = grid.get_algo().get_controller_q();
    REQUIRE(q.size() == 2);
    const real_type w0 = 2000. / 100., w1 = 1000. / 100.;  // pu, sn_mva = 100
    CHECK(q(0) / w0 == Approx(q(1) / w1).epsilon(1e-6));
}


// ===========================================================================
// Family 2: direct drives of the NRSystem protocol
// ===========================================================================

namespace {

// Snapshot of everything compute_pf hands an NRSystem, taken from a grid that
// has already solved once (so its solver-side labelling exists). Held by value:
// NRSystem caches Ybus by reference and Sbus by raw pointer, so these must
// outlive every phase call -- which is exactly the lifetime contract the
// production caller (LSGrid, holding Ybus_ac_ / acSbus_ as members) satisfies.
struct SolverInputs
{
    Eigen::SparseMatrix<cplx_type> Ybus;
    CplxVect V;
    CplxVect Sbus;
    IntVect slack_ids;
    RealVect slack_weights;
    IntVect pv;
    IntVect pq;

    explicit SolverInputs(const LSGrid & grid):
        Ybus(grid.get_Ybus_solver()),
        V(grid.get_V_solver()),
        Sbus(grid.get_Sbus_solver()),
        slack_ids(grid.get_slack_ids_solver_numpy()),
        slack_weights(grid.get_slack_weights_solver()),
        pv(grid.get_pv_solver_numpy()),
        pq(grid.get_pq_solver_numpy())
    {
        Ybus.makeCompressed();
    }
};

// the full rebuild path, in the order NRAlgo::compute_pf runs it
void run_all_phases(MultiSlackNRSystem & sys, const LSGrid & grid, const SolverInputs & in)
{
    sys.update_state(&grid, in.Ybus, in.V, in.Sbus, in.slack_weights, all_changed());
    sys.init_topology(in.slack_ids, in.slack_weights, in.pv, in.pq);
    sys.build_J_sparsity();
    sys.validate_registration();
}

}  // anonymous namespace

TEST_CASE("NRSystem: the 3-phase protocol builds a square, fully resolved J", "[nr][nrsystem]")
{
    LSGrid grid = make_six_bus_skeleton();
    add_droop_hvdc(grid);
    add_svc(grid, 1.01);
    add_remote_voltage_control(grid);
    grid.add_gen_slackbus(1, 0.4);
    grid.tell_solver_need_reset();
    REQUIRE(solve_ac(grid).size() == 6);

    const SolverInputs in(grid);
    MultiSlackNRSystem sys;
    REQUIRE_NOTHROW(run_all_phases(sys, grid, in));

    const int dim = static_cast<int>(sys.total_state_variables());
    CHECK(sys.J().rows() == dim);
    CHECK(sys.J().cols() == dim);

    // filling J must not touch anything outside it: under ASan a stale feature
    // handle or an unresolved (-1) position would fire right here
    REQUIRE_NOTHROW(sys.fill_internal_variables());
    REQUIRE_NOTHROW(sys.fill_J());
    for (int k = 0; k < sys.J().nonZeros(); ++k)
        CHECK(std::isfinite(sys.J().valuePtr()[k]));

    // the mismatch is one residual per unknown, and a zero step is a no-op
    const RealVect F = sys.mismatch();
    REQUIRE(F.size() == dim);
    CHECK(F.allFinite());
    const RealVect zero_step = RealVect::Zero(dim);
    CHECK(sys.mismatch_sq_norm_at(zero_step) == Approx(F.squaredNorm()).epsilon(1e-12));
    CHECK(sys.max_abs_dtheta(zero_step) == Approx(0.));
    CHECK(sys.max_abs_dvm(zero_step) == Approx(0.));
}

TEST_CASE("NRSystem: refreshing the data without a rebuild is caught", "[nr][nrsystem][validation]")
{
    // THE regression this validation exists for. update_state runs on every
    // solve and re-pulls the extension data from the grid; register_in only
    // re-derives the cached row/column arrays on a topology rebuild. If a grid
    // change alters the number of hvdc droop lines or voltage controllers
    // without signalling that rebuild through AlgoControl, the two drift apart
    // and the per-iteration loops index the stale arrays out of bounds --
    // through FeatureWriter, an out-of-bounds WRITE into J.valuePtr().
    //
    // No LSGrid API reaches that state today (every setter that changes a
    // controller count calls tell_pv_changed), which is exactly why it has to
    // be reproduced by driving the protocol by hand: the guard has to hold even
    // if some future setter forgets.
    LSGrid grid = make_six_bus_skeleton();
    add_droop_hvdc(grid);
    add_remote_voltage_control(grid);
    grid.tell_solver_need_reset();
    REQUIRE(solve_ac(grid).size() == 6);

    const SolverInputs in(grid);

    SECTION("a voltage controller appears") {
        MultiSlackNRSystem sys;
        REQUIRE_NOTHROW(run_all_phases(sys, grid, in));

        // add an SVC: one more controller, one more group. A real rebuild would
        // follow; here we deliberately skip it, as a missing tell_pv_changed()
        // would.
        add_svc(grid, 1.01);
        sys.update_state(&grid, in.Ybus, in.V, in.Sbus, in.slack_weights, all_changed());
        REQUIRE_THROWS_MATCHES(
            sys.validate_registration(), std::runtime_error,
            MessageMatches(ContainsSubstring("VoltageControl")));
    }

    SECTION("a voltage controller disappears") {
        MultiSlackNRSystem sys;
        REQUIRE_NOTHROW(run_all_phases(sys, grid, in));

        grid.set_gen_regulated_bus(1, 1);  // back to local PV: no controller left
        sys.update_state(&grid, in.Ybus, in.V, in.Sbus, in.slack_weights, all_changed());
        REQUIRE_THROWS_MATCHES(
            sys.validate_registration(), std::runtime_error,
            MessageMatches(ContainsSubstring("VoltageControl")));
    }

    SECTION("an hvdc droop line disappears") {
        MultiSlackNRSystem sys;
        REQUIRE_NOTHROW(run_all_phases(sys, grid, in));

        grid.deactivate_dcline(0);  // no longer connected: dropped from the data
        sys.update_state(&grid, in.Ybus, in.V, in.Sbus, in.slack_weights, all_changed());
        REQUIRE_THROWS_MATCHES(
            sys.validate_registration(), std::runtime_error,
            MessageMatches(ContainsSubstring("Hvdc")));
    }
}

TEST_CASE("NRSystem: extension bus ids are checked against the solver size", "[nr][nrsystem][validation]")
{
    // The extension data never travels through compute_pf's arguments, so
    // BaseAlgo::check_pf_inputs -- which range-checks Ybus / Sbus / pv / pq --
    // never sees it. Handing the system a Jacobian smaller than the grid makes
    // the grid's own (valid) bus ids out of range for it, which is the shape of
    // every way this can go wrong: it must raise, not index past the ledger.
    //
    // This is the DEEP validation: per-element, and therefore opt-in
    // (set_deep_validation), so that it runs on the python-facing entry point
    // and not inside a ContingencyAnalysis / TimeSeries loop. The staleness
    // check that does run on every solve is exercised by the test above.
    LSGrid grid = make_six_bus_skeleton();
    add_droop_hvdc(grid);       // buses 2 and 5
    add_svc(grid, 1.01);        // bus 4
    grid.tell_solver_need_reset();
    REQUIRE(solve_ac(grid).size() == 6);

    const SolverInputs full(grid);

    // a 3-bus solver view: buses 3, 4 and 5 no longer exist for it
    const int small = 3;
    SolverInputs in(grid);
    in.Ybus = full.Ybus.topLeftCorner(small, small);
    in.Ybus.makeCompressed();
    in.V = full.V.head(small);
    in.Sbus = full.Sbus.head(small);
    in.slack_weights = full.slack_weights.head(small);
    in.slack_ids = IntVect::Zero(1);          // bus 0 is the slack
    in.pv = IntVect();
    in.pq = IntVect(2);
    in.pq << 1, 2;

    MultiSlackNRSystem sys;
    sys.set_deep_validation(true);
    sys.update_state(&grid, in.Ybus, in.V, in.Sbus, in.slack_weights, all_changed());
    // init_topology validates before register_in indexes the ledger with them
    REQUIRE_THROWS_MATCHES(
        sys.init_topology(in.slack_ids, in.slack_weights, in.pv, in.pq),
        std::runtime_error,
        MessageMatches(ContainsSubstring("outside the valid range")));
}

TEST_CASE("NRSystem: mismatched voltage / Sbus sizes are rejected", "[nr][nrsystem][validation]")
{
    LSGrid grid = make_six_bus_skeleton();
    REQUIRE(solve_ac(grid).size() == 6);
    const SolverInputs full(grid);

    SECTION("V shorter than Ybus") {
        SolverInputs in(grid);
        const CplxVect V_short = full.V.head(4);
        MultiSlackNRSystem sys;
        sys.set_deep_validation(true);
        sys.update_state(&grid, in.Ybus, V_short, in.Sbus, in.slack_weights, all_changed());
        REQUIRE_THROWS_MATCHES(
            sys.init_topology(in.slack_ids, in.slack_weights, in.pv, in.pq),
            std::runtime_error,
            MessageMatches(ContainsSubstring("voltage vectors")));
    }
    SECTION("Sbus shorter than Ybus") {
        SolverInputs in(grid);
        const CplxVect S_short = full.Sbus.head(4);
        MultiSlackNRSystem sys;
        sys.set_deep_validation(true);
        sys.update_state(&grid, in.Ybus, in.V, S_short, in.slack_weights, all_changed());
        REQUIRE_THROWS_MATCHES(
            sys.init_topology(in.slack_ids, in.slack_weights, in.pv, in.pq),
            std::runtime_error,
            MessageMatches(ContainsSubstring("Sbus")));
    }
}

TEST_CASE("NRSystem: a system with no extension data is a strict no-op", "[nr][nrsystem]")
{
    // The documented promise of both extensions: with zero droop lines and zero
    // controllers, the ledger / J / results are identical to a build without
    // them. Compare the single- and multi-slack systems on a single-slack grid.
    LSGrid grid = make_six_bus_skeleton();
    REQUIRE(solve_ac(grid).size() == 6);
    const SolverInputs in(grid);

    SingleSlackNRSystem sys;
    sys.update_state(&grid, in.Ybus, in.V, in.Sbus, in.slack_weights, all_changed());
    sys.init_topology(in.slack_ids, in.slack_weights, in.pv, in.pq);
    sys.build_J_sparsity();
    REQUIRE_NOTHROW(sys.validate_registration());

    // 5 non-slack buses: bus 1 is PV (theta only), buses 2..5 PQ (theta + vm)
    CHECK(sys.total_state_variables() == 9u);
    CHECK(sys.controller_q().size() == 0);
    CHECK(sys.slack_col() == -1);
}

TEST_CASE("NRSystem: bus masking rewrites rows without touching the pattern", "[nr][nrsystem][masking]")
{
    // set_masked_buses is a pure value-level edit used by ContingencyAnalysis to
    // keep J non-singular when a contingency isolates buses. It resolves
    // J.valuePtr() positions from the (fixed) sparsity, so a stale mask would
    // write at positions of a matrix that no longer exists.
    LSGrid grid = make_six_bus_skeleton();
    REQUIRE(solve_ac(grid).size() == 6);
    const SolverInputs in(grid);

    MultiSlackNRSystem sys;
    sys.update_state(&grid, in.Ybus, in.V, in.Sbus, in.slack_weights, all_changed());
    sys.init_topology(in.slack_ids, in.slack_weights, in.pv, in.pq);
    sys.build_J_sparsity();
    sys.fill_internal_variables();
    sys.fill_J();
    const Eigen::SparseMatrix<real_type> J_unmasked = sys.J();

    // mask bus 4 (never the reference slack, per set_masked_buses' contract)
    sys.set_masked_buses({4});
    sys.fill_J();
    const Eigen::SparseMatrix<real_type> J_masked = sys.J();

    REQUIRE(J_masked.nonZeros() == J_unmasked.nonZeros());
    for (int k = 0; k <= J_unmasked.rows(); ++k)
        CHECK(J_masked.outerIndexPtr()[k] == J_unmasked.outerIndexPtr()[k]);
    for (int k = 0; k < J_masked.nonZeros(); ++k)
        CHECK(std::isfinite(J_masked.valuePtr()[k]));

    // bus 4's P row is now the identity: a single 1 on its theta column
    const std::vector<int> & p_row_of_bus = sys.p_to_J_row();
    const std::vector<int> & theta_col_of_bus = sys.theta_to_J_col();
    const int masked_row = p_row_of_bus[4];
    const int masked_col = theta_col_of_bus[4];
    REQUIRE(masked_row >= 0);
    REQUIRE(masked_col >= 0);
    const Eigen::MatrixXd dense = Eigen::MatrixXd(J_masked);
    for (int col = 0; col < dense.cols(); ++col)
        CHECK(dense(masked_row, col) == Approx(col == masked_col ? 1. : 0.));

    // an empty mask restores the unmasked values bit-for-bit
    sys.set_masked_buses({});
    sys.fill_J();
    for (int k = 0; k < J_unmasked.nonZeros(); ++k)
        CHECK(sys.J().valuePtr()[k] == J_unmasked.valuePtr()[k]);
}

TEST_CASE("NRSystem: clear_jacobian leaves a reusable, empty system", "[nr][nrsystem]")
{
    // NRAlgo::reset() goes through this. Everything the extensions cached must
    // be dropped, so the next solve rebuilds from scratch rather than reusing
    // indices of a Jacobian that no longer exists.
    LSGrid grid = make_six_bus_skeleton();
    add_droop_hvdc(grid);
    add_svc(grid, 1.01);
    grid.tell_solver_need_reset();
    REQUIRE(solve_ac(grid).size() == 6);
    const SolverInputs in(grid);

    MultiSlackNRSystem sys;
    REQUIRE_NOTHROW(run_all_phases(sys, grid, in));
    REQUIRE(sys.total_state_variables() > 0u);

    sys.clear_jacobian();
    CHECK(sys.total_state_variables() == 0u);
    CHECK(sys.J().nonZeros() == 0);
    CHECK(sys.controller_q().size() == 0);

    // and it can be driven again from scratch
    REQUIRE_NOTHROW(run_all_phases(sys, grid, in));
    CHECK(sys.total_state_variables() > 0u);
    REQUIRE_NOTHROW(sys.fill_internal_variables());
    REQUIRE_NOTHROW(sys.fill_J());
}

TEST_CASE("NRSystem: a null grid pointer yields empty extension data", "[nr][nrsystem]")
{
    // External / batched solvers drive an NRSystem without an LSGrid. Both
    // extensions must then behave as absent rather than dereference the pointer.
    LSGrid grid = make_six_bus_skeleton();
    REQUIRE(solve_ac(grid).size() == 6);
    const SolverInputs in(grid);

    MultiSlackNRSystem sys;
    sys.update_state(nullptr, in.Ybus, in.V, in.Sbus, in.slack_weights, all_changed());
    REQUIRE_NOTHROW(sys.init_topology(in.slack_ids, in.slack_weights, in.pv, in.pq));
    REQUIRE_NOTHROW(sys.build_J_sparsity());
    REQUIRE_NOTHROW(sys.validate_registration());
    CHECK(sys.controller_q().size() == 0);
    REQUIRE_NOTHROW(sys.fill_internal_variables());
    REQUIRE_NOTHROW(sys.fill_J());
}
