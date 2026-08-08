// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Tests of the Newton-Raphson *system* (src/core/powerflow_algorithm/NRSystem*)
// rather than of the answers a grid produces: test_lsgrid.cpp already checks
// that a grid with an hvdc droop / a remote-regulating generator / a
// distributed slack converges to the right numbers. What is checked here is
// that the machinery underneath stays internally consistent -- the augmented
// Jacobian is square and matches the ledger, every bus -> row / column map
// stays in range, the feature (non dS-derived) entries are all resolved, and
// the sparsity pattern only changes when the topology does (so a linear
// solver may reuse its symbolic factorization).
//
// The systems are driven directly, phase by phase, exactly the way
// NRAlgo::compute_pf drives them (update_state -> init_topology ->
// build_J_sparsity -> fill_internal_variables -> fill_J), on the solver-
// labelled data (Ybus / Sbus / pv / pq / slack) a real LSGrid produced. Grids
// are built by hand through the init_* API, like test_lsgrid.cpp: no python,
// no grid2op, no external file. C++14 only (project policy).

#include <cmath>
#include <complex>
#include <set>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "powerflow_algorithm/NRSystem.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::CplxVect;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::MultiSlackNRSystem;
using ls2g::RealVect;
using ls2g::SingleSlackNRSystem;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using SpMatReal = Eigen::SparseMatrix<real_type, Eigen::ColMajor>;

// ---------------------------------------------------------------------------
// hand-built grids
// ---------------------------------------------------------------------------

// `n_bus` buses in a radial feeder 0-1-2-...-(n_bus - 1), identical lines
// (r = 0.01, x = 0.1 pu), one 50 MW / 10 MVar load at the last bus, sn_mva
// = 100. No generator: each test declares the generation / slack it needs.
LSGrid make_radial_skeleton(int n_bus)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    const RealVect bus_vn_kv = RealVect::Constant(n_bus, 138.);
    grid.init_bus(static_cast<unsigned int>(n_bus), 1, bus_vn_kv, 0, 0);

    const int n_line = n_bus - 1;
    const RealVect branch_r = RealVect::Constant(n_line, 0.01);
    const RealVect branch_x = RealVect::Constant(n_line, 0.1);
    const CplxVect branch_h = CplxVect::Zero(n_line);
    Eigen::VectorXi from_id(n_line), to_id(n_line);
    for (int i = 0; i < n_line; ++i) {
        from_id(i) = i;
        to_id(i) = i + 1;
    }
    grid.init_powerlines(branch_r, branch_x, branch_h, from_id, to_id);

    RealVect load_p(1), load_q(1);
    load_p << 50.;
    load_q << 10.;
    Eigen::VectorXi load_bus(1);
    load_bus << n_bus - 1;
    grid.init_loads(load_p, load_q, load_bus);
    return grid;
}

// adds `nb` generators, all at 0 MW, `v_pu` setpoint and +/- `q_range_mvar`/2
// reactive range (the sharing key of the VoltageControl extension is
// qmax - qmin), at the given buses
void add_generators(LSGrid & grid,
                    const std::vector<int> & buses,
                    const std::vector<real_type> & v_pu,
                    const std::vector<real_type> & q_range_mvar)
{
    const int nb = static_cast<int>(buses.size());
    RealVect gen_p(nb), gen_v(nb), gen_min_q(nb), gen_max_q(nb);
    Eigen::VectorXi gen_bus(nb);
    for (int i = 0; i < nb; ++i) {
        gen_p(i) = 0.;
        gen_v(i) = v_pu[i];
        gen_min_q(i) = -0.5 * q_range_mvar[i];
        gen_max_q(i) = 0.5 * q_range_mvar[i];
        gen_bus(i) = buses[i];
    }
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
}

// the radial skeleton plus a single slack generator (1.02 pu) at bus 0
LSGrid make_radial_grid(int n_bus = 3)
{
    LSGrid grid = make_radial_skeleton(n_bus);
    add_generators(grid, {0}, {1.02}, {2000.});
    grid.add_gen_slackbus(0, 1.);
    return grid;
}

// one VSC-VSC hvdc line between `bus_1` and `bus_2`, angle-droop enabled
// (p0 MW, k MW/deg), lossless stations and lossless dc line so the droop
// slopes are exactly +/- k on both sides (loss_mult == 1).
void add_droop_hvdc(LSGrid & grid, int bus_1, int bus_2,
                    real_type droop_p0_mw, real_type droop_k_mw_per_deg,
                    real_type pmax_mw = 300.)
{
    Eigen::VectorXi bus1(1), bus2(1);
    bus1 << bus_1;
    bus2 << bus_2;
    // ConverterStationContainer::ConverterType: VSC = 0
    const std::vector<int> type1{0}, type2{0};
    // HvdcLineContainer::ConvertersMode: SIDE_1_RECTIFIER = 0
    const std::vector<int> mode{0};
    const std::vector<bool> vreg1{false}, vreg2{false}, droop_on{true};
    const RealVect zero = RealVect::Zero(1);
    RealVect vm1_pu(1), vm2_pu(1), q_min(1), q_max(1);
    RealVect pf1(1), pf2(1), p_set(1), p_max(1), p0(1), k_deg(1);
    vm1_pu << 1.0;
    vm2_pu << 1.0;
    q_min << -1e6;
    q_max << 1e6;
    pf1 << 1.;
    pf2 << 1.;
    p_set << 0.;
    p_max << pmax_mw;
    p0 << droop_p0_mw;
    k_deg << droop_k_mw_per_deg;
    grid.init_hvdc_lines(bus1, bus2, type1, type2,
                         zero, zero,           // no station losses
                         vreg1, vreg2, vm1_pu, vm2_pu,
                         zero, zero,           // q setpoints
                         q_min, q_max, q_min, q_max,
                         pf1, pf2, mode, p_set,
                         zero, zero,           // r_ohm, nominal_v: no dc line loss
                         droop_on, p0, k_deg, p_max, p_max);
    grid.tell_solver_need_reset();
}

// one voltage-mode SVC at `svc_bus` regulating `reg_bus` at `target_vm_pu`
// with the given (pu) slope
void add_voltage_svc(LSGrid & grid, int svc_bus, int reg_bus,
                     real_type target_vm_pu, real_type slope_pu)
{
    RealVect target_vm(1), q_set(1), slope(1), b_min(1), b_max(1);
    target_vm << target_vm_pu;
    q_set << 0.;
    slope << slope_pu;
    b_min << -100.;
    b_max << 100.;
    Eigen::VectorXi reg(1), bus(1);
    reg << reg_bus;
    bus << svc_bus;
    // SvcContainer::RegulationMode: VOLTAGE = 1
    grid.init_svcs({1}, target_vm, q_set, slope, b_min, b_max, reg, bus);
    grid.tell_solver_need_reset();
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
}

// ---------------------------------------------------------------------------
// driving an NRSystem the way NRAlgo::compute_pf does
// ---------------------------------------------------------------------------

// Owns local copies of everything an NRSystem keeps a pointer / Eigen::Ref
// into (Ybus and Sbus in particular: NRSystem::update_state caches
// Ybus.data() / Sbus.data() for the whole solve), so a probe can safely
// outlive the expression that built it.
struct SolverInputs
{
    Eigen::SparseMatrix<cplx_type> Ybus;
    CplxVect Sbus;
    CplxVect V;
    IntVect pv;
    IntVect pq;
    IntVect slack_ids;
    RealVect slack_weights;

    // the solver-labelled state of `grid` after a converged AC powerflow
    explicit SolverInputs(const LSGrid & grid):
        Ybus(grid.get_Ybus_solver()),
        Sbus(grid.get_Sbus_solver()),
        V(grid.get_V_solver()),
        pv(grid.get_pv_solver_numpy()),
        pq(grid.get_pq_solver_numpy()),
        slack_ids(grid.get_slack_ids_solver_numpy()),
        slack_weights(grid.get_slack_weights_solver())
    {}

    int nb_bus() const { return static_cast<int>(Ybus.rows()); }
};

// same call order as NRAlgo::compute_pf on a "need_rebuild" solve
template <class System>
void build_system(System & sys, const SolverInputs & in, const LSGrid * lsgrid_ptr)
{
    sys.update_state(lsgrid_ptr, in.Ybus, in.V, in.Sbus, in.slack_weights);
    sys.init_topology(in.slack_ids, in.slack_weights, in.pv, in.pq);
    sys.build_J_sparsity();
    sys.fill_internal_variables();
    sys.fill_J();
}

// solves `grid` with the requested NR flavour and returns its solver inputs
SolverInputs solved_inputs(LSGrid & grid, AlgorithmType algo)
{
    grid.change_algorithm(algo);
    const CplxVect V = grid.ac_pf(flat_start(grid), 30, 1e-10);
    REQUIRE(V.size() == static_cast<Eigen::Index>(grid.total_bus()));
    return SolverInputs(grid);
}

// ---------------------------------------------------------------------------
// invariant checks
// ---------------------------------------------------------------------------

// a bus-keyed map: one slot per solver bus, either -1 or a valid J index
void check_bus_map(const std::vector<int> & map, int nb_bus, int dim, const std::string & what)
{
    INFO("bus-keyed map: " << what);
    REQUIRE(static_cast<int>(map.size()) == nb_bus);
    for (int bus = 0; bus < nb_bus; ++bus) {
        INFO("bus " << bus << " -> " << map[bus] << " (dim " << dim << ")");
        CHECK(map[bus] >= -1);
        CHECK(map[bus] < dim);
    }
}

// a compact (bus, index) registration pair list, and its agreement with the
// bus-keyed map above: NRLedger documents "last registration wins" there, and
// indices are handed out by an increasing counter, so the map must hold the
// LARGEST index registered for that bus.
void check_pair_list(const std::vector<int> & buses,
                     const std::vector<int> & indices,
                     const std::vector<int> & map,
                     int nb_bus, int dim, const std::string & what)
{
    INFO("registration pair list: " << what);
    REQUIRE(buses.size() == indices.size());
    std::vector<int> largest(nb_bus, -1);
    for (size_t k = 0; k < buses.size(); ++k) {
        INFO("pair #" << k << ": bus " << buses[k] << " -> " << indices[k]);
        REQUIRE(buses[k] >= 0);
        REQUIRE(buses[k] < nb_bus);
        CHECK(indices[k] >= 0);
        CHECK(indices[k] < dim);
        if (indices[k] > largest[buses[k]]) largest[buses[k]] = indices[k];
    }
    for (int bus = 0; bus < nb_bus; ++bus) {
        INFO("bus " << bus);
        CHECK(map[bus] == largest[bus]);
    }
}

// A structurally empty row (or column) means a singular Jacobian: the equation
// has no unknown (resp. the unknown appears in no equation). It is also how an
// unresolved feature entry shows up -- the bordered VoltageControl rows and the
// MultiSlack slack_absorbed column carry NOTHING but feature entries, so they
// are empty here iff those entries were never declared / resolved.
void check_no_empty_row_or_col(const SpMatReal & J)
{
    const int dim = static_cast<int>(J.cols());
    std::vector<int> per_row(dim, 0);
    for (int col = 0; col < dim; ++col) {
        const int n_in_col = J.outerIndexPtr()[col + 1] - J.outerIndexPtr()[col];
        INFO("column " << col << " of " << dim);
        CHECK(n_in_col > 0);
        for (int p = J.outerIndexPtr()[col]; p < J.outerIndexPtr()[col + 1]; ++p) {
            const int row = J.innerIndexPtr()[p];
            REQUIRE(row >= 0);
            REQUIRE(row < dim);
            ++per_row[row];
        }
    }
    for (int row = 0; row < dim; ++row) {
        INFO("row " << row << " of " << dim);
        CHECK(per_row[row] > 0);
    }
}

template <class System>
void check_system_invariants(const System & sys, int nb_bus)
{
    const SpMatReal J = sys.J();
    const int dim = static_cast<int>(sys.total_state_variables());

    REQUIRE(dim > 0);
    CHECK(static_cast<int>(J.rows()) == dim);
    CHECK(static_cast<int>(J.cols()) == dim);
    CHECK(J.isCompressed());

    check_bus_map(sys.theta_to_J_col(), nb_bus, dim, "theta_to_J_col");
    check_bus_map(sys.vm_to_J_col(), nb_bus, dim, "vm_to_J_col");
    check_bus_map(sys.q_to_J_col(), nb_bus, dim, "q_to_J_col");
    check_bus_map(sys.p_to_J_row(), nb_bus, dim, "p_to_J_row");
    check_bus_map(sys.q_to_J_row(), nb_bus, dim, "q_to_J_row");

    check_pair_list(sys.p_buses(), sys.p_rows(), sys.p_to_J_row(), nb_bus, dim, "P equations");
    check_pair_list(sys.q_buses(), sys.q_rows(), sys.q_to_J_row(), nb_bus, dim, "Q equations");
    check_pair_list(sys.theta_buses(), sys.theta_cols(), sys.theta_to_J_col(), nb_bus, dim, "theta unknowns");
    check_pair_list(sys.vm_buses(), sys.vm_cols(), sys.vm_to_J_col(), nb_bus, dim, "vm unknowns");

    check_no_empty_row_or_col(J);

    // every filled value is a real number
    for (int k = 0; k < J.nonZeros(); ++k) {
        INFO("J nonzero #" << k);
        REQUIRE(std::isfinite(J.valuePtr()[k]));
    }

    // the residual spans the same space as the unknowns
    const RealVect F = sys.mismatch();
    CHECK(F.size() == dim);
    CHECK(F.allFinite());

    // the VoltageControl / MultiSlack result accessors agree with each other
    CHECK(sys.controller_q().size() == sys.controller_kind().size());
    CHECK(sys.controller_q().size() == sys.controller_elem_id().size());
    CHECK(sys.controller_q().size() == sys.controller_q_col().size());
    const IntVect q_cols = sys.controller_q_col();
    for (int c = 0; c < q_cols.size(); ++c) {
        INFO("controller " << c);
        CHECK(q_cols(c) >= 0);
        CHECK(q_cols(c) < dim);
    }
    CHECK(sys.slack_col() < dim);
}

// The augmented dimension, derived from the topology alone -- the ledger
// bookkeeping the components share. Base claims one theta unknown + one P
// equation per pv/pq bus, one vm unknown + one Q equation per pq bus and per
// slack bus that no local generator voltage-pins; VoltageControl adds one
// column and one row per controller (N columns vs 1 + (N-1) rows per group);
// MultiSlack, when present, adds one theta unknown per non-ref slack, one P
// equation per slack and the slack_absorbed column.
int expected_dim(const SolverInputs & in, const LSGrid & grid, int n_controllers, bool multi_slack)
{
    const int n_free_vm_slack = static_cast<int>(grid.get_free_vm_slack_solver_buses().size());
    int dim = static_cast<int>(in.pv.size()) + 2 * static_cast<int>(in.pq.size())
              + n_free_vm_slack + n_controllers;
    if (multi_slack) dim += static_cast<int>(in.slack_ids.size());
    return dim;
}

// the exact sparsity pattern, as a linear solver's symbolic phase sees it
std::vector<int> pattern_of(const SpMatReal & J)
{
    std::vector<int> res;
    res.push_back(static_cast<int>(J.rows()));
    res.push_back(static_cast<int>(J.cols()));
    res.push_back(static_cast<int>(J.nonZeros()));
    for (int col = 0; col <= static_cast<int>(J.cols()); ++col) res.push_back(J.outerIndexPtr()[col]);
    for (int k = 0; k < J.nonZeros(); ++k) res.push_back(J.innerIndexPtr()[k]);
    return res;
}

std::vector<real_type> values_of(const SpMatReal & J)
{
    return std::vector<real_type>(J.valuePtr(), J.valuePtr() + J.nonZeros());
}

// MW/deg -> pu/rad, the unit the Hvdc extension stores its droop slope in
real_type droop_k_pu_per_rad(real_type k_mw_per_deg, real_type sn_mva)
{
    return k_mw_per_deg * 180. / std::acos(real_type(-1.)) / sn_mva;
}

}  // anonymous namespace


// ===========================================================================
// the shared NRSystem machinery
// ===========================================================================

TEST_CASE("the NR system of a plain grid is square and consistent", "[NRSystem]")
{
    SECTION("multi-slack") {
        LSGrid grid = make_radial_grid();
        const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);
        MultiSlackNRSystem sys;
        build_system(sys, in, &grid);

        check_system_invariants(sys, in.nb_bus());
        CHECK(static_cast<int>(sys.total_state_variables()) == expected_dim(in, grid, 0, true));
        // the distributed-slack unknown is a real column of the augmented J
        CHECK(sys.slack_col() >= 0);
        // ... and it is the LAST column claimed (MultiSlack registers it last,
        // and neither VoltageControl nor Hvdc claims anything on this grid)
        CHECK(sys.slack_col() == static_cast<int>(sys.total_state_variables()) - 1);
        // no voltage controller here
        CHECK(sys.controller_q().size() == 0);
    }

    SECTION("single-slack") {
        LSGrid grid = make_radial_grid();
        const SolverInputs in = solved_inputs(grid, AlgorithmType::NRSing_SparseLU);
        SingleSlackNRSystem sys;
        build_system(sys, in, &grid);

        check_system_invariants(sys, in.nb_bus());
        CHECK(static_cast<int>(sys.total_state_variables()) == expected_dim(in, grid, 0, false));
        // no MultiSlack extension in this instantiation
        CHECK(sys.slack_col() == -1);
        CHECK(sys.slack_absorbed() == 0.);
        // the slack bus owns neither an equation nor an unknown here
        const int slack = in.slack_ids(0);
        CHECK(sys.p_to_J_row()[slack] == -1);
        CHECK(sys.q_to_J_row()[slack] == -1);
        CHECK(sys.theta_to_J_col()[slack] == -1);
        CHECK(sys.vm_to_J_col()[slack] == -1);
        // the slack generator pins its own bus, so no free-Vm slack bus
        CHECK(grid.get_free_vm_slack_solver_buses().empty());
    }
}

TEST_CASE("the single-slack residual vanishes at the converged voltages", "[NRSystem]")
{
    // The whole residual assembly (Ybus product, Sbus, the ledger's per-bus
    // row scatter) evaluated at the solution the algorithm itself reached: it
    // must reproduce the convergence. Single-slack only -- the MultiSlack /
    // VoltageControl extensions restart their own state (slack_absorbed,
    // controller Q) from a per-solve initial guess in update_state, so their
    // rows are legitimately non-zero at V_converged.
    LSGrid grid = make_radial_grid(4);
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NRSing_SparseLU);
    SingleSlackNRSystem sys;
    build_system(sys, in, &grid);

    check_system_invariants(sys, in.nb_bus());
    CHECK(sys.mismatch().lpNorm<Eigen::Infinity>() < 1e-9);
    // the cached V/Va/Vm views agree with the input
    CHECK((sys.V() - in.V).norm() < 1e-12);
    CHECK((sys.Vm() - in.V.array().abs().matrix()).norm() < 1e-12);
}

TEST_CASE("clear_jacobian resets the system for reuse", "[NRSystem]")
{
    LSGrid grid = make_radial_grid();
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);

    MultiSlackNRSystem sys;
    build_system(sys, in, &grid);
    const SpMatReal J_first = sys.J();
    const std::vector<int> pattern_first = pattern_of(J_first);
    const std::vector<real_type> values_first = values_of(J_first);
    const RealVect F_first = sys.mismatch();

    sys.clear_jacobian();
    CHECK(sys.total_state_variables() == 0);
    CHECK(sys.J().nonZeros() == 0);
    CHECK(sys.slack_col() == -1);

    // rebuilding from scratch must give back exactly the same system
    build_system(sys, in, &grid);
    const SpMatReal J_second = sys.J();
    CHECK(pattern_of(J_second) == pattern_first);
    CHECK(values_of(J_second) == values_first);
    CHECK((sys.mismatch() - F_first).norm() == 0.);
    check_system_invariants(sys, in.nb_bus());
}

TEST_CASE("a null grid pointer is a no-op for the grid-backed extensions", "[NRSystem]")
{
    // A solver used standalone (the `newtonpf` entry point) never sees an
    // LSGrid: BaseAlgo::lsgrid_ptr_ stays nullptr and is forwarded as such.
    // On a grid with no hvdc droop, no voltage controller and a locally
    // pinned slack there is nothing to pull, so the two builds must agree
    // bit for bit.
    LSGrid grid = make_radial_grid();
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);
    REQUIRE(grid.get_free_vm_slack_solver_buses().empty());

    MultiSlackNRSystem with_grid;
    build_system(with_grid, in, &grid);
    MultiSlackNRSystem without_grid;
    build_system(without_grid, in, nullptr);

    check_system_invariants(without_grid, in.nb_bus());
    const SpMatReal J_with = with_grid.J();
    const SpMatReal J_without = without_grid.J();
    CHECK(pattern_of(J_without) == pattern_of(J_with));
    CHECK(values_of(J_without) == values_of(J_with));
    CHECK((without_grid.mismatch() - with_grid.mismatch()).norm() == 0.);
}

TEST_CASE("bus masking replaces the masked rows by the identity", "[NRSystem][mask]")
{
    // ContingencyAnalysis' "handle disconnected grid" mode. The documented
    // contract is that this is a pure VALUE-level edit: the J pattern (hence
    // the symbolic factorization) is untouched.
    LSGrid grid = make_radial_grid(4);
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);

    MultiSlackNRSystem plain;
    build_system(plain, in, &grid);
    const std::vector<int> pattern_plain = pattern_of(SpMatReal(plain.J()));

    const int masked_bus = in.pq(in.pq.size() - 1);  // never the reference slack
    MultiSlackNRSystem masked;
    masked.set_masked_buses(std::vector<int>{masked_bus});
    build_system(masked, in, &grid);

    check_system_invariants(masked, in.nb_bus());
    const SpMatReal J = masked.J();
    CHECK(pattern_of(J) == pattern_plain);

    const int p_row = masked.p_to_J_row()[masked_bus];
    const int q_row = masked.q_to_J_row()[masked_bus];
    const int theta_col = masked.theta_to_J_col()[masked_bus];
    const int vm_col = masked.vm_to_J_col()[masked_bus];
    REQUIRE(p_row >= 0);
    REQUIRE(q_row >= 0);
    REQUIRE(theta_col >= 0);
    REQUIRE(vm_col >= 0);

    // every stored coefficient of the two masked rows is 0, except the single
    // diagonal-like entry (p_row, theta_col) / (q_row, vm_col), which is 1
    int nb_ones = 0;
    for (int col = 0; col < static_cast<int>(J.cols()); ++col) {
        for (int p = J.outerIndexPtr()[col]; p < J.outerIndexPtr()[col + 1]; ++p) {
            const int row = J.innerIndexPtr()[p];
            if ((row != p_row) && (row != q_row)) continue;
            const bool is_one = ((row == p_row) && (col == theta_col)) ||
                                ((row == q_row) && (col == vm_col));
            INFO("masked entry (" << row << ", " << col << ")");
            if (is_one) {
                CHECK(J.valuePtr()[p] == 1.);
                ++nb_ones;
            } else {
                CHECK(J.valuePtr()[p] == 0.);
            }
        }
    }
    CHECK(nb_ones == 2);

    // ... and the matching residual entries are forced to 0, so the step on
    // that bus is exactly 0
    const RealVect F = masked.mismatch();
    CHECK(F(p_row) == 0.);
    CHECK(F(q_row) == 0.);

    // masking is part of the state clear_jacobian resets
    masked.clear_jacobian();
    build_system(masked, in, &grid);
    const SpMatReal J_after = masked.J();
    CHECK(values_of(J_after) == values_of(SpMatReal(plain.J())));
}


// ===========================================================================
// MultiSlack
// ===========================================================================

TEST_CASE("MultiSlack borders the system with one column and one row per slack", "[NRSystem][slack]")
{
    LSGrid grid = make_radial_skeleton(4);
    add_generators(grid, {0, 1}, {1.02, 1.02}, {2000., 2000.});
    grid.add_gen_slackbus(0, 0.7);
    grid.add_gen_slackbus(1, 0.3);
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);
    REQUIRE(in.slack_ids.size() == 2);

    MultiSlackNRSystem sys;
    build_system(sys, in, &grid);
    check_system_invariants(sys, in.nb_bus());
    CHECK(static_cast<int>(sys.total_state_variables()) == expected_dim(in, grid, 0, true));

    // every slack bus owns a P equation; the non-reference ones also own a
    // theta unknown, the reference one does not (its angle is the reference)
    const int ref_slack = in.slack_ids(0);
    CHECK(sys.theta_to_J_col()[ref_slack] == -1);
    for (int k = 1; k < in.slack_ids.size(); ++k) {
        INFO("non-reference slack bus " << in.slack_ids(k));
        CHECK(sys.theta_to_J_col()[in.slack_ids(k)] >= 0);
    }

    // the slack_absorbed column carries exactly one feature entry per slack
    // bus, at that bus' P row, worth its participation weight. Nothing else
    // writes into that column (it is bus-less, so no dS entry lands there),
    // so the value is the weight itself.
    const SpMatReal J = sys.J();
    const int slack_col = sys.slack_col();
    REQUIRE(slack_col >= 0);
    int nb_in_col = 0;
    for (int p = J.outerIndexPtr()[slack_col]; p < J.outerIndexPtr()[slack_col + 1]; ++p) {
        ++nb_in_col;
        const int row = J.innerIndexPtr()[p];
        // find which slack bus owns this P row
        int owner = -1;
        for (int k = 0; k < in.slack_ids.size(); ++k)
            if (sys.p_to_J_row()[in.slack_ids(k)] == row) owner = in.slack_ids(k);
        INFO("slack column entry at row " << row);
        REQUIRE(owner >= 0);
        CHECK(J.valuePtr()[p] == Approx(in.slack_weights(owner)).epsilon(1e-12));
    }
    CHECK(nb_in_col == static_cast<int>(in.slack_ids.size()));

    // the converged distributed-slack state is the initial guess (Sbus sum)
    // only until apply_step runs; here it is exactly the per-solve initial
    // value NRAlgo starts from
    CHECK(sys.slack_absorbed() == Approx(std::real(in.Sbus.sum())).epsilon(1e-12));
}

TEST_CASE("a PQ distributed-slack participant keeps a free Vm unknown", "[NRSystem][slack]")
{
    // A slack bus that no LOCAL generator voltage-pins needs a free Vm
    // unknown + Q equation, exactly like an ordinary PQ bus. Base owns that
    // (not MultiSlack), so it must happen in the single-slack instantiation
    // too. Here the generator at bus 1 regulates bus 2 remotely, so it does
    // NOT pin bus 1, which is nonetheless a slack participant.
    LSGrid grid = make_radial_skeleton(4);
    add_generators(grid, {0, 1}, {1.02, 1.03}, {2000., 2000.});
    grid.add_gen_slackbus(0, 0.5);
    grid.add_gen_slackbus(1, 0.5);
    grid.set_gen_regulated_bus(1, 2);
    grid.tell_solver_need_reset();

    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);
    const std::set<int> free_vm = grid.get_free_vm_slack_solver_buses();
    REQUIRE(free_vm.size() == 1);
    const int free_bus = *free_vm.begin();

    MultiSlackNRSystem sys;
    build_system(sys, in, &grid);
    check_system_invariants(sys, in.nb_bus());
    CHECK(static_cast<int>(sys.total_state_variables()) == expected_dim(in, grid, 1, true));

    // that slack bus has the full set: P equation + theta (or reference),
    // AND the extra Vm unknown / Q equation
    CHECK(sys.vm_to_J_col()[free_bus] >= 0);
    CHECK(sys.q_to_J_row()[free_bus] >= 0);
    CHECK(sys.p_to_J_row()[free_bus] >= 0);
}


// ===========================================================================
// VoltageControl
// ===========================================================================

TEST_CASE("VoltageControl borders the system squarely, group by group", "[NRSystem][vctrl]")
{
    // two generators (bus 1 and bus 2) remotely regulating bus 3 at the same
    // setpoint: ONE control group of two controllers, i.e. 2 reactive-injection
    // columns against 1 voltage row + 1 sharing row.
    LSGrid grid = make_radial_skeleton(4);
    add_generators(grid, {0, 1, 2}, {1.02, 1.04, 1.04}, {2000., 2000., 500.});
    grid.add_gen_slackbus(0, 1.);
    grid.set_gen_regulated_bus(1, 3);
    grid.set_gen_regulated_bus(2, 3);
    grid.tell_solver_need_reset();

    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);
    MultiSlackNRSystem sys;
    build_system(sys, in, &grid);

    check_system_invariants(sys, in.nb_bus());
    REQUIRE(sys.controller_q().size() == 2);
    CHECK(static_cast<int>(sys.total_state_variables()) == expected_dim(in, grid, 2, true));

    // both controllers are generators, reported in registration order
    const IntVect kind = sys.controller_kind();
    const IntVect elem = sys.controller_elem_id();
    for (int c = 0; c < kind.size(); ++c) {
        INFO("controller " << c);
        CHECK(kind(c) == static_cast<int>(ls2g::VoltageControlSolverData::GEN));
    }
    CHECK(((elem(0) == 1 && elem(1) == 2) || (elem(0) == 2 && elem(1) == 1)));

    // Two controllers on two DIFFERENT buses get two distinct columns. The
    // bus-keyed q_to_J_col map only keeps the last registration per bus, so
    // controller_q_col() is the authoritative one; here both agree.
    const IntVect q_cols = sys.controller_q_col();
    CHECK(q_cols(0) != q_cols(1));

    // each controller's own Q unknown appears in its own bus' Q equation with
    // the -1 coefficient of the bordered formulation (nothing else lands in
    // that column at that row: the column is bus-less for the dS pass)
    const SpMatReal J = sys.J();
    for (int c = 0; c < q_cols.size(); ++c) {
        const int col = q_cols(c);
        int nb_minus_one = 0;
        for (int p = J.outerIndexPtr()[col]; p < J.outerIndexPtr()[col + 1]; ++p) {
            if (J.valuePtr()[p] == Approx(-1.).epsilon(1e-12)) ++nb_minus_one;
        }
        INFO("controller " << c << ", J column " << col);
        CHECK(nb_minus_one >= 1);
    }

    // the physical outcome of the sharing rows: reactive output proportional
    // to the sharing key (qmax - qmin), which is 2000 MVAr vs 500 MVAr here
    const RealVect q = grid.get_algo().get_controller_q();
    REQUIRE(q.size() == 2);
    const IntVect elem_ids = grid.get_algo().get_controller_elem_id();
    real_type q_gen1 = 0., q_gen2 = 0.;
    for (int c = 0; c < q.size(); ++c) {
        if (elem_ids(c) == 1) q_gen1 = q(c);
        if (elem_ids(c) == 2) q_gen2 = q(c);
    }
    CHECK(std::abs(q_gen1) > 1e-6);
    CHECK(q_gen1 / 2000. == Approx(q_gen2 / 500.).epsilon(1e-8));
    // and the regulated bus sits at the setpoint
    CHECK(std::abs(grid.get_V()(3)) == Approx(1.04).epsilon(1e-8));
}

TEST_CASE("a sloped SVC trades voltage against reactive output", "[NRSystem][vctrl][svc]")
{
    // A voltage-mode SVC is ALWAYS a VoltageControl controller (even local and
    // non-sloped). With a slope the bordered voltage row is
    //     Vm(reg) + s.Q - v_set = 0
    // so the regulated bus settles BELOW the setpoint when the SVC injects.
    const real_type v_set = 1.05;
    const real_type slope_pu = 0.02;

    LSGrid grid = make_radial_grid();
    add_voltage_svc(grid, /*svc_bus=*/2, /*reg_bus=*/2, v_set, slope_pu);
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);

    MultiSlackNRSystem sys;
    build_system(sys, in, &grid);
    check_system_invariants(sys, in.nb_bus());
    REQUIRE(sys.controller_q().size() == 1);
    CHECK(sys.controller_kind()(0) == static_cast<int>(ls2g::VoltageControlSolverData::SVC));
    CHECK(static_cast<int>(sys.total_state_variables()) == expected_dim(in, grid, 1, true));

    // the slope law, on the converged solution of the grid itself
    const RealVect q = grid.get_algo().get_controller_q();
    REQUIRE(q.size() == 1);
    CHECK(q(0) > 0.);  // it has to raise the bus, so it injects
    CHECK(std::abs(grid.get_V()(2)) == Approx(v_set - slope_pu * q(0)).epsilon(1e-8));
    CHECK(std::abs(grid.get_V()(2)) < v_set);

    // The slope entry is declared for EVERY SVC controller, slope or not, so
    // that changing the slope never changes the J pattern (symbolic reuse).
    LSGrid flat = make_radial_grid();
    add_voltage_svc(flat, 2, 2, v_set, /*slope_pu=*/0.);
    const SolverInputs in_flat = solved_inputs(flat, AlgorithmType::NR_SparseLU);
    MultiSlackNRSystem sys_flat;
    build_system(sys_flat, in_flat, &flat);
    CHECK(pattern_of(SpMatReal(sys_flat.J())) == pattern_of(SpMatReal(sys.J())));
    // a non-sloped SVC does hold its bus exactly at the setpoint
    CHECK(std::abs(flat.get_V()(2)) == Approx(v_set).epsilon(1e-8));
}


// ===========================================================================
// Hvdc (angle droop)
// ===========================================================================

TEST_CASE("the hvdc droop claims no row or column of its own", "[NRSystem][hvdc]")
{
    // The Hvdc extension only writes into rows / columns other components
    // already own, so adding a droop line must not change the dimension.
    LSGrid plain = make_radial_grid();
    const SolverInputs in_plain = solved_inputs(plain, AlgorithmType::NR_SparseLU);
    MultiSlackNRSystem sys_plain;
    build_system(sys_plain, in_plain, &plain);

    LSGrid grid = make_radial_grid();
    add_droop_hvdc(grid, /*bus_1=*/1, /*bus_2=*/2, /*p0=*/10., /*k=*/2.);
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);
    MultiSlackNRSystem sys;
    build_system(sys, in, &grid);

    check_system_invariants(sys, in.nb_bus());
    CHECK(sys.total_state_variables() == sys_plain.total_state_variables());
    CHECK(static_cast<int>(sys.total_state_variables()) == expected_dim(in, grid, 0, true));
    // buses 1 and 2 are adjacent on the feeder, so all four droop entries land
    // on positions the dS pass already created: the pattern is unchanged too
    CHECK(pattern_of(SpMatReal(sys.J())) == pattern_of(SpMatReal(sys_plain.J())));
}

TEST_CASE("the single-slack residual absorbs the droop flows", "[NRSystem][hvdc]")
{
    LSGrid grid = make_radial_grid();
    add_droop_hvdc(grid, 1, 2, 10., 2.);
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NRSing_SparseLU);

    SingleSlackNRSystem sys;
    build_system(sys, in, &grid);
    check_system_invariants(sys, in.nb_bus());
    // the theta-dependent droop injections are NOT in Sbus: only
    // Hvdc::adjust_mismatch puts them back, so a zero residual at the
    // converged voltages is a check of that path
    CHECK(sys.mismatch().lpNorm<Eigen::Infinity>() < 1e-9);
}

TEST_CASE("a status_droop flip keeps the J sparsity pattern", "[NRSystem][hvdc]")
{
    // The four dP/dtheta entries are declared whatever the regime, so a
    // saturation flip between two solves changes VALUES only -- the symbolic
    // factorization of the linear solver stays valid.
    const real_type p0 = 10., k_mw_per_deg = 2.;
    LSGrid grid = make_radial_grid();
    add_droop_hvdc(grid, 1, 2, p0, k_mw_per_deg);
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);
    REQUIRE(grid.get_status_droop_hvdc(0) == 0);  // linear regime

    MultiSlackNRSystem linear;
    build_system(linear, in, &grid);
    const SpMatReal J_linear = linear.J();

    // saturate the line 1 -> 2 (a pure input, as documented on HvdcDroopSolverData)
    grid.set_status_droop_hvdc(0, 1);
    MultiSlackNRSystem saturated;
    build_system(saturated, in, &grid);
    const SpMatReal J_saturated = saturated.J();

    check_system_invariants(saturated, in.nb_bus());
    CHECK(pattern_of(J_saturated) == pattern_of(J_linear));
    // ... but the slopes are gone: a saturated line is a constant injection
    CHECK(values_of(J_saturated) != values_of(J_linear));

    // and the difference is exactly the droop slopes. Stations are lossless
    // and so is the dc line, so loss_mult == 1 and dp1 = k, dp2 = -k.
    const real_type k_pu = droop_k_pu_per_rad(k_mw_per_deg, 100.);
    const int p1 = linear.p_to_J_row()[1], p2 = linear.p_to_J_row()[2];
    const int t1 = linear.theta_to_J_col()[1], t2 = linear.theta_to_J_col()[2];
    REQUIRE(p1 >= 0);
    REQUIRE(p2 >= 0);
    REQUIRE(t1 >= 0);
    REQUIRE(t2 >= 0);
    CHECK(J_linear.coeff(p1, t1) - J_saturated.coeff(p1, t1) == Approx(k_pu).epsilon(1e-10));
    CHECK(J_linear.coeff(p1, t2) - J_saturated.coeff(p1, t2) == Approx(-k_pu).epsilon(1e-10));
    CHECK(J_linear.coeff(p2, t1) - J_saturated.coeff(p2, t1) == Approx(-k_pu).epsilon(1e-10));
    CHECK(J_linear.coeff(p2, t2) - J_saturated.coeff(p2, t2) == Approx(k_pu).epsilon(1e-10));

    // the residual changes too: a saturated line carries pmax, not the droop law
    CHECK((saturated.mismatch() - linear.mismatch()).lpNorm<Eigen::Infinity>() > 1e-6);
}

TEST_CASE("a droop end on the reference slack drops its dropped entries", "[NRSystem][hvdc]")
{
    // Side 1 of the hvdc sits on the reference slack bus, which owns no theta
    // unknown: the (P row, theta column) entries that would need it simply are
    // not declared. Buses 0 and 2 are NOT adjacent on the feeder, so the
    // surviving (P row of bus 0, theta column of bus 2) entry is created by
    // the droop alone -- its value is the droop slope and nothing else.
    const real_type p0 = 10., k_mw_per_deg = 2.;
    LSGrid grid = make_radial_grid();
    add_droop_hvdc(grid, /*bus_1=*/0, /*bus_2=*/2, p0, k_mw_per_deg);
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);
    const int slack = in.slack_ids(0);
    REQUIRE(slack == 0);

    MultiSlackNRSystem sys;
    build_system(sys, in, &grid);
    check_system_invariants(sys, in.nb_bus());

    const int p_slack = sys.p_to_J_row()[0];
    const int t_slack = sys.theta_to_J_col()[0];
    const int t_2 = sys.theta_to_J_col()[2];
    // the reference slack owns a P equation (MultiSlack) but no theta unknown
    REQUIRE(p_slack >= 0);
    CHECK(t_slack == -1);
    REQUIRE(t_2 >= 0);

    const SpMatReal J = sys.J();
    const real_type k_pu = droop_k_pu_per_rad(k_mw_per_deg, 100.);
    CHECK(J.coeff(p_slack, t_2) == Approx(-k_pu).epsilon(1e-10));

    // without the grid the extension pulls no droop data at all, so that
    // entry is not even declared: the pattern is strictly smaller
    MultiSlackNRSystem no_grid;
    build_system(no_grid, in, nullptr);
    check_system_invariants(no_grid, in.nb_bus());
    CHECK(no_grid.total_state_variables() == sys.total_state_variables());
    CHECK(no_grid.J().nonZeros() < J.nonZeros());
    CHECK(SpMatReal(no_grid.J()).coeff(p_slack, t_2) == 0.);
}

TEST_CASE("a droop end on the slack of a single-slack system drops three entries", "[NRSystem][hvdc]")
{
    // Same grid, single-slack instantiation: the slack bus owns NEITHER a P
    // equation NOR a theta unknown, so only the (P row of bus 2, theta column
    // of bus 2) entry of the four survives -- and it lands on a position the
    // dS pass already owns.
    LSGrid grid = make_radial_grid();
    add_droop_hvdc(grid, 0, 2, 10., 2.);
    const SolverInputs in = solved_inputs(grid, AlgorithmType::NRSing_SparseLU);

    SingleSlackNRSystem sys;
    build_system(sys, in, &grid);
    check_system_invariants(sys, in.nb_bus());
    CHECK(sys.p_to_J_row()[0] == -1);
    CHECK(sys.theta_to_J_col()[0] == -1);

    LSGrid plain = make_radial_grid();
    const SolverInputs in_plain = solved_inputs(plain, AlgorithmType::NRSing_SparseLU);
    SingleSlackNRSystem sys_plain;
    build_system(sys_plain, in_plain, &plain);
    CHECK(pattern_of(SpMatReal(sys.J())) == pattern_of(SpMatReal(sys_plain.J())));

    // the whole droop law still reaches the residual through both end buses
    CHECK(sys.mismatch().lpNorm<Eigen::Infinity>() < 1e-9);
}


// ===========================================================================
// all three extensions at once
// ===========================================================================

TEST_CASE("a grid exercising every NR extension at once stays consistent", "[NRSystem]")
{
    // distributed slack (bus 0 + bus 1) + a remote-regulating generator
    // (bus 2 -> bus 4) + an angle-droop hvdc (bus 3 <-> bus 4): MultiSlack,
    // VoltageControl and Hvdc all register into the same ledger, in that order.
    LSGrid grid = make_radial_skeleton(5);
    add_generators(grid, {0, 1, 2}, {1.02, 1.02, 1.03}, {2000., 2000., 800.});
    grid.add_gen_slackbus(0, 0.6);
    grid.add_gen_slackbus(1, 0.4);
    grid.set_gen_regulated_bus(2, 4);
    add_droop_hvdc(grid, /*bus_1=*/3, /*bus_2=*/4, /*p0=*/5., /*k=*/1.5);
    grid.tell_solver_need_reset();

    const SolverInputs in = solved_inputs(grid, AlgorithmType::NR_SparseLU);
    MultiSlackNRSystem sys;
    build_system(sys, in, &grid);

    check_system_invariants(sys, in.nb_bus());
    REQUIRE(sys.controller_q().size() == 1);
    CHECK(static_cast<int>(sys.total_state_variables()) == expected_dim(in, grid, 1, true));
    CHECK(sys.slack_col() >= 0);
    // the slack column is claimed by MultiSlack, BEFORE VoltageControl's
    // reactive-injection columns
    CHECK(sys.slack_col() < sys.controller_q_col()(0));
    // the regulated bus is held at the setpoint of its controller
    CHECK(std::abs(grid.get_V()(4)) == Approx(1.03).epsilon(1e-8));

    // rebuilding the very same system twice in a row is idempotent
    const std::vector<int> pattern = pattern_of(SpMatReal(sys.J()));
    const std::vector<real_type> values = values_of(SpMatReal(sys.J()));
    build_system(sys, in, &grid);
    CHECK(pattern_of(SpMatReal(sys.J())) == pattern);
    CHECK(values_of(SpMatReal(sys.J())) == values);
}
