// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Tests for LSGrid::check_grid() (the whole-grid consistency validator) and for
// the solver-output sanity check at the algorithm->grid seam. Everything is
// driven from C++, no python / grid2op / external grid file.
//
// The concern: release wheels are built -O3 -DNDEBUG, so an out-of-range index
// (bus / substation / topology-vector) restored from a well-formed but poisoned
// pickle or binary file would cause an out-of-bounds read/write on the next
// powerflow rather than a clean error. check_grid() (called at the end of
// set_state) must turn every such case into a std::out_of_range / runtime_error.
// These cases also run under valgrind in the cpp_unit_tests workflow, which is
// what proves no out-of-bounds access happens on the poisoned inputs.

#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "AlgorithmRegistry.hpp"
#include "powerflow_algorithm/BaseAlgo.hpp"

using ls2g::LSGrid;
using ls2g::CplxVect;
using ls2g::RealVect;
using ls2g::IntVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

// Three 138 kV buses (3 substations, 1 busbar each), two lines 0-1 and 1-2,
// one 50 MW / 10 MVar load at bus 2 and one slack generator (1.02 pu) at bus 0.
// A valid, internally-consistent grid: check_grid() must accept it.
LSGrid make_valid_grid()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    RealVect bus_vn_kv(3);
    bus_vn_kv << 138., 138., 138.;
    grid.init_bus(3, 1, bus_vn_kv, 0, 0);

    RealVect branch_r(2), branch_x(2);
    branch_r << 0.01, 0.01;
    branch_x << 0.1, 0.1;
    const CplxVect branch_h = CplxVect::Zero(2);
    Eigen::VectorXi from_id(2), to_id(2);
    from_id << 0, 1;
    to_id << 1, 2;
    grid.init_powerlines(branch_r, branch_x, branch_h, from_id, to_id);

    RealVect load_p(1), load_q(1);
    load_p << 50.;
    load_q << 10.;
    Eigen::VectorXi load_bus(1);
    load_bus << 2;
    grid.init_loads(load_p, load_q, load_bus);

    RealVect gen_p(1), gen_v(1), gen_min_q(1), gen_max_q(1);
    gen_p << 0.;
    gen_v << 1.02;
    gen_min_q << -1000.;
    gen_max_q << 1000.;
    Eigen::VectorXi gen_bus(1);
    gen_bus << 0;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    return grid;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<int>(grid.total_bus()), cplx_type(1., 0.));
}

// --- accessors into the (deeply nested) StateRes tuple -----------------------
// LOAD_ID -> LoadContainer::StateRes = tuple<OneSideContainer_PQ::StateRes>
//         -> [0] OneSideContainer_PQ::StateRes = tuple<OSC::StateRes, p_mw, q_mvar>
//         -> [0] OneSideContainer::StateRes    = tuple<names, bus_id, status,
//                                                       has_subid, subid,
//                                                       has_topo, pos_topo_vect>
std::vector<int>& load_bus_id(LSGrid::StateRes & st)
{
    return std::get<1>(std::get<0>(std::get<0>(std::get<LSGrid::LOAD_ID>(st))));
}
// GEN_ID -> GeneratorContainer::StateRes = tuple<OSC_PQ::StateRes, bool,
//   voltage_reg, vm_pu, min_q, max_q, [6] slack_bus, [7] slack_weight,
//   [8] regulated_bus_id>
std::vector<real_type>& gen_slack_weight(LSGrid::StateRes & st)
{
    return std::get<7>(std::get<LSGrid::GEN_ID>(st));
}
// AC_ALGO_CONFIG_ID -> tuple<int_params, real_params> (an AlgoConfig)
// int_params  = {scaling_policy, refactor_policy, ls_max_iter, refactor_every_n}
// real_params = {max_dVa, max_dVm, ls_c, ls_rho, iw_mu_min, iw_mu_max}
std::vector<int>& ac_algo_int_params(LSGrid::StateRes & st)
{
    return std::get<0>(std::get<LSGrid::AC_ALGO_CONFIG_ID>(st));
}
std::vector<double>& ac_algo_real_params(LSGrid::StateRes & st)
{
    return std::get<1>(std::get<LSGrid::AC_ALGO_CONFIG_ID>(st));
}
std::vector<int>& gen_regulated_bus(LSGrid::StateRes & st)
{
    return std::get<8>(std::get<LSGrid::GEN_ID>(st));
}

// A well-behaved AC solver (returns the input voltages, claims convergence).
// Registered under a non-built-in name it stands in for a plugin solver, whose
// AlgorithmType is the catch-all AlgorithmType::Custom.
class PluginLikeAlgo : public ls2g::BaseAlgo {
public:
    PluginLikeAlgo() : ls2g::BaseAlgo(/*is_ac=*/true) {}
    bool compute_pf(const ls2g::EigenRefConstCplxSpMat & /*Ybus*/,
                    const Eigen::Ref<const CplxVect> & V,
                    const Eigen::Ref<const CplxVect> & /*Sbus*/,
                    const Eigen::Ref<const ls2g::IntVect> & /*slack_ids*/,
                    const Eigen::Ref<const RealVect> & /*slack_weights*/,
                    const Eigen::Ref<const ls2g::IntVect> & /*pv*/,
                    const Eigen::Ref<const ls2g::IntVect> & /*pq*/,
                    int /*max_iter*/, real_type /*tol*/) override
    {
        V_ = V;
        Va_ = V.array().arg();
        Vm_ = V.array().abs();
        n_ = static_cast<int>(V.size());
        nr_iter_ = 1;
        err_ = ls2g::ErrorType::NoError;
        return true;
    }
};

// The DC counterpart of PluginLikeAlgo: IS_AC == false, so LSGrid routes it to
// the DC slot (see LSGrid::change_algorithm, which dispatches on IS_AC).
class DcPluginLikeAlgo : public ls2g::BaseAlgo {
public:
    DcPluginLikeAlgo() : ls2g::BaseAlgo(/*is_ac=*/false) {}
    bool compute_pf_dc(const ls2g::EigenRefConstRealSpMat & /*Bbus*/,
                       const Eigen::Ref<const CplxVect> & V,
                       const Eigen::Ref<const RealVect> & /*Pbus*/,
                       const Eigen::Ref<const ls2g::IntVect> & /*slack_ids*/,
                       const Eigen::Ref<const RealVect> & /*slack_weights*/,
                       const Eigen::Ref<const ls2g::IntVect> & /*pv*/,
                       const Eigen::Ref<const ls2g::IntVect> & /*pq*/) override
    {
        V_ = V;
        Va_ = V.array().arg();
        Vm_ = V.array().abs();
        n_ = static_cast<int>(V.size());
        nr_iter_ = 1;
        err_ = ls2g::ErrorType::NoError;
        return true;
    }
};

// A BaseAlgo that claims convergence but returns V/Va/Vm one entry too long.
class WrongSizeAlgo : public ls2g::BaseAlgo {
public:
    WrongSizeAlgo() : ls2g::BaseAlgo(/*is_ac=*/true) {}
    bool compute_pf(const ls2g::EigenRefConstCplxSpMat & /*Ybus*/,
                    const Eigen::Ref<const CplxVect> & V,
                    const Eigen::Ref<const CplxVect> & /*Sbus*/,
                    const Eigen::Ref<const ls2g::IntVect> & /*slack_ids*/,
                    const Eigen::Ref<const RealVect> & /*slack_weights*/,
                    const Eigen::Ref<const ls2g::IntVect> & /*pv*/,
                    const Eigen::Ref<const ls2g::IntVect> & /*pq*/,
                    int /*max_iter*/, real_type /*tol*/) override
    {
        const int bad = static_cast<int>(V.size()) + 1;
        V_ = CplxVect::Constant(bad, cplx_type(1., 0.));
        Va_ = RealVect::Zero(bad);
        Vm_ = RealVect::Constant(bad, 1.);
        n_ = bad;
        nr_iter_ = 1;
        err_ = ls2g::ErrorType::NoError;  // wrongly claims success
        return true;
    }
};

// A BaseAlgo that claims convergence but injects a NaN into the result.
class NaNAlgo : public ls2g::BaseAlgo {
public:
    NaNAlgo() : ls2g::BaseAlgo(/*is_ac=*/true) {}
    bool compute_pf(const ls2g::EigenRefConstCplxSpMat & /*Ybus*/,
                    const Eigen::Ref<const CplxVect> & V,
                    const Eigen::Ref<const CplxVect> & /*Sbus*/,
                    const Eigen::Ref<const ls2g::IntVect> & /*slack_ids*/,
                    const Eigen::Ref<const RealVect> & /*slack_weights*/,
                    const Eigen::Ref<const ls2g::IntVect> & /*pv*/,
                    const Eigen::Ref<const ls2g::IntVect> & /*pq*/,
                    int /*max_iter*/, real_type /*tol*/) override
    {
        V_ = V;
        if (V_.size() > 0) V_(0) = cplx_type(std::numeric_limits<real_type>::quiet_NaN(), 0.);
        Va_ = V_.array().arg();
        Vm_ = V_.array().abs();
        n_ = static_cast<int>(V_.size());
        nr_iter_ = 1;
        err_ = ls2g::ErrorType::NoError;  // wrongly claims success
        return true;
    }
};


} // namespace

TEST_CASE("check_grid accepts a valid grid and its state round-trip", "[check_grid]")
{
    LSGrid grid = make_valid_grid();
    CHECK_NOTHROW(grid.check_grid());

    LSGrid::StateRes st = grid.get_state();
    LSGrid restored;
    CHECK_NOTHROW(restored.set_state(st));  // set_state runs check_grid internally
}

TEST_CASE("check_grid rejects an out-of-range bus id on load (set_state)", "[check_grid]")
{
    LSGrid grid = make_valid_grid();
    LSGrid::StateRes st = grid.get_state();
    REQUIRE(load_bus_id(st).size() == 1);
    load_bus_id(st)[0] = 9999;  // >> nb_bus (3)

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::out_of_range);
}

TEST_CASE("check_grid rejects a negative (non-sentinel) bus id", "[check_grid]")
{
    LSGrid grid = make_valid_grid();
    LSGrid::StateRes st = grid.get_state();
    load_bus_id(st)[0] = -5;  // -1 would mean "disconnected"; -5 is invalid

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::out_of_range);
}

TEST_CASE("check_grid rejects a connected element sitting on the -1 (disconnected) bus", "[check_grid]")
{
    LSGrid grid = make_valid_grid();
    LSGrid::StateRes st = grid.get_state();
    load_bus_id(st)[0] = -1;  // sentinel, but the load's status stays connected

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

TEST_CASE("check_grid rejects a slack generator with a non-positive weight", "[check_grid]")
{
    LSGrid grid = make_valid_grid();
    LSGrid::StateRes st = grid.get_state();
    REQUIRE(gen_slack_weight(st).size() == 1);
    gen_slack_weight(st)[0] = 0.;  // gen 0 is the (only) slack

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

TEST_CASE("check_grid rejects an out-of-range remote-regulated bus id", "[check_grid]")
{
    LSGrid grid = make_valid_grid();
    LSGrid::StateRes st = grid.get_state();
    gen_regulated_bus(st)[0] = 9999;  // >> nb_bus

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::out_of_range);
}

TEST_CASE("check_grid accepts a valid pos_topo_vect permutation", "[check_grid]")
{
    LSGrid grid = make_valid_grid();
    // only the load and the generator carry a pos_topo_vect here: {0, 1} is a
    // valid permutation of [0, 2).
    IntVect load_pos(1); load_pos << 0;
    IntVect gen_pos(1);  gen_pos  << 1;
    grid.set_load_pos_topo_vect(load_pos);
    grid.set_gen_pos_topo_vect(gen_pos);
    CHECK_NOTHROW(grid.check_grid());
}

TEST_CASE("check_grid rejects a duplicated pos_topo_vect entry", "[check_grid]")
{
    LSGrid grid = make_valid_grid();
    IntVect load_pos(1); load_pos << 0;
    IntVect gen_pos(1);  gen_pos  << 0;  // collides with the load's slot
    grid.set_load_pos_topo_vect(load_pos);
    grid.set_gen_pos_topo_vect(gen_pos);
    CHECK_THROWS_AS(grid.check_grid(), std::runtime_error);
}

TEST_CASE("check_grid rejects an out-of-range pos_topo_vect entry", "[check_grid]")
{
    LSGrid grid = make_valid_grid();
    IntVect load_pos(1); load_pos << 0;
    IntVect gen_pos(1);  gen_pos  << 7;  // only 2 topo slots exist -> out of range
    grid.set_load_pos_topo_vect(load_pos);
    grid.set_gen_pos_topo_vect(gen_pos);
    CHECK_THROWS_AS(grid.check_grid(), std::out_of_range);
}

// --- solver policy parameters carried by the serialized state ----------------
// AlgoConfig is part of StateRes, so these values arrive from a pickle / binary
// file. refactor_every_n == 0 used to reach `iter % refactor_every_n` and kill
// the process with SIGFPE (integer division by zero); it must be rejected on the
// way in instead.

TEST_CASE("a state with refactor_every_n = 0 is rejected, not divided by", "[check_grid][algo_config]")
{
    LSGrid grid = make_valid_grid();
    LSGrid::StateRes st = grid.get_state();
    REQUIRE(ac_algo_int_params(st).size() >= 4);
    ac_algo_int_params(st)[3] = 0;  // refactor_every_n

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

TEST_CASE("a state with a negative refactor_every_n is rejected", "[check_grid][algo_config]")
{
    LSGrid grid = make_valid_grid();
    LSGrid::StateRes st = grid.get_state();
    ac_algo_int_params(st)[3] = -3;

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

TEST_CASE("a state with an out-of-range line-search parameter is rejected", "[check_grid][algo_config]")
{
    LSGrid grid = make_valid_grid();

    SECTION("ls_rho outside (0, 1)") {
        LSGrid::StateRes st = grid.get_state();
        REQUIRE(ac_algo_real_params(st).size() >= 6);
        ac_algo_real_params(st)[3] = 1.5;  // ls_rho
        LSGrid restored;
        CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
    }
    SECTION("non-finite max_dVa") {
        LSGrid::StateRes st = grid.get_state();
        ac_algo_real_params(st)[0] = std::numeric_limits<double>::quiet_NaN();
        LSGrid restored;
        CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
    }
}

TEST_CASE("a valid algo config still round-trips", "[check_grid][algo_config]")
{
    LSGrid grid = make_valid_grid();
    LSGrid::StateRes st = grid.get_state();
    LSGrid restored;
    CHECK_NOTHROW(restored.set_state(st));
    CHECK(restored.get_ac_algo_config().int_params[3] == grid.get_ac_algo_config().int_params[3]);
}

TEST_CASE("a solver returning a wrong-sized voltage vector is rejected, not indexed", "[check_grid][plugin]")
{
    ls2g::AlgorithmRegistry::instance().register_solver(
        "__check_grid_wrong_size__",
        [] { return std::unique_ptr<ls2g::BaseAlgo>(new WrongSizeAlgo()); });

    LSGrid grid = make_valid_grid();
    grid.change_algorithm("__check_grid_wrong_size__");
    // the name is not one of the built-ins, so this solver is an *external* one
    // (AlgorithmType::Custom) -- which is exactly the case the output check
    // guards; built-in solvers deliberately skip it.
    REQUIRE(grid.get_algo_type() == ls2g::AlgorithmType::Custom);
    CHECK_THROWS_AS(grid.ac_pf(flat_start(grid), 20, 1e-8), std::runtime_error);
}

TEST_CASE("a solver returning non-finite voltages is reported as non-converged", "[check_grid][plugin]")
{
    ls2g::AlgorithmRegistry::instance().register_solver(
        "__check_grid_nan__",
        [] { return std::unique_ptr<ls2g::BaseAlgo>(new NaNAlgo()); });

    LSGrid grid = make_valid_grid();
    grid.change_algorithm("__check_grid_nan__");
    REQUIRE(grid.get_algo_type() == ls2g::AlgorithmType::Custom);  // external solver
    CplxVect res;
    CHECK_NOTHROW(res = grid.ac_pf(flat_start(grid), 20, 1e-8));
    CHECK(res.size() == 0);                       // diverged => empty result
    CHECK_FALSE(grid.get_algo().converged());     // marked non-converged, no NaN leaked
}

// --- a grid using an external (plugin) solver survives serialization ---------
// The algorithm is stored by registry name, so a plugin solver round-trips (and
// a grid using one can be copied). Storing AlgorithmType instead collapsed every
// plugin onto AlgorithmType::Custom, which made both operations throw.

TEST_CASE("a grid using a plugin solver round-trips through get_state/set_state", "[check_grid][plugin]")
{
    ls2g::AlgorithmRegistry::instance().register_solver(
        "__roundtrip_plugin_algo__",
        [] { return std::unique_ptr<ls2g::BaseAlgo>(new PluginLikeAlgo()); });

    LSGrid grid = make_valid_grid();
    grid.change_algorithm("__roundtrip_plugin_algo__");
    REQUIRE(grid.get_algo_type() == ls2g::AlgorithmType::Custom);

    LSGrid::StateRes st = grid.get_state();
    LSGrid restored;
    REQUIRE_NOTHROW(restored.set_state(st));
    // the *same* plugin solver is selected again, not a silent fallback
    CHECK(restored.get_algo_type() == ls2g::AlgorithmType::Custom);
    CHECK(restored.get_algo().get_name() == "__roundtrip_plugin_algo__");
}

TEST_CASE("a grid using a plugin solver can be copied, and the copy uses that solver", "[check_grid][plugin]")
{
    ls2g::AlgorithmRegistry::instance().register_solver(
        "__copy_plugin_algo__",
        [] { return std::unique_ptr<ls2g::BaseAlgo>(new PluginLikeAlgo()); });

    LSGrid grid = make_valid_grid();
    grid.change_algorithm("__copy_plugin_algo__");

    LSGrid copied = grid.copy();

    // the copy is on the same solver...
    CHECK(copied.get_algo().get_name() == "__copy_plugin_algo__");
    CHECK(copied.get_algo_type() == ls2g::AlgorithmType::Custom);
    // ... and it is a real, working instance of it, not just the name carried
    // over: the factory has to have been called for this powerflow to run.
    // PluginLikeAlgo returns the input voltages, so convergence means it ran.
    const CplxVect V = copied.ac_pf(flat_start(copied), 20, 1e-8);
    REQUIRE(V.size() > 0);
    CHECK(copied.get_algo().converged());

    // the copy is independent: re-selecting on the copy leaves the original be
    copied.change_algorithm("NR_SparseLU");
    CHECK(copied.get_algo().get_name() == "NR_SparseLU");
    CHECK(grid.get_algo().get_name() == "__copy_plugin_algo__");
    CHECK(grid.get_algo_type() == ls2g::AlgorithmType::Custom);
}

TEST_CASE("a grid using a DC plugin solver can be copied, and the copy uses that solver", "[check_grid][plugin]")
{
    // DC plugins land in the other selector (LSGrid routes on BaseAlgo::IS_AC),
    // so copy() has to carry that one across by name too.
    ls2g::AlgorithmRegistry::instance().register_solver(
        "__copy_dc_plugin_algo__",
        [] { return std::unique_ptr<ls2g::BaseAlgo>(new DcPluginLikeAlgo()); });

    LSGrid grid = make_valid_grid();
    grid.change_algorithm("__copy_dc_plugin_algo__");
    REQUIRE(grid.get_dc_algo().get_name() == "__copy_dc_plugin_algo__");
    REQUIRE(grid.get_algo().get_name() == "NR_SparseLU");  // AC slot untouched

    LSGrid copied = grid.copy();
    CHECK(copied.get_dc_algo().get_name() == "__copy_dc_plugin_algo__");
    CHECK(copied.get_dc_algo_type() == ls2g::AlgorithmType::Custom);
    // again: a live instance, proven by running a DC powerflow with it
    const CplxVect V = copied.dc_pf(flat_start(copied), 1, 1e-8);
    REQUIRE(V.size() > 0);
    CHECK(copied.get_dc_algo().converged());
}

TEST_CASE("a grid using a DC plugin solver round-trips through get_state/set_state", "[check_grid][plugin]")
{
    ls2g::AlgorithmRegistry::instance().register_solver(
        "__roundtrip_dc_plugin_algo__",
        [] { return std::unique_ptr<ls2g::BaseAlgo>(new DcPluginLikeAlgo()); });

    LSGrid grid = make_valid_grid();
    grid.change_algorithm("__roundtrip_dc_plugin_algo__");

    LSGrid::StateRes st = grid.get_state();
    LSGrid restored;
    REQUIRE_NOTHROW(restored.set_state(st));
    CHECK(restored.get_dc_algo().get_name() == "__roundtrip_dc_plugin_algo__");
    CHECK(restored.get_dc_algo_type() == ls2g::AlgorithmType::Custom);
}

TEST_CASE("loading a grid whose solver is not registered names it and how to get it", "[check_grid][plugin]")
{
    ls2g::AlgorithmRegistry::instance().register_solver(
        "__vanishing_plugin_algo__",
        [] { return std::unique_ptr<ls2g::BaseAlgo>(new PluginLikeAlgo()); });

    LSGrid grid = make_valid_grid();
    grid.change_algorithm("__vanishing_plugin_algo__");
    LSGrid::StateRes st = grid.get_state();

    // simulate loading in a process where that plugin was never loaded
    std::get<LSGrid::AC_ALGO_NAME_ID>(st) = "__a_plugin_nobody_loaded__";

    LSGrid restored;
    REQUIRE_THROWS_AS(restored.set_state(st), std::runtime_error);
    try {
        LSGrid again;
        again.set_state(st);
    } catch (const std::runtime_error & e) {
        const std::string msg = e.what();
        CHECK(msg.find("__a_plugin_nobody_loaded__") != std::string::npos);   // names the solver
        CHECK(msg.find("load_algorithm_plugin") != std::string::npos);        // says what to do
    }
}

TEST_CASE("a grid whose solver is unavailable still loads without the algorithm", "[check_grid][plugin]")
{
    ls2g::AlgorithmRegistry::instance().register_solver(
        "__escape_hatch_plugin_algo__",
        [] { return std::unique_ptr<ls2g::BaseAlgo>(new PluginLikeAlgo()); });

    LSGrid grid = make_valid_grid();
    grid.change_algorithm("__escape_hatch_plugin_algo__");
    LSGrid::StateRes st = grid.get_state();
    std::get<LSGrid::AC_ALGO_NAME_ID>(st) = "__not_loaded_here__";

    // strict load refuses...
    LSGrid strict;
    CHECK_THROWS_AS(strict.set_state(st), std::runtime_error);

    // ... but the grid data itself loads fine when the solver is not restored,
    // keeping this object's default solver.
    LSGrid lenient;
    REQUIRE_NOTHROW(lenient.set_state(st, /*restore_algorithm=*/false));
    CHECK(lenient.get_algo().get_name() == "NR_SparseLU");   // the default
    CHECK(lenient.get_loads().nb() == grid.get_loads().nb());  // data is there
    CHECK_NOTHROW(lenient.check_grid());
    // and it can still run a powerflow with that default solver
    const CplxVect V = lenient.ac_pf(flat_start(lenient), 20, 1e-8);
    CHECK(V.size() > 0);
}

TEST_CASE("a corrupted solver name is escaped in the error message", "[check_grid][plugin]")
{
    LSGrid grid = make_valid_grid();
    LSGrid::StateRes st = grid.get_state();
    // raw bytes, as a corrupted file would hold: they must not reach what() as-is
    // (pybind11 would raise UnicodeDecodeError instead of our RuntimeError).
    std::get<LSGrid::AC_ALGO_NAME_ID>(st) = std::string("\xff\xfe bad", 7);

    LSGrid restored;
    try {
        restored.set_state(st);
        FAIL("expected set_state to throw");
    } catch (const std::runtime_error & e) {
        const std::string msg = e.what();
        CHECK(msg.find("\xff") == std::string::npos);   // escaped, not raw
        CHECK(msg.find("\\xff") != std::string::npos);
    }
}
