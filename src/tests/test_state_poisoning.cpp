// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Second round of "poisoned state" tests, on the element types the first one
// (test_check_grid.cpp) could not reach with its 3-bus / 2-line / 1-load /
// 1-gen fixture: transformers with a shift-dependent r/x table, SVCs, HVDC
// lines. The grid comes from case_exotic_elements.hpp (IEEE14 + SVC + storage +
// 3 HVDC lines + a phase-shifting transformer), so every container involved is
// non-empty.
//
// Same concern as test_check_grid.cpp: release wheels are built -O3 -DNDEBUG, so
// an inconsistency restored from a well-formed-but-poisoned pickle / binary file
// turns into an out-of-bounds read or write rather than an error. Every case
// below segfaulted before the corresponding fix; they also run under
// valgrind/ASan in the cpp_unit_tests workflow, which is what proves the
// rejection happens *before* anything is indexed.

#include <stdexcept>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "case_exotic_elements.hpp"

using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

LSGrid exotic() { return ls2g_test::make_exotic_elements_grid(); }

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
}

// --- accessors into the (deeply nested) StateRes tuple -----------------------
// TRAFO_ID -> TrafoContainer::StateRes = tuple<[0] TwoSidesContainer_rxh_A::StateRes,
//   ratio, is_tap_side1, shift, ignore_tap_side_for_shift, [5] shift_dependent_rx,
//   base_r, base_x, [8] rx_corr_alpha, [9] rx_corr_pct>
bool& trafo_shift_dependent_rx(LSGrid::StateRes & st)
{
    return std::get<5>(std::get<LSGrid::TRAFO_ID>(st));
}
std::vector<std::vector<real_type> >& trafo_rx_corr_alpha(LSGrid::StateRes & st)
{
    return std::get<8>(std::get<LSGrid::TRAFO_ID>(st));
}
std::vector<std::vector<real_type> >& trafo_rx_corr_pct(LSGrid::StateRes & st)
{
    return std::get<9>(std::get<LSGrid::TRAFO_ID>(st));
}
// SVC_ID -> SvcContainer::StateRes = tuple<OSC_PQ::StateRes, regulation_mode,
//   target_vm_pu, slope_pu, b_min, b_max, [6] regulated_bus_id>
std::vector<int>& svc_regulated_bus(LSGrid::StateRes & st)
{
    return std::get<6>(std::get<LSGrid::SVC_ID>(st));
}
// LINE_ID / TRAFO_ID -> [0] TwoSidesContainer_rxh_A::StateRes
//   -> [0] TwoSidesContainer::StateRes = tuple<ignore_status_global,
//          synch_status_both_side, names, [3] status_global, side_1, side_2>
template<std::size_t BRANCH_ID>
std::vector<bool>& branch_status_global(LSGrid::StateRes & st)
{
    return std::get<3>(std::get<0>(std::get<0>(std::get<BRANCH_ID>(st))));
}
// HVDC_ID -> HvdcLineContainer::StateRes = tuple<[0] TwoSidesContainer::StateRes, ...>
std::vector<bool>& hvdc_status_global(LSGrid::StateRes & st)
{
    return std::get<3>(std::get<0>(std::get<LSGrid::HVDC_ID>(st)));
}
// SHUNT_ID -> ShuntContainer::StateRes -> [0] OneSideContainer_PQ::StateRes
//   -> [0] OneSideContainer::StateRes = tuple<names, bus_id, status, has_subid,
//          subid, [5] has_pos_topo_vect, [6] pos_topo_vect>
bool& shunt_has_pos_topo_vect(LSGrid::StateRes & st)
{
    return std::get<5>(std::get<0>(std::get<0>(std::get<LSGrid::SHUNT_ID>(st))));
}
std::vector<int>& shunt_pos_topo_vect(LSGrid::StateRes & st)
{
    return std::get<6>(std::get<0>(std::get<0>(std::get<LSGrid::SHUNT_ID>(st))));
}

} // namespace

TEST_CASE("the exotic-elements grid state round-trips and stays valid", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    LSGrid restored;
    REQUIRE_NOTHROW(restored.set_state(st));
    CHECK_NOTHROW(restored.check_grid());
    const CplxVect V = restored.ac_pf(flat_start(restored), 20, 1e-7);
    CHECK(V.size() == 14);
}

// ---------------------------------------------------------------------------
// TrafoContainer: the (alpha -> r/x correction) tables
// ---------------------------------------------------------------------------
// `_update_model_coeffs()` runs at the very end of TrafoContainer::set_state --
// long before check_grid() -- and reads rx_corr_alpha_[el_id] / rx_corr_pct_[el_id]
// for el_id in [0, nb()) with an unchecked std::vector::operator[]. A state
// declaring fewer tables than transformers therefore builds a std::vector out of
// heap memory past the end of the array and dereferences its pointers.

TEST_CASE("set_state rejects shift-dependent r/x with no correction table at all", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    REQUIRE(trafo_shift_dependent_rx(st));            // the fixture uses this feature
    REQUIRE(trafo_rx_corr_alpha(st).size() == 3u);
    trafo_rx_corr_alpha(st).clear();
    trafo_rx_corr_pct(st).clear();

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

TEST_CASE("set_state rejects a correction table shorter than the number of trafos", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    REQUIRE(trafo_rx_corr_alpha(st).size() == 3u);
    trafo_rx_corr_alpha(st).pop_back();               // 2 tables for 3 transformers

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

TEST_CASE("set_state rejects alpha / correction tables of different lengths", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    // one alpha sample without its correction: _shift_rx_corr_pct() indexes `ys`
    // with positions derived from `xs.size()`
    trafo_rx_corr_alpha(st)[0] = std::vector<real_type>{-0.5, 0., 0.5};
    trafo_rx_corr_pct(st)[0] = std::vector<real_type>{1., 2.};

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

TEST_CASE("set_state accepts a grid with no correction table and the feature off", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    trafo_shift_dependent_rx(st) = false;
    trafo_rx_corr_alpha(st).clear();
    trafo_rx_corr_pct(st).clear();

    LSGrid restored;
    REQUIRE_NOTHROW(restored.set_state(st));
    CHECK_NOTHROW(restored.check_grid());
}

// ---------------------------------------------------------------------------
// SvcContainer: the regulated bus id
// ---------------------------------------------------------------------------
// `regulated_bus_id_` is a gridmodel bus id used directly as an index in
// SvcContainer::set_vm (`id_grid_to_solver[...]`, then a write into V) and in
// LSGrid::fill_voltage_control_solver_data. Nothing bounded it: SvcContainer::init
// only checks the vector's length, and set_state copies it verbatim. It is the
// exact counterpart of GeneratorContainer's field of the same name, which
// check_grid() has always validated.

TEST_CASE("check_grid rejects an out-of-range svc regulated bus id (set_state)", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    REQUIRE(svc_regulated_bus(st).size() == 1u);
    svc_regulated_bus(st)[0] = 10000000;

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::out_of_range);
}

TEST_CASE("check_grid rejects a negative (non-sentinel) svc regulated bus id", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    svc_regulated_bus(st)[0] = -5;

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::out_of_range);
}

TEST_CASE("check_grid rejects an out-of-range svc regulated bus set through init_svcs", "[check_grid][state_poisoning]")
{
    // same hole, reached without any serialization: a grid loader (init_from_pypowsybl
    // builds the SVCs) passing a bus id its source file got wrong.
    LSGrid grid = exotic();
    const std::vector<int> mode{1};                    // VOLTAGE
    RealVect target_vm_pu(1), q_setpoint(1), slope(1), b_min(1), b_max(1);
    target_vm_pu << 1.0;
    q_setpoint << 0.;
    slope << 0.;
    b_min << -0.03;
    b_max << 0.03;
    Eigen::VectorXi regulated_bus(1), svc_bus(1);
    regulated_bus << 10000000;                         // out of range
    svc_bus << 5;
    grid.init_svcs(mode, target_vm_pu, q_setpoint, slope, b_min, b_max, regulated_bus, svc_bus);

    CHECK_THROWS_AS(grid.check_grid(), std::out_of_range);
}

TEST_CASE("check_grid accepts the -1 sentinel as an svc regulated bus", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    svc_regulated_bus(st)[0] = -1;                     // regulates no bus

    LSGrid restored;
    REQUIRE_NOTHROW(restored.set_state(st));
    CHECK_NOTHROW(restored.check_grid());
}

// ---------------------------------------------------------------------------
// TwoSidesContainer: status_global_
// ---------------------------------------------------------------------------
// `nb()` is side_1_.nb(), NOT status_global_.size(), and set_tsc_state tied the
// two together nowhere -- while status_global_[el_id] (unchecked operator[], with
// el_id bounded by nb()) is read all over the class and the batch solvers.

TEST_CASE("set_state rejects an empty status_global on the powerlines", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    REQUIRE(!branch_status_global<LSGrid::LINE_ID>(st).empty());
    branch_status_global<LSGrid::LINE_ID>(st).clear();

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

TEST_CASE("set_state rejects a short status_global on the powerlines", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    branch_status_global<LSGrid::LINE_ID>(st).pop_back();

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

TEST_CASE("set_state rejects a wrong-sized status_global on the trafos", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    branch_status_global<LSGrid::TRAFO_ID>(st).clear();

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

TEST_CASE("set_state rejects a wrong-sized status_global on the hvdc lines", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    REQUIRE(!hvdc_status_global(st).empty());
    hvdc_status_global(st).push_back(true);            // one too many

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}

// ---------------------------------------------------------------------------
// pos_topo_vect on an element that is not in the grid2op topology vector
// ---------------------------------------------------------------------------
// check_grid() proves the collected positions are a permutation of [0, K) where
// K is how many it collected. update_topo() sizes its caller arrays as
// `nb loads + nb gens + nb storages + 2 * nb lines + 2 * nb trafos` and indexes
// them with those very positions. Letting a shunt / sgen / svc / hvdc position
// into the same pot makes K bigger than that, so a *validated* position could
// still be past the end of the arrays update_topo() reads.

TEST_CASE("check_grid rejects a pos_topo_vect on a shunt", "[check_grid][state_poisoning]")
{
    LSGrid grid = exotic();
    LSGrid::StateRes st = grid.get_state();
    const std::size_t nb_shunt = std::get<1>(std::get<0>(std::get<0>(std::get<LSGrid::SHUNT_ID>(st)))).size();
    REQUIRE(nb_shunt > 0u);
    shunt_has_pos_topo_vect(st) = true;
    std::vector<int> pos(nb_shunt);
    for(std::size_t i = 0; i < nb_shunt; ++i) pos[i] = static_cast<int>(i);
    shunt_pos_topo_vect(st) = pos;

    LSGrid restored;
    CHECK_THROWS_AS(restored.set_state(st), std::runtime_error);
}
