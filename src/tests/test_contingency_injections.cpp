// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Contingencies that trip shunts, loads and generators, not only branches:
// ContingencyAnalysis' add_contingency / add_contingency_shunt / _load / _gen, and
// ScenarioSweep's set_contingency_shunts / set_contingency_loads. The oracle is a one-off
// powerflow on a copy of the grid with the elements really deactivated.

#include <complex>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "powerflow_algorithm/outer_loop/ShuntVoltageControlLoop.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::ContingencyAnalysis;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::ScenarioSweep;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using BoolMat = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const real_type TOL = 1e-10;

// buses 0-1-2-3 in a line, a second line 1-2; generators: 0 the slack at bus 0, 1 a PV
// one at bus 2 (the only one regulating it), 2 a PQ injection at bus 3; loads at 2 and
// 3; shunts at 1 and 3. `distributed`: generators 0 and 1 share the slack.
LSGrid make_grid(bool distributed = false)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(4, 1, RealVect::Constant(4, 138.), 0, 0);

    const int n_line = 4;
    RealVect branch_r = RealVect::Constant(n_line, 0.01);
    RealVect branch_x = RealVect::Constant(n_line, 0.1);
    Eigen::VectorXi from_id(n_line), to_id(n_line);
    from_id << 0, 1, 2, 1;
    to_id << 1, 2, 3, 2;
    branch_x(3) = 0.15;
    grid.init_powerlines(branch_r, branch_x, CplxVect::Constant(n_line, {0., 0.02}), from_id, to_id);

    RealVect load_p(2), load_q(2);
    load_p << 30., 50.;
    load_q << 8., 12.;
    Eigen::VectorXi load_bus(2);
    load_bus << 2, 3;
    grid.init_loads(load_p, load_q, load_bus);

    RealVect gen_p(3), gen_v(3), gen_q(3), gen_min_q(3), gen_max_q(3);
    Eigen::VectorXi gen_bus(3);
    gen_p << 0., 40., 15.;
    gen_v << 1.03, 1.02, 1.0;
    gen_q << 0., 0., 4.;
    gen_min_q << -1000., -1000., -1000.;
    gen_max_q << 1000., 1000., 1000.;
    gen_bus << 0, 2, 3;
    const std::vector<bool> regulates = {true, true, false};
    grid.init_generators_full(gen_p, gen_v, gen_q, regulates, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    if(distributed) grid.add_gen_slackbus(1, 1.);

    RealVect shunt_p(2), shunt_q(2);
    shunt_p << 0.5, 1.;
    shunt_q << 10., -25.;   // an absorbing one, a producing one
    Eigen::VectorXi shunt_bus(2);
    shunt_bus << 1, 3;
    grid.init_shunt(shunt_p, shunt_q, shunt_bus);

    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    grid.tell_solver_need_reset();
    return grid;
}

CplxVect flat(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

struct Trip {
    std::vector<int> branches, shunts, loads, gens;
};

// the one-off oracle: a copy of the grid with the elements really out
CplxVect reference(const LSGrid & base, const Trip & trip, bool ac = true)
{
    LSGrid grid(base);
    for(int b : trip.branches) grid.deactivate_powerline(b);
    for(int s : trip.shunts) grid.deactivate_shunt(s);
    for(int l : trip.loads) grid.deactivate_load(l);
    for(int g : trip.gens) grid.deactivate_gen(g);
    grid.tell_solver_need_reset();
    return ac ? grid.ac_pf(flat(grid), 30, TOL) : grid.dc_pf(flat(grid), 30, TOL);
}

void require_close(const CplxVect & ref, const Eigen::Ref<const CplxVect> & got, real_type tol = 1e-8)
{
    REQUIRE(ref.size() == got.size());
    for(Eigen::Index b = 0; b < ref.size(); ++b){
        CHECK(std::abs(ref(b) - got(b)) <= tol);
    }
}

// the row a ContingencyAnalysis gave a contingency (rows follow my_contingencies())
Eigen::Index row_of(const ContingencyAnalysis & ca, const Trip & trip)
{
    const auto conts = ca.my_contingencies();
    for(size_t row = 0; row < conts.size(); ++row){
        if(std::get<0>(conts[row]) == trip.branches && std::get<1>(conts[row]) == trip.shunts &&
           std::get<2>(conts[row]) == trip.loads && std::get<3>(conts[row]) == trip.gens) return static_cast<Eigen::Index>(row);
    }
    FAIL("contingency not registered");
    return -1;
}

}  // namespace

TEST_CASE("ContingencyAnalysis trips shunts, loads and generators like a one-off powerflow", "[batch][contingency_injections]")
{
    for(bool distributed : {false, true}){
        DYNAMIC_SECTION("distributed slack " << distributed){
            const LSGrid grid = make_grid(distributed);
            ContingencyAnalysis ca(grid);
            ca.change_algorithm(AlgorithmType::NR_SparseLU);
            const std::vector<Trip> trips = {
                {{}, {0}, {}, {}},          // an absorbing shunt
                {{}, {1}, {}, {}},          // a producing one
                {{}, {}, {0}, {}},          // a load
                {{}, {}, {}, {2}},          // a PQ generator: its P and Q leave
                {{}, {}, {}, {1}},          // the last generator regulating bus 2: PV -> PQ
                {{1}, {0}, {1}, {}},        // a branch, a shunt and a load at once
                {{3}, {1}, {0, 1}, {2}},    // everything together
            };
            ca.add_contingency_shunt({0});
            ca.add_contingency_shunt({1});
            ca.add_contingency_load({0});
            ca.add_contingency_gen({2});
            ca.add_contingency_gen({1});
            ca.add_contingency({1}, {0}, {1}, {});
            ca.add_contingency({3}, {1}, {0, 1}, {2});
            REQUIRE(ca.my_defaults().size() == trips.size());
            ca.compute(flat(grid), 30, TOL);
            const auto & V = ca.get_voltages();
            const CplxVect V_n = reference(grid, Trip());
            for(const Trip & trip : trips){
                const Eigen::Index row = row_of(ca, trip);
                REQUIRE(ca.converged_mask()[static_cast<size_t>(row)]);
                const CplxVect ref = reference(grid, trip);
                require_close(ref, V.row(row).transpose());
                // the trip does change the state: the match above is not the base case's
                CHECK((ref - V_n).cwiseAbs().maxCoeff() > 1e-4);
            }
        }
    }
}

TEST_CASE("a contingency is its whole set of elements", "[batch][contingency_injections]")
{
    const LSGrid grid = make_grid();
    ContingencyAnalysis ca(grid);
    ca.add_nk({1});
    ca.add_contingency({1}, {0}, {}, {});
    ca.add_contingency({1}, {0}, {}, {});   // registered once
    ca.add_contingency({}, {}, {1}, {});
    ca.add_contingency_load({1});           // the same one
    REQUIRE(ca.my_defaults().size() == 3);
    // branch-only first (the set orders field by field, branches first)
    const auto conts = ca.my_contingencies();
    CHECK(std::get<0>(conts[0]).empty());
    CHECK(std::get<2>(conts[0]) == std::vector<int>({1}));
    CHECK(std::get<0>(conts[1]) == std::vector<int>({1}));
    CHECK(std::get<1>(conts[1]).empty());
    CHECK(std::get<1>(conts[2]) == std::vector<int>({0}));
    // my_defaults_vect lists the branches only
    const auto branches = ca.my_defaults_vect();
    CHECK(branches[0].empty());
    CHECK(branches[1] == std::vector<int>({1}));
    CHECK(branches[2] == std::vector<int>({1}));

    CHECK(ca.remove_contingency({1}, {0}, {}, {}));
    CHECK_FALSE(ca.remove_contingency({1}, {0}, {}, {}));
    CHECK(ca.remove_nk({1}));
    REQUIRE(ca.my_defaults().size() == 1);

    CHECK_THROWS_AS(ca.add_contingency({}, {}, {}, {}), std::runtime_error);
    CHECK_THROWS_AS(ca.add_contingency_shunt({2}), std::runtime_error);
    CHECK_THROWS_AS(ca.add_contingency_load({-1}), std::runtime_error);
    CHECK_THROWS_AS(ca.add_contingency_gen({3}), std::runtime_error);
    CHECK_THROWS_AS(ca.add_contingency({4}, {}, {}, {}), std::runtime_error);
}

TEST_CASE("ContingencyAnalysis trips shunts, loads and generators in DC", "[batch][contingency_injections][dc]")
{
    const LSGrid grid = make_grid(true);
    ContingencyAnalysis ca(grid);
    ca.change_algorithm(AlgorithmType::DC_SparseLU);
    const std::vector<Trip> trips = {
        {{}, {1}, {}, {}},
        {{}, {}, {0}, {}},
        {{}, {}, {}, {1}},
        {{3}, {0}, {1}, {2}},
    };
    for(const Trip & trip : trips) ca.add_contingency(trip.branches, trip.shunts, trip.loads, trip.gens);
    ca.compute(flat(grid), 30, TOL);
    const auto & V = ca.get_voltages();
    for(const Trip & trip : trips){
        const Eigen::Index row = row_of(ca, trip);
        REQUIRE(ca.converged_mask()[static_cast<size_t>(row)]);
        const CplxVect ref = reference(grid, trip, false);
        REQUIRE(ref.size() == V.cols());
        for(Eigen::Index b = 0; b < ref.size(); ++b){
            CHECK(std::arg(V(row, b)) == Approx(std::arg(ref(b))).margin(1e-8));
        }
    }
}

TEST_CASE("ScenarioSweep trips shunts and loads per row like a one-off powerflow", "[batch][scenario_sweep][contingency_injections]")
{
    const LSGrid grid = make_grid(true);
    ScenarioSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    BoolMat shunt_off = BoolMat::Constant(4, 2, false);
    BoolMat load_off = BoolMat::Constant(4, 2, false);
    BoolMat line_off = BoolMat::Constant(4, 4, false);
    shunt_off(0, 1) = true;
    load_off(1, 0) = true;
    shunt_off(2, 0) = true;
    load_off(2, 1) = true;
    line_off(2, 3) = true;
    // row 3 trips nothing
    sweep.set_contingency_shunts(shunt_off);
    sweep.set_contingency_loads(load_off);
    sweep.set_contingency_lines(line_off);
    sweep.compute(flat(grid), 30, TOL);
    const std::vector<Trip> trips = {
        {{}, {1}, {}, {}},
        {{}, {}, {0}, {}},
        {{3}, {0}, {1}, {}},
        {{}, {}, {}, {}},
    };
    const auto & V = sweep.get_voltages();
    for(size_t row = 0; row < trips.size(); ++row){
        REQUIRE(sweep.converged_mask()[row]);
        require_close(reference(grid, trips[row]), V.row(static_cast<Eigen::Index>(row)).transpose());
    }
}

TEST_CASE("a ScenarioSweep load trip uses the row's own load", "[batch][scenario_sweep][contingency_injections]")
{
    // the row's own load_p / load_q is what leaves, not the grid's: a trip of a load the
    // row also gives a value to is the same as the trip on the grid's value
    const LSGrid grid = make_grid();
    ScenarioSweep sweep(grid);
    RealMat load_p(1, 2), load_q(1, 2);
    load_p << 300., 50.;    // load 0 changed by the row
    load_q << 80., 12.;
    BoolMat load_off = BoolMat::Constant(1, 2, false);
    load_off(0, 0) = true;
    sweep.modify_load_p(load_p);
    sweep.modify_load_q(load_q);
    sweep.set_contingency_loads(load_off);
    sweep.compute(flat(grid), 30, TOL);
    REQUIRE(sweep.converged_mask()[0]);
    require_close(reference(grid, {{}, {}, {0}, {}}), sweep.get_voltages().row(0).transpose());
}

TEST_CASE("ContingencyAnalysis and ScenarioSweep agree on the same trips, slack redistributed", "[batch][contingency_injections][redistribute]")
{
    // the slack pre-pass counts what a tripped load / shunt / generator takes out of the
    // balance (term (e) / (a)): both classes, on the same trip, solve the same system
    const LSGrid grid = make_grid(true);
    ContingencyAnalysis ca(grid);
    ca.change_algorithm(AlgorithmType::NR_SparseLU);
    ca.set_redistribute_slack(true);
    ca.add_contingency({}, {1}, {1}, {});
    ca.add_contingency({}, {}, {}, {2});
    ca.compute(flat(grid), 30, TOL);

    ScenarioSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    sweep.set_redistribute_slack(true);
    BoolMat shunt_off = BoolMat::Constant(2, 2, false);
    BoolMat load_off = BoolMat::Constant(2, 2, false);
    BoolMat gen_off = BoolMat::Constant(2, 3, false);
    shunt_off(0, 1) = true;
    load_off(0, 1) = true;
    gen_off(1, 2) = true;
    sweep.set_contingency_shunts(shunt_off);
    sweep.set_contingency_loads(load_off);
    sweep.set_contingency_gens(gen_off);
    sweep.compute(flat(grid), 30, TOL);

    const Eigen::Index ca_row0 = row_of(ca, {{}, {1}, {1}, {}});
    const Eigen::Index ca_row1 = row_of(ca, {{}, {}, {}, {2}});
    require_close(ca.get_voltages().row(ca_row0).transpose(), sweep.get_voltages().row(0).transpose(), 1e-9);
    require_close(ca.get_voltages().row(ca_row1).transpose(), sweep.get_voltages().row(1).transpose(), 1e-9);
}

TEST_CASE("tripping a regulating shunt is refused with outer loops", "[batch][contingency_injections][outer]")
{
    LSGrid grid = make_grid();
    std::vector<real_type> p(4, 0.), q;
    for(int k = 1; k <= 4; ++k) q.push_back(-5. * k);
    grid.set_shunt_sections(1, 2, p, q);
    grid.set_shunt_section_regulation(1, true, 1.0, 0., 3);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ls2g::ShuntVoltageControlLoop>());

    ContingencyAnalysis ca(grid);
    ca.add_contingency_shunt({1});
    CHECK_THROWS_AS(ca.compute(flat(grid), 30, TOL), std::runtime_error);

    // the non-regulating one is fine
    ContingencyAnalysis ca_ok(grid);
    ca_ok.add_contingency_shunt({0});
    ca_ok.compute(flat(grid), 30, TOL);
    CHECK(ca_ok.converged_mask()[0]);

    ScenarioSweep sweep(grid);
    BoolMat shunt_off = BoolMat::Constant(1, 2, false);
    shunt_off(0, 1) = true;
    sweep.set_contingency_shunts(shunt_off);
    CHECK_THROWS_AS(sweep.compute(flat(grid), 30, TOL), std::runtime_error);
}

TEST_CASE("ContingencyAnalysis trips loads and generators with outer loops", "[batch][contingency_injections][outer]")
{
    // the rows run the grid's outer loops on the edited system (the tripped load leaves
    // the row's injection, the tripped generator its injection and its bus' pinning)
    LSGrid grid = make_grid(true);
    grid.change_algorithm("NROuter_SparseLU");
    ContingencyAnalysis ca(grid);
    const std::vector<Trip> trips = {
        {{}, {0}, {}, {}},
        {{}, {}, {1}, {}},
        {{}, {}, {}, {1}},
        {{1}, {1}, {0}, {2}},
    };
    for(const Trip & trip : trips) ca.add_contingency(trip.branches, trip.shunts, trip.loads, trip.gens);
    ca.compute(flat(grid), 30, TOL);
    const auto & V = ca.get_voltages();
    for(const Trip & trip : trips){
        const Eigen::Index row = row_of(ca, trip);
        REQUIRE(ca.converged_mask()[static_cast<size_t>(row)]);
        require_close(reference(grid, trip), V.row(row).transpose());
    }
}
