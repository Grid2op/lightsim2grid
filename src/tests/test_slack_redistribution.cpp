// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// slack_redistribution::distribute is OpenLoadFlow's GenerationActivePowerDistributionStep:
// share a mismatch by weight, clamp each unit to [min_p, max_p], a clamped unit leaves the
// pool and what it could not take is shared again. Pure function: tested on its own.

#include <cmath>
#include <limits>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "element_container/SlackRedistribution.hpp"
#include "LSGrid.hpp"

using ls2g::real_type;
using ls2g::slack_redistribution::Participant;
using ls2g::slack_redistribution::Report;
using ls2g::slack_redistribution::UnitKind;
using ls2g::slack_redistribution::distribute;
using ls2g::slack_redistribution::default_eps_mw;

namespace {

const real_type NaN = std::numeric_limits<real_type>::quiet_NaN();

Participant unit(int id, real_type inj, real_type w, real_type min_p, real_type max_p){
    Participant p;
    p.kind = UnitKind::GENERATOR;
    p.el_id = id;
    p.bus = id;
    p.injection_mw = inj;
    p.weight = w;
    p.min_p_mw = min_p;
    p.max_p_mw = max_p;
    return p;
}

real_type total(const std::vector<real_type> & v){
    real_type s = 0.;
    for(real_type x : v) s += x;
    return s;
}

}  // namespace

TEST_CASE("distribute: no bound, shares by weight in one round", "[slack_redistribution]"){
    std::vector<Participant> units = {unit(0, 10., 1., NaN, NaN), unit(1, 20., 3., NaN, NaN)};
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    const Report rep = distribute(units, 40., default_eps_mw, new_inj, sat);
    REQUIRE(rep.nb_rounds == 1);
    REQUIRE(rep.nb_saturated == 0);
    REQUIRE_FALSE(rep.all_saturated);
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(20., 1e-12));
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(50., 1e-12));
    CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(0., 1e-12));
}

TEST_CASE("distribute: a saturated unit leaves the pool, the rest is shared again", "[slack_redistribution]"){
    // equal weights, 40 MW to place: 10 each, but unit 1 can only take 3 and unit 2 only 5
    std::vector<Participant> units = {unit(0, 0., 1., NaN, NaN), unit(1, 0., 1., NaN, 3.),
                                      unit(2, 0., 1., NaN, 5.), unit(3, 0., 1., NaN, NaN)};
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    const Report rep = distribute(units, 40., default_eps_mw, new_inj, sat);
    REQUIRE(rep.nb_saturated == 2);
    REQUIRE_FALSE(rep.all_saturated);
    CHECK(sat[1] == 1);
    CHECK(sat[2] == 1);
    CHECK(sat[0] == 0);
    CHECK(sat[3] == 0);
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(3., 1e-12));
    CHECK_THAT(new_inj[2], Catch::Matchers::WithinAbs(5., 1e-12));
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(16., 1e-12));
    CHECK_THAT(new_inj[3], Catch::Matchers::WithinAbs(16., 1e-12));
    CHECK_THAT(total(new_inj), Catch::Matchers::WithinAbs(40., 1e-12));
    CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(0., 1e-12));
}

TEST_CASE("distribute: a negative mismatch clamps at min_p", "[slack_redistribution]"){
    std::vector<Participant> units = {unit(0, 10., 1., 8., NaN), unit(1, 10., 1., NaN, NaN)};
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    const Report rep = distribute(units, -10., default_eps_mw, new_inj, sat);
    REQUIRE(rep.nb_saturated == 1);
    CHECK(sat[0] == 1);
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(8., 1e-12));
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(2., 1e-12));
}

TEST_CASE("distribute: a unit already beyond its bound is never moved back", "[slack_redistribution]"){
    std::vector<Participant> units = {unit(0, 12., 1., NaN, 10.), unit(1, 0., 1., NaN, NaN)};
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    distribute(units, 4., default_eps_mw, new_inj, sat);
    CHECK(sat[0] == 1);
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(12., 1e-12));  // not clamped DOWN to 10
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(4., 1e-12));
}

TEST_CASE("distribute: every unit saturated keeps them all in the slack", "[slack_redistribution]"){
    std::vector<Participant> units = {unit(0, 0., 1., NaN, 1.), unit(1, 0., 1., NaN, 2.)};
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    const Report rep = distribute(units, 10., default_eps_mw, new_inj, sat);
    REQUIRE(rep.all_saturated);
    CHECK(rep.nb_saturated == 2);
    CHECK(sat[0] == 0);  // cleared: nobody leaves the slack
    CHECK(sat[1] == 0);
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(1., 1e-12));
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(2., 1e-12));
    CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(7., 1e-12));
}

TEST_CASE("distribute: a unit injecting never crosses 0 MW on the way down", "[slack_redistribution]"){
    // unit 1 has min_p < 0 < max_p: its share (-10) would take it to -8, OLF stops it at 0
    // and the rest goes to unit 0
    std::vector<Participant> units = {unit(0, 30., 1., NaN, NaN), unit(1, 2., 1., -50., 50.)};
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    const Report rep = distribute(units, -20., default_eps_mw, new_inj, sat);
    REQUIRE(rep.nb_rounds == 2);
    REQUIRE(rep.nb_saturated == 1);
    REQUIRE_FALSE(rep.all_saturated);
    CHECK(sat[1] == 1);
    CHECK(sat[0] == 0);
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(0., 1e-12));
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(12., 1e-12));
    CHECK_THAT(total(new_inj), Catch::Matchers::WithinAbs(12., 1e-12));
    CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(0., 1e-12));
}

TEST_CASE("distribute: a unit drawing power never crosses 0 MW on the way up", "[slack_redistribution]"){
    // a charging storage unit (injection -5, range [-50, 50]): a positive mismatch stops it at 0
    std::vector<Participant> units = {unit(0, 30., 1., NaN, NaN), unit(1, -5., 1., -50., 50.)};
    units[1].kind = UnitKind::STORAGE;
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    const Report rep = distribute(units, 20., default_eps_mw, new_inj, sat);
    REQUIRE(rep.nb_rounds == 2);
    REQUIRE(rep.nb_saturated == 1);
    CHECK(sat[1] == 1);
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(0., 1e-12));
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(45., 1e-12));
    CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(0., 1e-12));
}

TEST_CASE("distribute: a unit moves freely on its own side of 0 MW", "[slack_redistribution]"){
    // drawing 5 MW, pushed down: only its own min_p bounds it
    std::vector<Participant> units = {unit(0, -5., 1., -50., 50.)};
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    Report rep = distribute(units, -20., default_eps_mw, new_inj, sat);
    REQUIRE(rep.nb_saturated == 0);
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(-25., 1e-12));
    // ... and its min_p still holds
    rep = distribute(units, -60., default_eps_mw, new_inj, sat);
    REQUIRE(rep.all_saturated);
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(-50., 1e-12));
    CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(-15., 1e-12));
}

TEST_CASE("distribute: a unit at exactly 0 MW only moves up", "[slack_redistribution]"){
    std::vector<Participant> units = {unit(0, 30., 1., NaN, NaN), unit(1, 0., 1., -50., 50.)};
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    Report rep = distribute(units, -10., default_eps_mw, new_inj, sat);
    CHECK(sat[1] == 1);
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(0., 1e-12));
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(20., 1e-12));
    rep = distribute(units, 10., default_eps_mw, new_inj, sat);
    REQUIRE(rep.nb_saturated == 0);
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(5., 1e-12));
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(35., 1e-12));
}

TEST_CASE("distribute: nothing to share is a no-op", "[slack_redistribution]"){
    std::vector<Participant> units = {unit(0, 5., 1., NaN, NaN)};
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    const Report rep = distribute(units, 1e-9, default_eps_mw, new_inj, sat);
    CHECK(rep.nb_rounds == 0);
    CHECK(rep.nb_participants == 1);
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(5., 1e-12));
    std::vector<Participant> nobody;
    const Report rep2 = distribute(nobody, 100., default_eps_mw, new_inj, sat);
    CHECK(rep2.nb_rounds == 0);
    CHECK(new_inj.empty());
}

namespace {

// buses 0-1-2-3 in a row plus a leaf bus 4 off bus 1 (line 3); 60 MW of load on bus 3 and
// 20 MW on bus 4. Gen 0 (bus 0) is the slack; gen 1 (bus 2) sits at its max_p, out of the
// slack, flagged (or not) "can participate in the slack" with the same weight as gen 0.
ls2g::LSGrid make_capped_grid(bool flagged)
{
    ls2g::LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(5u, 1, ls2g::RealVect::Constant(5, 138.), 0, 0);
    Eigen::VectorXi fr(4), to(4);
    fr << 0, 1, 2, 1;
    to << 1, 2, 3, 4;
    grid.init_powerlines(ls2g::RealVect::Constant(4, 0.01), ls2g::RealVect::Constant(4, 0.1),
                         ls2g::CplxVect::Zero(4), fr, to);
    ls2g::RealVect load_p(2), load_q(2);
    load_p << 60., 20.;
    load_q << 10., 5.;
    Eigen::VectorXi load_bus(2);
    load_bus << 3, 4;
    grid.init_loads(load_p, load_q, load_bus);
    ls2g::RealVect gen_p(2), gen_v(2), gen_q(2), min_q(2), max_q(2), min_p(2), max_p(2);
    gen_p << 30., 40.;
    gen_v << 1.02, 1.02;
    gen_q << 0., 0.;
    min_q << -1e3, -1e3;
    max_q << 1e3, 1e3;
    min_p << 0., 0.;
    max_p << 500., 40.;
    Eigen::VectorXi gen_bus(2);
    gen_bus << 0, 2;
    grid.init_generators_full(gen_p, gen_v, gen_q, std::vector<bool>{true, true}, min_q, max_q, gen_bus);
    grid.set_gen_p_limits(min_p, max_p);
    grid.add_gen_slackbus(0, 0.5);
    if(flagged){
        ls2g::RealVect w(2);
        w << 0., 0.5;
        grid.set_gen_can_participate_slack(std::vector<bool>{false, true}, w);
    }
    grid.tell_solver_need_reset();
    return grid;
}

}  // namespace

TEST_CASE("a unit flagged can_participate_slack takes a pre-pass share away from its limit", "[slack_redistribution]"){
    SECTION("flagged: it takes half of what the island takes out, and stays out of the slack"){
        ls2g::LSGrid grid = make_capped_grid(true);
        grid.deactivate_powerline(3);
        const Report rep = grid.consider_only_main_component(true);
        CHECK(rep.nb_participants == 2);
        CHECK_THAT(rep.mismatch_mw, Catch::Matchers::WithinAbs(-20., 1e-9));
        CHECK_THAT(grid.get_gen_target_p()(0), Catch::Matchers::WithinAbs(20., 1e-9));
        CHECK_THAT(grid.get_gen_target_p()(1), Catch::Matchers::WithinAbs(30., 1e-9));
        CHECK(grid.get_generators().is_slack(0));
        CHECK_FALSE(grid.get_generators().is_slack(1));
    }
    SECTION("not flagged: the slack unit takes it all"){
        ls2g::LSGrid grid = make_capped_grid(false);
        grid.deactivate_powerline(3);
        const Report rep = grid.consider_only_main_component(true);
        CHECK(rep.nb_participants == 1);
        CHECK_THAT(grid.get_gen_target_p()(0), Catch::Matchers::WithinAbs(10., 1e-9));
        CHECK_THAT(grid.get_gen_target_p()(1), Catch::Matchers::WithinAbs(40., 1e-9));
    }
    SECTION("the flag is kept by a copy"){
        ls2g::LSGrid grid = make_capped_grid(true);
        const ls2g::LSGrid other = grid.copy();
        CHECK(other.get_generators().get_can_participate_slack(1));
        CHECK_THAT(other.get_generators().get_can_participate_slack_weight(1), Catch::Matchers::WithinAbs(0.5, 1e-15));
    }
}
