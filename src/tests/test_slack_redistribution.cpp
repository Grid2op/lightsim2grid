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
