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

TEST_CASE("distribute: every slack unit saturated keeps them in the slack, a flagged unit still moving", "[slack_redistribution]"){
    // unit 0 is the only slack unit and can take 2 MW; unit 1 is only flagged "can
    // participate": it takes the rest, but the solve cannot distribute on it
    std::vector<Participant> units = {unit(0, 30., 1., 0., 32.), unit(1, 40., 1., 0., 100.)};
    units[1].in_slack = false;
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    const Report rep = distribute(units, 10., default_eps_mw, new_inj, sat);
    REQUIRE(rep.all_saturated);
    CHECK(sat[0] == 0);  // cleared: the slack is not left empty
    CHECK(sat[1] == 0);
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(32., 1e-12));
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(48., 1e-12));
    CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(0., 1e-12));
    SECTION("a second slack unit with room: the saturated one leaves the slack as usual"){
        units.push_back(unit(2, 20., 1., 0., 100.));
        const Report rep2 = distribute(units, 10., default_eps_mw, new_inj, sat);
        REQUIRE_FALSE(rep2.all_saturated);
        CHECK(sat[0] == 1);
        CHECK(sat[2] == 0);
    }
}

TEST_CASE("distribute with an overshoot: every slack unit saturated keeps them in the slack", "[slack_redistribution]"){
    // the same as above, but the flagged unit carries an overshoot (it sat at its min_p,
    // 5 MW below it): the bisection path has to keep the slack from being emptied too
    std::vector<Participant> units = {unit(0, 30., 1., 0., 32.), unit(1, 40., 1., 40., 100.)};
    units[1].in_slack = false;
    units[1].overshoot_mw = -5.;
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    const Report rep = distribute(units, 10., default_eps_mw, new_inj, sat);
    REQUIRE(rep.all_saturated);
    CHECK(sat[0] == 0);  // cleared: the slack is not left empty
    CHECK(sat[1] == 0);
    CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(32., 1e-6));
    CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(48., 1e-6));
    CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(0., 1e-6));
    SECTION("a second slack unit with room: the saturated one leaves the slack as usual"){
        units.push_back(unit(2, 20., 1., 0., 100.));
        const Report rep2 = distribute(units, 10., default_eps_mw, new_inj, sat);
        REQUIRE_FALSE(rep2.all_saturated);
        CHECK(sat[0] == 1);
        CHECK(sat[2] == 0);
    }
}

TEST_CASE("distribute: no weight to share on, the whole mismatch is not distributed", "[slack_redistribution]"){
    // both paths agree: nothing moved, so nothing was distributed
    std::vector<Participant> units = {unit(0, 30., 0., 0., 100.), unit(1, 40., 0., 0., 40.)};
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    SECTION("round-based path"){
        const Report rep = distribute(units, -10., default_eps_mw, new_inj, sat);
        CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(30., 1e-12));
        CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(40., 1e-12));
        CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(-10., 1e-12));
    }
    SECTION("overshoot path"){
        units[1].in_slack = false;
        units[1].overshoot_mw = 5.;
        const Report rep = distribute(units, -10., default_eps_mw, new_inj, sat);
        CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(30., 1e-12));
        CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(40., 1e-12));
        CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(-10., 1e-12));
    }
}

TEST_CASE("distribute: the overshoot left is what the shift did not use up", "[slack_redistribution]"){
    // unit 1 sat 25 MW beyond its max_p in the reference solve: carried into the next
    // distribution, what is left of it makes sharing d1 then d2 end where sharing d1 + d2 at
    // once does. A distribution never makes an overshoot of its own.
    std::vector<real_type> new_inj, left;
    std::vector<char> sat;
    std::vector<Participant> units = {unit(0, 30., 1., 0., 500.), unit(1, 40., 1., 0., 40.)};
    units[1].in_slack = false;
    units[1].overshoot_mw = 25.;
    SECTION("used up in part, then in full"){
        // -20: shift -40, unit 1 stays at 40 with 5 MW left
        distribute(units, -20., default_eps_mw, new_inj, sat, &left);
        CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(40., 1e-6));
        CHECK_THAT(left[1], Catch::Matchers::WithinAbs(5., 1e-6));
        CHECK_THAT(left[0], Catch::Matchers::WithinAbs(0., 1e-12));
        // -30: shift -55 > 50, unit 1 leaves its max_p: nothing left
        distribute(units, -30., default_eps_mw, new_inj, sat, &left);
        CHECK(new_inj[1] < 40.);
        CHECK_THAT(left[1], Catch::Matchers::WithinAbs(0., 1e-12));
    }
    SECTION("chained then at once, both signs"){
        for(const real_type d2 : {-10., -25., 5.}){
            std::vector<Participant> chained = units;
            distribute(chained, -20., default_eps_mw, new_inj, sat, &left);
            for(std::size_t k = 0; k < chained.size(); ++k){
                chained[k].injection_mw = new_inj[k];
                chained[k].overshoot_mw = left[k];
            }
            distribute(chained, d2, default_eps_mw, new_inj, sat, &left);
            const std::vector<real_type> chained_inj = new_inj;
            distribute(units, -20. + d2, default_eps_mw, new_inj, sat, &left);
            CHECK_THAT(chained_inj[0], Catch::Matchers::WithinAbs(new_inj[0], 1e-6));
            CHECK_THAT(chained_inj[1], Catch::Matchers::WithinAbs(new_inj[1], 1e-6));
        }
    }
    SECTION("a unit a distribution saturates carries nothing"){
        std::vector<Participant> plain = {unit(0, 30., 1., 0., 500.), unit(1, 30., 1., 0., 40.)};
        distribute(plain, 30., default_eps_mw, new_inj, sat, &left);
        CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(40., 1e-12));
        CHECK_THAT(left[1], Catch::Matchers::WithinAbs(0., 1e-12));
    }
}

TEST_CASE("distribute with an overshoot: a drawing unit capped at 0 MW keeps its side", "[slack_redistribution]"){
    // a charging storage unit (range [-50, 50]) a positive mismatch pushed up to 0 MW, 5 MW
    // beyond it: the overshoot is above 0 MW (> 0), and 0 MW stays its UPPER bound
    std::vector<Participant> units = {unit(0, 30., 1., NaN, NaN), unit(1, 0., 1., -50., 50.)};
    units[1].kind = UnitKind::STORAGE;
    units[1].in_slack = false;
    units[1].overshoot_mw = 5.;
    std::vector<real_type> new_inj;
    std::vector<char> sat;
    SECTION("pushed down: it moves once the shift used up its overshoot"){
        // equal weights, shift -25: the slack unit gives 12.5, the storage unit 5 - 12.5
        const Report rep = distribute(units, -20., default_eps_mw, new_inj, sat);
        REQUIRE_FALSE(rep.all_saturated);
        CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(17.5, 1e-6));
        CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(-7.5, 1e-6));
        CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(0., 1e-6));
    }
    SECTION("pushed up: it stays at 0 MW, never crossing it"){
        const Report rep = distribute(units, 10., default_eps_mw, new_inj, sat);
        CHECK(sat[1] == 1);
        CHECK_THAT(new_inj[1], Catch::Matchers::WithinAbs(0., 1e-6));
        CHECK_THAT(new_inj[0], Catch::Matchers::WithinAbs(40., 1e-6));
        CHECK_THAT(rep.not_distributed_mw, Catch::Matchers::WithinAbs(0., 1e-6));
    }
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

TEST_CASE("a unit that can participate rejoins the slack when moved off its limit", "[slack_redistribution]"){
    ls2g::LSGrid grid = make_capped_grid(true);
    CHECK(grid.get_generators().get_can_participate_slack(0));   // a slack unit carries the flag
    grid.change_p_gen(1, 40.);                                   // still at its max_p
    CHECK_FALSE(grid.get_generators().is_slack(1));
    grid.change_p_gen(1, 35.);
    REQUIRE(grid.get_generators().is_slack(1));
    CHECK_THAT(grid.get_generators().get_slack_weight(1), Catch::Matchers::WithinAbs(0.5, 1e-15));
    SECTION("stranded: it stays in the slack, and takes no share until reconnected"){
        grid.deactivate_powerline(1);   // bus 2 (gen 1) and bus 3 cut off
        grid.consider_only_main_component(false);
        CHECK_FALSE(grid.get_generators().get_status()[1]);
        CHECK(grid.get_generators().is_slack(1));
        const ls2g::CplxVect V = grid.ac_pf(ls2g::CplxVect::Constant(5, ls2g::cplx_type(1., 0.)), 20, 1e-8);
        REQUIRE(V.size() > 0);
        CHECK(std::get<0>(grid.get_gen_res())(1) == 0.);
        grid.reactivate_powerline(1);
        grid.reactivate_gen(1);
        const ls2g::CplxVect V2 = grid.ac_pf(ls2g::CplxVect::Constant(5, ls2g::cplx_type(1., 0.)), 20, 1e-8);
        REQUIRE(V2.size() > 0);
        CHECK(std::abs(std::get<0>(grid.get_gen_res())(1) - 35.) > 1e-3);   // its share again
    }
    SECTION("removed by hand: it stays out"){
        grid.remove_gen_slackbus(1);
        CHECK_FALSE(grid.get_generators().get_can_participate_slack(1));
        grid.change_p_gen(1, 30.);
        CHECK_FALSE(grid.get_generators().is_slack(1));
    }
}

namespace {

// the AC then the DC solve of `cached` (which has solved before, and reuses what it
// built) must be those of `fresh` (which never solved): a slack change the cache missed
// would have `cached` solve with its former slack
void check_same_as_fresh(ls2g::LSGrid & cached, ls2g::LSGrid & fresh){
    const ls2g::CplxVect V0 = ls2g::CplxVect::Constant(5, ls2g::cplx_type(1., 0.));
    const ls2g::CplxVect Vc = cached.ac_pf(V0, 20, 1e-10);
    const ls2g::CplxVect Vf = fresh.ac_pf(V0, 20, 1e-10);
    REQUIRE(Vc.size() > 0);
    REQUIRE(Vc.size() == Vf.size());
    CHECK((Vc - Vf).cwiseAbs().maxCoeff() < 1e-12);
    const ls2g::CplxVect Vc_dc = cached.dc_pf(V0, 20, 1e-10);
    const ls2g::CplxVect Vf_dc = fresh.dc_pf(V0, 20, 1e-10);
    REQUIRE(Vc_dc.size() > 0);
    REQUIRE(Vc_dc.size() == Vf_dc.size());
    CHECK((Vc_dc - Vf_dc).cwiseAbs().maxCoeff() < 1e-12);
}

}  // namespace

TEST_CASE("a unit coming back to the slack invalidates the cached slack", "[slack_redistribution][cache_reuse]"){
    const ls2g::CplxVect V0 = ls2g::CplxVect::Constant(5, ls2g::cplx_type(1., 0.));
    ls2g::LSGrid grid = make_capped_grid(true);
    REQUIRE(grid.ac_pf(V0, 20, 1e-10).size() > 0);   // builds and marks the AC cache
    REQUIRE(grid.dc_pf(V0, 20, 1e-10).size() > 0);   // ... and the DC one
    REQUIRE_FALSE(grid.get_ac_algo_controler().has_slack_participate_changed());
    REQUIRE_FALSE(grid.get_dc_algo_controler().has_slack_participate_changed());

    grid.change_p_gen(1, 40.);   // still at its max_p: it stays out, the slack is unchanged
    CHECK_FALSE(grid.get_generators().is_slack(1));
    CHECK_FALSE(grid.get_ac_algo_controler().has_slack_participate_changed());
    CHECK_FALSE(grid.get_dc_algo_controler().has_slack_participate_changed());

    grid.change_p_gen(1, 35.);   // off it: back into the slack
    REQUIRE(grid.get_generators().is_slack(1));
    CHECK(grid.get_ac_algo_controler().has_slack_participate_changed());
    CHECK(grid.get_dc_algo_controler().has_slack_participate_changed());
    CHECK(grid.get_ac_algo_controler().has_slack_weight_changed());
    CHECK(grid.get_dc_algo_controler().has_slack_weight_changed());
    {
        ls2g::LSGrid fresh = make_capped_grid(true);
        fresh.change_p_gen(1, 35.);
        check_same_as_fresh(grid, fresh);
        grid.ac_pf(V0, 20, 1e-10);
        CHECK(std::abs(std::get<0>(grid.get_gen_res())(1) - 35.) > 1e-3);   // it takes a share
    }

    SECTION("stranded, then reconnected"){
        grid.deactivate_powerline(1);   // bus 2 (gen 1) and bus 3 cut off
        grid.consider_only_main_component(false);
        CHECK(grid.get_generators().is_slack(1));
        CHECK(grid.get_ac_algo_controler().has_slack_participate_changed());
        CHECK(grid.get_dc_algo_controler().has_slack_participate_changed());
        ls2g::LSGrid fresh = make_capped_grid(true);
        fresh.change_p_gen(1, 35.);
        fresh.deactivate_powerline(1);
        fresh.consider_only_main_component(false);
        check_same_as_fresh(grid, fresh);

        grid.reactivate_powerline(1);
        grid.reactivate_gen(1);
        CHECK(grid.get_ac_algo_controler().has_slack_participate_changed());
        CHECK(grid.get_dc_algo_controler().has_slack_participate_changed());
        ls2g::LSGrid fresh2 = make_capped_grid(true);
        fresh2.change_p_gen(1, 35.);
        fresh2.deactivate_powerline(1);
        fresh2.consider_only_main_component(false);
        fresh2.reactivate_powerline(1);
        fresh2.reactivate_gen(1);
        check_same_as_fresh(grid, fresh2);
    }
}
