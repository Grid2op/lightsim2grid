// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// LSGrid::set_dc_distribute_slack_on_can_participate: a DC powerflow distributing the imbalance
// on the units flagged "can participate in the slack" rather than on the grid's own slack
// participants (OpenLoadFlow's DC_VALUES start around a single-slack Newton). C++14 only.

#include <cmath>
#include <tuple>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"

using Catch::Approx;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::real_type;

namespace {

// buses 0-1-2-3 in a line, a 200 MW load at bus 3, generators at buses 0 (the slack), 1 and
// 2 producing 150 MW in total: a 50 MW imbalance in DC (no losses)
LSGrid make_grid()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(4, 1, RealVect::Constant(4, 138.), 0, 0);
    Eigen::VectorXi from_id(3), to_id(3);
    from_id << 0, 1, 2;
    to_id << 1, 2, 3;
    grid.init_powerlines(RealVect::Constant(3, 0.01), RealVect::Constant(3, 0.05),
                         CplxVect::Zero(3), from_id, to_id);
    RealVect load_p(1), load_q(1);
    load_p << 200.;
    load_q << 20.;
    Eigen::VectorXi load_bus(1);
    load_bus << 3;
    grid.init_loads(load_p, load_q, load_bus);
    RealVect gen_p(3), gen_v(3), gen_min_q(3), gen_max_q(3);
    gen_p << 50., 40., 60.;
    gen_v << 1.02, 1.01, 1.01;
    gen_min_q << -500., -500., -500.;
    gen_max_q << 500., 500., 500.;
    Eigen::VectorXi gen_bus(3);
    gen_bus << 0, 1, 2;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    return grid;
}

CplxVect dc(LSGrid & grid)
{
    return grid.dc_pf(CplxVect::Constant(grid.total_bus(), ls2g::cplx_type(1., 0.)), 10, 1e-10);
}

RealVect gen_p(const LSGrid & grid)
{
    const auto res_p = std::get<0>(grid.get_generators().get_res());
    RealVect res(res_p.size());
    for (int i = 0; i < res.size(); ++i) res(i) = static_cast<real_type>(res_p(i));
    return res;
}

}  // namespace

TEST_CASE("off by default: the DC keeps the grid's own slack", "[dc][slack]")
{
    LSGrid grid = make_grid();
    RealVect w(3);
    w << 1., 2., 3.;
    grid.set_gen_can_participate_slack({true, true, true}, w);
    CHECK_FALSE(grid.get_dc_distribute_slack_on_can_participate());
    REQUIRE(dc(grid).size() == 4);
    const RealVect p = gen_p(grid);
    CHECK(p(0) == Approx(100.));  // the slack generator took the whole imbalance
    CHECK(p(1) == Approx(40.));
    CHECK(p(2) == Approx(60.));
}

TEST_CASE("on: the units flagged to share the slack do, with their weights", "[dc][slack]")
{
    // the reference: the same units as the grid's own distributed slack
    LSGrid ref = make_grid();
    ref.add_gen_slackbus(0, 1.);
    ref.add_gen_slackbus(1, 2.);
    ref.add_gen_slackbus(2, 3.);
    const CplxVect V_ref = dc(ref);
    REQUIRE(V_ref.size() == 4);

    LSGrid grid = make_grid();
    RealVect w(3);
    w << 1., 2., 3.;
    grid.set_gen_can_participate_slack({true, true, true}, w);
    REQUIRE(dc(grid).size() == 4);  // a first solve with the option off
    grid.set_dc_distribute_slack_on_can_participate(true);
    const CplxVect V = dc(grid);
    REQUIRE(V.size() == 4);
    for (int i = 0; i < 4; ++i) CHECK(std::arg(V(i)) == Approx(std::arg(V_ref(i))).margin(1e-12));
    const RealVect p = gen_p(grid);
    CHECK(p(0) == Approx(50. + 50. / 6.));
    CHECK(p(1) == Approx(40. + 100. / 6.));
    CHECK(p(2) == Approx(60. + 150. / 6.));

    // new weights reach the next DC solve
    w << 0., 1., 1.;
    grid.set_gen_can_participate_slack({false, true, true}, w);
    REQUIRE(dc(grid).size() == 4);
    const RealVect p2 = gen_p(grid);
    CHECK(p2(0) == Approx(50.));
    CHECK(p2(1) == Approx(65.));
    CHECK(p2(2) == Approx(85.));

    // and switching it off brings the grid's own slack back
    grid.set_dc_distribute_slack_on_can_participate(false);
    REQUIRE(dc(grid).size() == 4);
    CHECK(gen_p(grid)(0) == Approx(100.));
}

TEST_CASE("on, without any unit flagged: the grid's own slack", "[dc][slack]")
{
    LSGrid grid = make_grid();
    grid.set_dc_distribute_slack_on_can_participate(true);
    REQUIRE(dc(grid).size() == 4);
    CHECK(gen_p(grid)(0) == Approx(100.));
}

TEST_CASE("a copied grid keeps the option", "[dc][slack]")
{
    LSGrid grid = make_grid();
    grid.set_dc_distribute_slack_on_can_participate(true);
    LSGrid copy(grid);
    CHECK(copy.get_dc_distribute_slack_on_can_participate());
}
