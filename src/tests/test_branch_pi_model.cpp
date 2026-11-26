// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// The pi model of a branch computed for other parameters without changing it
// (LineContainer / TrafoContainer::pi_model_with) and the stamped block of a branch open at
// one end (BranchContainer::kron_block): what a batch row that changes branch parameters
// stamps. C++14 only.

#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"

using ls2g::BranchContainer;
using ls2g::BranchPiModel;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

const real_type NaN = std::numeric_limits<real_type>::quiet_NaN();
const cplx_type cNaN(NaN, NaN);

// buses 0-1-2: a line 0-1, a transformer 1-2, a load at bus 2, the slack at bus 0
LSGrid make_grid()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(3, 1, RealVect::Constant(3, 138.), 0, 0);
    Eigen::VectorXi line_from(1), line_to(1);
    line_from << 0;
    line_to << 1;
    CplxVect line_h(1);
    line_h << cplx_type(0.001, 0.02);
    grid.init_powerlines_full(RealVect::Constant(1, 0.01), RealVect::Constant(1, 0.1), line_h, 0.5 * line_h,
                              line_from, line_to);
    RealVect trafo_r(1), trafo_x(1), trafo_ratio(1), trafo_shift(1);
    CplxVect trafo_b(1);
    trafo_r << 0.01;
    trafo_x << 0.2;
    trafo_b << cplx_type(0.002, 0.04);
    trafo_ratio << 1.05;
    trafo_shift << 0.;
    Eigen::VectorXi trafo_bus1(1), trafo_bus2(1);
    trafo_bus1 << 1;
    trafo_bus2 << 2;
    grid.init_trafo(trafo_r, trafo_x, trafo_b, trafo_ratio, trafo_shift, {false}, trafo_bus1, trafo_bus2, false);
    Eigen::VectorXi load_bus(1);
    load_bus << 2;
    grid.init_loads(RealVect::Constant(1, 50.), RealVect::Constant(1, 20.), load_bus);
    Eigen::VectorXi gen_bus(1);
    gen_bus << 0;
    grid.init_generators(RealVect::Zero(1), RealVect::Constant(1, 1.02), RealVect::Constant(1, -1000.),
                         RealVect::Constant(1, 1000.), gen_bus);
    grid.add_gen_slackbus(0, 1.);
    return grid;
}

// a ratio tap changer (the transformer at position 1), a phase one, and a shift-dependent
// impedance: every correction the neutral values go through
void add_corrections(LSGrid & grid)
{
    grid.set_trafo_ratio_tap_changer(0, 0, 1, {0.95, 1.05, 1.1}, {0., 10., 20.}, {0., 5., 10.},
                                     {0., 2., 4.}, {0., 50., 100.});
    grid.set_trafo_shift_dependent_rx(true, {{-10. / 180. * M_PI, 0., 10. / 180. * M_PI}}, {{-20., 0., 20.}});
    grid.set_trafo_phase_tap_changer(0, -1, 0, {1., 1., 1.}, {-10., 0., 10.}, {-20., 0., 20.}, {-20., 0., 20.},
                                     {0., 0., 0.}, {0., 0., 0.});
}

// the model `branch` holds for element `el`, read back off its coefficients
BranchPiModel model_of(const BranchContainer & branch, int el, real_type dc_x_tau_shift)
{
    BranchPiModel res;
    res.yac = {branch.yac_11()(el), branch.yac_12()(el), branch.yac_21()(el), branch.yac_22()(el)};
    res.ydc = {branch.ydc_11()(el), branch.ydc_12()(el), branch.ydc_21()(el), branch.ydc_22()(el)};
    res.dc_x_tau_shift = dc_x_tau_shift;
    return res;
}

// bit for bit
bool same(const BranchPiModel & a, const BranchPiModel & b)
{
    for (int k = 0; k < 4; ++k) {
        if (a.yac[k] != b.yac[k] || a.ydc[k] != b.ydc[k]) return false;
    }
    return a.dc_x_tau_shift == b.dc_x_tau_shift;
}

}  // namespace

TEST_CASE("a line's pi model for other parameters is the one it gets with them", "[branch_pi_model]")
{
    LSGrid grid = make_grid();
    const auto & lines = grid.get_powerlines_as_data();
    // nothing given: exactly its own
    CHECK(same(lines.pi_model_with(0, NaN, NaN, cNaN, cNaN), model_of(lines, 0, 0.)));

    // given: what update_powerlines_parameters gives it; a NaN keeps the line's value
    const BranchPiModel model = lines.pi_model_with(0, 0.02, NaN, cplx_type(0., 0.03), cplx_type(NaN, 1.));
    RealVect r(1), x(1);
    CplxVect h1(1), h2(1);
    r << 0.02;
    x << 0.1;
    h1 << cplx_type(0., 0.03);
    h2 << cplx_type(0.0005, 0.01);
    grid.update_powerlines_parameters(r, x, h1, h2);
    CHECK(same(model, model_of(grid.get_powerlines_as_data(), 0, 0.)));

    CHECK_THROWS(lines.pi_model_with(1, NaN, NaN, cNaN, cNaN));
}

TEST_CASE("a transformer's pi model for other parameters is the one it gets with them", "[branch_pi_model]")
{
    for (bool corrected : {false, true}) {
        LSGrid grid = make_grid();
        if (corrected) add_corrections(grid);
        const auto & trafos = grid.get_trafos_as_data();
        CHECK(same(trafos.pi_model_with(0, NaN, NaN, cNaN, cNaN, NaN, NaN),
                   model_of(trafos, 0, trafos.dc_x_tau_shift()(0))));

        // neutral r, x, h as update_trafos_parameters_full takes them, the ratio and the
        // shift as change_ratio_trafo / change_shift_trafo do
        const real_type shift = 5. / 180. * M_PI;
        const BranchPiModel model = trafos.pi_model_with(0, 0.015, 0.25, cplx_type(0.001, 0.01), cplx_type(0., 0.02),
                                                         1.02, shift);
        RealVect r(1), x(1);
        CplxVect h1(1), h2(1);
        r << 0.015;
        x << 0.25;
        h1 << cplx_type(0.001, 0.01);
        h2 << cplx_type(0., 0.02);
        grid.update_trafos_parameters_full(r, x, h1, h2);
        grid.change_ratio_trafo(0, 1.02);
        grid.change_shift_trafo(0, shift);
        const auto & after = grid.get_trafos_as_data();
        CHECK(same(model, model_of(after, 0, after.dc_x_tau_shift()(0))));
        CHECK(model.dc_x_tau_shift != 0.);
    }
}

TEST_CASE("kron_block is what a branch open at one end stamps", "[branch_pi_model]")
{
    const std::array<cplx_type, 4> y = {cplx_type(1., -10.), cplx_type(-1., 9.), cplx_type(-1.1, 9.5), cplx_type(1.2, -11.)};
    const cplx_type zero(0., 0.);
    CHECK(BranchContainer::kron_block(y, true, true, true) == y);
    for (bool s1 : {false, true}) {
        for (bool s2 : {false, true}) {
            CHECK(BranchContainer::kron_block(y, false, s1, s2) == std::array<cplx_type, 4>{zero, zero, zero, zero});
        }
    }
    CHECK(BranchContainer::kron_block(y, true, false, false) == std::array<cplx_type, 4>{zero, zero, zero, zero});
    // the reduced end (to the last bits: the library and this file may contract the
    // arithmetic differently)
    const auto close = [](cplx_type a, cplx_type b){ return std::abs(a - b) <= 1e-14 * std::abs(b); };
    const std::array<cplx_type, 4> open1 = BranchContainer::kron_block(y, true, false, true);
    CHECK((open1[0] == zero && open1[1] == zero && open1[2] == zero));
    CHECK(close(open1[3], y[3] - y[2] * y[1] / y[0]));
    const std::array<cplx_type, 4> open2 = BranchContainer::kron_block(y, true, true, false);
    CHECK((open2[1] == zero && open2[2] == zero && open2[3] == zero));
    CHECK(close(open2[0], y[0] - y[2] * y[1] / y[3]));
}
