// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// ---------------------------------------------------------------------------
// KIRCHHOFF'S CURRENT LAW, CHECKED ON THE PUBLISHED RESULTS
//
// Every other powerflow test in this suite checks the SOLVER: the voltages it
// returns, the mismatch it drove to zero, the iterations it took. This one
// checks what the library HANDS BACK. After a converged solve, the per-element
// results are the whole answer as far as a user is concerned -- p, q at every
// load, generator, shunt, storage, SVC, HVDC station and at both ends of every
// line and transformer -- and they only mean anything if they balance.
//
// So the sum here is deliberately naive. It uses no formula, no Ybus, no Sbus,
// no mismatch vector: it walks the containers, reads the numbers a python user
// would read through get_loads_res() / get_gen_res() / get_line_res1() and so
// on, adds them up per bus with their documented sign convention, and requires
// the total to be zero. Anything the solver knows and the results do not is
// exactly what this catches.
//
// The sign convention, taken from each container's own fillSbus (which is what
// defines what an element injects at its bus):
//   generators, static generators, SVCs,
//   HVDC converter stations              -> INJECT   (Sbus += p + i.q)
//   loads, storages, shunts              -> DRAW     (Sbus -= p + i.q)
//   line / transformer sides             -> DRAW, the flow LEAVING the bus into
//                                           the branch (res_p_side_1 is
//                                           V1 . conj(y11.V1 + y12.V2))
// ---------------------------------------------------------------------------

#include <complex>
#include <sstream>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "LSGrid.hpp"
#include "case_exotic_elements.hpp"

using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

/// Per-bus (grid numbering) complex power that does not balance, in MW / MVAr.
/// Zero everywhere iff the published results satisfy KCL.
CplxVect kcl_residual(const LSGrid & grid)
{
    const int nb_bus = static_cast<int>(grid.total_bus());
    CplxVect res = CplxVect::Zero(nb_bus);

    // ---- one-sided elements ------------------------------------------------
    // `sign` is +1 for an element that injects (generator convention) and -1 for
    // one that draws (load convention); see the note at the top of the file.
    auto add_one_sided = [&res](const auto & container, real_type sign) {
        const auto resu = container.get_res();
        const RealVect & p = std::get<0>(resu);
        const RealVect & q = std::get<1>(resu);
        const std::vector<bool> & status = container.get_status();
        const auto & buses = container.get_buses();
        const int nb_el = static_cast<int>(container.nb());
        for (int el_id = 0; el_id < nb_el; ++el_id) {
            if (!status[el_id]) continue;
            const int bus_id = buses(el_id).cast_int();
            if (bus_id < 0) continue;
            res(bus_id) += sign * cplx_type(p(el_id), q(el_id));
        }
    };

    add_one_sided(grid.get_generators(), 1.);
    add_one_sided(grid.get_static_generators(), 1.);
    add_one_sided(grid.get_svcs(), 1.);
    add_one_sided(grid.get_loads(), -1.);
    add_one_sided(grid.get_storages(), -1.);
    add_one_sided(grid.get_shunts(), -1.);

    // ---- two-sided elements ------------------------------------------------
    // Each side's (p, q) is the flow leaving that bus INTO the branch, so it is
    // subtracted at the bus it leaves. A side that is individually open carries
    // no flow (the container writes 0 there), but skip it explicitly rather than
    // relying on that -- an open side's bus id is not meaningful.
    auto add_two_sided = [&res](const auto & container) {
        const auto res1 = container.get_res_side_1();
        const auto res2 = container.get_res_side_2();
        const RealVect & p1 = std::get<0>(res1);
        const RealVect & q1 = std::get<1>(res1);
        const RealVect & p2 = std::get<0>(res2);
        const RealVect & q2 = std::get<1>(res2);
        const std::vector<bool> & status_global = container.get_status_global();
        const std::vector<bool> & status1 = container.get_status_side_1();
        const std::vector<bool> & status2 = container.get_status_side_2();
        const int nb_el = static_cast<int>(container.nb());
        for (int el_id = 0; el_id < nb_el; ++el_id) {
            if (!status_global[el_id]) continue;
            if (status1[el_id]) {
                const int bus_id = container.get_bus_side_1(el_id).cast_int();
                if (bus_id >= 0) res(bus_id) -= cplx_type(p1(el_id), q1(el_id));
            }
            if (status2[el_id]) {
                const int bus_id = container.get_bus_side_2(el_id).cast_int();
                if (bus_id >= 0) res(bus_id) -= cplx_type(p2(el_id), q2(el_id));
            }
        }
    };

    add_two_sided(grid.get_lines());
    add_two_sided(grid.get_trafos());
    // An HVDC line is two converter stations, and a station INJECTS in the
    // generator convention -- ConverterStationContainer::fillSbus_station stamps
    // `Sbus += target_p_mw_` with target_p_mw_ already signed (negative on the
    // sending side), and HvdcLineContainer::compute_results says so in as many
    // words for the droop case: "res_p is the station injection (generator
    // convention) = -flow". So these are added, like a generator, not subtracted.
    {
        const auto res1 = grid.get_dcline_res1();
        const auto res2 = grid.get_dcline_res2();
        const RealVect & p1 = std::get<0>(res1);
        const RealVect & q1 = std::get<1>(res1);
        const RealVect & p2 = std::get<0>(res2);
        const RealVect & q2 = std::get<1>(res2);
        const auto & hvdc = grid.get_dclines();
        const std::vector<bool> & status_global = hvdc.get_status_global();
        const std::vector<bool> & status1 = hvdc.get_status_side_1();
        const std::vector<bool> & status2 = hvdc.get_status_side_2();
        const int nb_el = static_cast<int>(hvdc.nb());
        for (int el_id = 0; el_id < nb_el; ++el_id) {
            if (!status_global[el_id]) continue;
            if (status1[el_id]) {
                const int bus_id = hvdc.get_bus_side_1(el_id).cast_int();
                if (bus_id >= 0) res(bus_id) += cplx_type(p1(el_id), q1(el_id));
            }
            if (status2[el_id]) {
                const int bus_id = hvdc.get_bus_side_2(el_id).cast_int();
                if (bus_id >= 0) res(bus_id) += cplx_type(p2(el_id), q2(el_id));
            }
        }
    }

    return res;
}

/// The exotic case14 with every exotic element switched off: plain IEEE14 (14
/// buses, 5 generators, 11 loads, 17 lines, 3 transformers, 1 shunt). Used to
/// validate the sum above on a grid whose balance is not in doubt.
LSGrid make_plain_case14()
{
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    grid.deactivate_svc(0);
    grid.deactivate_storage(0);
    for (int hvdc_id = 0; hvdc_id < 3; ++hvdc_id) grid.deactivate_dcline(hvdc_id);
    return grid;
}

/// solve, and report the worst bus
void solve(LSGrid & grid)
{
    const CplxVect v0 = CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()),
                                           cplx_type(1., 0.));
    const CplxVect V = grid.ac_pf(v0, 30, 1e-8);
    REQUIRE(V.size() == static_cast<Eigen::Index>(grid.total_bus()));
    REQUIRE(grid.get_algo().converged());
}

/// worst |residual| over the buses, plus a readable dump for the failure message
std::string describe(const CplxVect & residual, real_type tol)
{
    std::ostringstream out;
    for (Eigen::Index bus_id = 0; bus_id < residual.size(); ++bus_id) {
        if (std::abs(residual(bus_id)) <= tol) continue;
        out << "\n  bus " << bus_id
            << ": dP = " << residual(bus_id).real() << " MW"
            << ", dQ = " << residual(bus_id).imag() << " MVAr";
    }
    return out.str();
}

// The powerflow is solved to 1e-8 pu on a 100 MVA base, so a bus that balances
// exactly still shows ~1e-6 MW of rounding. Anything above this is a real
// imbalance, not numerical noise.
const real_type KCL_TOL = 1e-5;

}  // namespace


TEST_CASE("plain case14: the published results satisfy KCL at every bus",
          "[LSGrid][kcl]")
{
    // This is the calibration: a grid with nothing exotic on it, where the
    // per-bus sum of the published results MUST be zero. If it is not, the sum
    // above is wrong and every other case in this file is meaningless.
    LSGrid grid = make_plain_case14();
    solve(grid);

    const CplxVect residual = kcl_residual(grid);
    INFO("KCL residual on plain case14:" << describe(residual, KCL_TOL));
    for (Eigen::Index bus_id = 0; bus_id < residual.size(); ++bus_id) {
        CHECK(std::abs(residual(bus_id)) < KCL_TOL);
    }
}

TEST_CASE("case14 + storage: the published results satisfy KCL at every bus",
          "[LSGrid][kcl]")
{
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    grid.deactivate_svc(0);
    for (int hvdc_id = 0; hvdc_id < 3; ++hvdc_id) grid.deactivate_dcline(hvdc_id);
    solve(grid);

    const CplxVect residual = kcl_residual(grid);
    INFO("KCL residual on case14 + storage:" << describe(residual, KCL_TOL));
    for (Eigen::Index bus_id = 0; bus_id < residual.size(); ++bus_id) {
        CHECK(std::abs(residual(bus_id)) < KCL_TOL);
    }
}

TEST_CASE("case14 + voltage-mode SVC: the published results satisfy KCL at every bus",
          "[LSGrid][kcl]")
{
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    grid.deactivate_storage(0);
    for (int hvdc_id = 0; hvdc_id < 3; ++hvdc_id) grid.deactivate_dcline(hvdc_id);
    solve(grid);

    const CplxVect residual = kcl_residual(grid);
    INFO("KCL residual on case14 + SVC:" << describe(residual, KCL_TOL));
    for (Eigen::Index bus_id = 0; bus_id < residual.size(); ++bus_id) {
        CHECK(std::abs(residual(bus_id)) < KCL_TOL);
    }
}

TEST_CASE("case14 + HVDC lines: the published results satisfy KCL at every bus",
          "[LSGrid][kcl]")
{
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    grid.deactivate_svc(0);
    grid.deactivate_storage(0);
    solve(grid);

    const CplxVect residual = kcl_residual(grid);
    INFO("KCL residual on case14 + HVDC:" << describe(residual, KCL_TOL));
    for (Eigen::Index bus_id = 0; bus_id < residual.size(); ++bus_id) {
        CHECK(std::abs(residual(bus_id)) < KCL_TOL);
    }
}

TEST_CASE("exotic case14: the published results satisfy KCL at every bus",
          "[LSGrid][kcl]")
{
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    solve(grid);

    const CplxVect residual = kcl_residual(grid);
    INFO("KCL residual on the full exotic case14:" << describe(residual, KCL_TOL));
    for (Eigen::Index bus_id = 0; bus_id < residual.size(); ++bus_id) {
        CHECK(std::abs(residual(bus_id)) < KCL_TOL);
    }
}

// ---------------------------------------------------------------------------
// The cases above all pass -- but none of them can FAIL for the reason that
// matters, and saying so is the point of this block.
//
// KCL is a per-bus law. It sees the TOTAL at a bus, never the split between the
// elements sitting on it. If a bus carries a generator and something whose
// injection the Newton solves for itself (a voltage-mode SVC, a remote
// controller, an angle-droop HVDC station), and the generator is credited with
// the other one's power, the bus still balances to the last digit and this
// check stays green.
//
// In the exotic fixture no bus carries both: the generators are at buses 0, 6,
// 7, 10 and 12, the HVDC stations at 1, 2, 3, 4, 5 and 13, the SVC at 5, the
// storage at 4. So the cases above establish that the published results are
// self-consistent, and nothing more. The three below move one element so that a
// generator and an NR-solved injection DO share a bus -- which is the only
// arrangement in which the attribution can be caught being wrong.
// ---------------------------------------------------------------------------

TEST_CASE("KCL holds with a generator on the same bus as an angle-droop HVDC station",
          "[LSGrid][kcl]")
{
    // hvdc line 1 has the droop enabled and sits on buses 1 and 2; in AC its
    // active power is NOT stamped into Sbus (HvdcLineContainer::fillSbus skips
    // it, the Hvdc extension of the NR system carries it), so the raw per-bus
    // mismatch compute_results() derives from `V .* conj(Ybus . V) - Sbus`
    // carries that power as if it were unexplained. Put a generator there and
    // the question becomes whether it gets credited with it.
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    grid.deactivate_svc(0);
    grid.deactivate_storage(0);
    grid.change_bus_gen_python(1, 1);   // generator 1: bus 6 -> bus 1
    solve(grid);

    // This used to publish NaN for the station and -0 for the generator. The
    // proportional split divided (max_q_me - min_q_me + eps) by
    // (max_q_bus - min_q_bus + n.eps), and the converter stations of hvdc lines 1
    // and 2 carry +/- DBL_MAX as their reactive limits (what the pypowsybl
    // converter writes for "unbounded", and this fixture is a capture of a real
    // conversion) -- so both spans overflowed to +inf and the ratio was inf / inf,
    // while the generator, measured against that same infinite denominator, got
    // 90 / inf = 0. The bus balanced in P and, in Q, only in the sense that NaN
    // absorbs everything.
    //
    // GenericContainer::_q_share now falls back to an equal split when a range is
    // not finite -- the same answer the formula already gives when every range is
    // zero, and the only defensible one when there is no information to weight two
    // unbounded machines against each other.
    const CplxVect residual = kcl_residual(grid);
    INFO("KCL residual, generator on a droop-HVDC bus:" << describe(residual, KCL_TOL));
    for (Eigen::Index bus_id = 0; bus_id < residual.size(); ++bus_id) {
        CHECK(std::abs(residual(bus_id)) < KCL_TOL);
    }
}

TEST_CASE("a generator cannot share the SVC's regulated bus (v1), so that path cannot be probed",
          "[LSGrid][kcl]")
{
    // The obvious third arrangement -- a generator on the bus a voltage-mode SVC
    // regulates -- cannot be built at all: v1 requires an SVC to be the only
    // controller of its bus and says so. Recorded here so the gap is explicit:
    // whether the generator would be credited with the SVC's reactive output is
    // a question this suite CANNOT ask yet, not one it has answered.
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    grid.deactivate_storage(0);
    for (int hvdc_id = 0; hvdc_id < 3; ++hvdc_id) grid.deactivate_dcline(hvdc_id);
    grid.change_bus_gen_python(1, 5);   // generator 1: bus 6 -> bus 5, the SVC's bus
    grid.change_v_gen(1, 1.0);          // agree with the SVC's setpoint

    const CplxVect v0 = CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()),
                                           cplx_type(1., 0.));
    CHECK_THROWS_WITH(grid.ac_pf(v0, 30, 1e-8),
                      Catch::Matchers::ContainsSubstring("must be the only controller of its bus"));
}


// KNOWN DEFECT, and a different one from the Q attribution above: this
// configuration does not solve AT ALL, so nothing here reaches the write-back --
// the REQUIRE in solve() is what fails.
//
// The trigger is narrow and has nothing to do with the droop: it is a
// VOLTAGE-REGULATING converter station sharing a bus with the slack generator.
// All four such buses of this fixture reproduce it (1, 2, 3 and 13, the stations
// of hvdc lines 1 and 2); the stations of line 0, whose voltage_regulator_on is
// false, do not (buses 4 and 5), and neither does a bus with no station at all.
// It is not a setpoint conflict -- making the two agree changes nothing -- not
// solver-specific (KLU, SparseLU and the single-slack Newton all fail), and not
// a starting-point problem: it diverges warm-started from the exact solution of
// the unmodified grid.
//
// What happens is that every residual falls to ~1e-7 by the third iteration
// except the REACTIVE balance at that bus, which then grows geometrically until
// the iterate runs away (voltages at 1e-4 pu and at 8.7 pu by iteration 30).
// A station with voltage_regulator_on excludes its reactive injection from Sbus
// (ConverterStationContainer::fillSbus_station), on the PV-bus contract: the bus
// owns a free Vm, and the reactive output is recovered from the residual
// afterwards. A slack bus pinned by a local generator honours neither half --
// LSGrid::get_free_vm_slack_solver_buses grants a free Vm and a Q equation only
// to a slack bus NO local generator pins, and it looks at generators alone, never
// at converter stations. So nothing in the assembled system determines that
// station's reactive injection, and the Newton is left with a residual it has no
// unknown to move.
//
// The SVC path rejects exactly this class of configuration with a clear message
// ("is at a bus with no reactive (Q) equation ... not supported in v1"); the
// converter-station path has no such guard and simply fails to converge.
//
// [!shouldfail] keeps the suite green while it is open and turns this test red
// the day it is fixed, which is the signal to drop the tag.
TEST_CASE("KCL holds with the slack generator on an angle-droop HVDC bus",
          "[LSGrid][kcl][!shouldfail]")
{
    // the sharpest arrangement: generator 0 IS the slack, so
    // GeneratorContainer::set_p_slack hands it the whole active mismatch of its
    // bus -- and that bus now also carries the droop station whose power never
    // reached Sbus.
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    grid.deactivate_svc(0);
    grid.deactivate_storage(0);
    grid.change_bus_gen_python(0, 1);   // slack generator: bus 0 -> bus 1
    solve(grid);

    const CplxVect residual = kcl_residual(grid);
    INFO("KCL residual, slack generator on a droop-HVDC bus:" << describe(residual, KCL_TOL));
    for (Eigen::Index bus_id = 0; bus_id < residual.size(); ++bus_id) {
        CHECK(std::abs(residual(bus_id)) < KCL_TOL);
    }
}
