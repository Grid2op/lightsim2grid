// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// `compute_physical_violations`: the batch-side checks whose violation means the row is not
// a state the grid can reach -- the reactive capability of a bus (batch_algorithm/
// BusQCheck.hpp) and the active power of an angle-droop hvdc line (HvdcPCheck.hpp).
//
// Two things are tested, and they are different:
//
//  * the VALUE. The reactive power a bus had to produce is re-derived from the algorithm's
//    per-bus mismatch and its controller list, because a batch publishes voltages and
//    nothing else. Every test therefore pins it against the one thing that cannot be
//    wrong: the sum, over that bus' machines, of what a single-shot LSGrid::ac_pf
//    publishes (GeneratorContainer's res_q_, through LSGrid::get_gen_res).
//  * the QUESTION. It is asked per BUS, against the SUM of its machines' limits -- not per
//    machine. "a bus is fine while one of its machines is not" below is the test that
//    tells the two apart, and it is the reason this is not a per-generator check.
//
// Note that a generator's reactive limits have NO influence on the solution (nothing
// enforces them), so a bus' reactive power can be measured on one limit configuration and
// the limits chosen afterwards to land either side of it.

#include <cmath>
#include <complex>
#include <limits>
#include <tuple>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::ContingencyAnalysis;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::LimitViolation;
using ls2g::LimitViolationType;
using ls2g::RealVect;
using ls2g::ScenarioSweep;
using ls2g::SvcContainer;
using ls2g::TimeSeries;
using ls2g::ViolationCategory;
using ls2g::ViolationElementType;
using ls2g::cplx_type;
using ls2g::real_type;
using ls2g::violation_category;

namespace {

using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using BoolMat = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const real_type V_SET = 1.05;
const real_type LOAD_P = 80.;
const real_type LOAD_Q = 60.;
const int NB_BUS = 4;
const real_type WIDE_Q = 1000.;
const int GEN_BUS = 1;
const int SVC_BUS = NB_BUS - 1;  // the load bus

struct GenSpec
{
    int bus;
    real_type vset;
    real_type p;
    real_type min_q;
    real_type max_q;
    int regulated_bus = -1;  // -1: regulates its own bus (the classical PV path)
};

// The 4-bus radial feeder 0-1-2-3 (r = 0.01, x = 0.1 pu, sn_mva = 100) shared by
// test_batch_voltage_control.cpp and test_scenario_sweep_violations.cpp, with the load on
// bus 3. `meshed` adds a fourth line 0--3, so that a single line outage leaves the load
// connected and the contingency rows actually converge. Generator 0 is the slack.
// `station_q_on_gen_bus`, when finite, also puts an hvdc VSC converter station on bus 1
// regulating the bus it stands on, with reactive limits +/- that value: bus 1's reactive
// power is then produced by the generators AND by the station, and so is its capability.
LSGrid make_grid(const std::vector<GenSpec> & gens, bool meshed = false,
                 real_type station_q_on_gen_bus = std::numeric_limits<real_type>::quiet_NaN())
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    const RealVect bus_vn_kv = RealVect::Constant(NB_BUS, 138.);
    grid.init_bus(static_cast<unsigned int>(NB_BUS), 1, bus_vn_kv, 0, 0);

    const int n_line = meshed ? NB_BUS : NB_BUS - 1;
    const RealVect branch_r = RealVect::Constant(n_line, 0.01);
    const RealVect branch_x = RealVect::Constant(n_line, 0.1);
    const CplxVect branch_h = CplxVect::Zero(n_line);
    Eigen::VectorXi from_id(n_line), to_id(n_line);
    for (int i = 0; i < NB_BUS - 1; ++i) { from_id(i) = i; to_id(i) = i + 1; }
    if (meshed) { from_id(NB_BUS - 1) = 0; to_id(NB_BUS - 1) = NB_BUS - 1; }
    grid.init_powerlines(branch_r, branch_x, branch_h, from_id, to_id);

    RealVect load_p(1), load_q(1);
    load_p << LOAD_P;
    load_q << LOAD_Q;
    Eigen::VectorXi load_bus(1);
    load_bus << NB_BUS - 1;
    grid.init_loads(load_p, load_q, load_bus);

    const int nb_gen = static_cast<int>(gens.size());
    RealVect gen_p(nb_gen), gen_v(nb_gen), gen_min_q(nb_gen), gen_max_q(nb_gen);
    Eigen::VectorXi gen_bus(nb_gen);
    for (int k = 0; k < nb_gen; ++k) {
        gen_p(k) = gens[k].p;
        gen_v(k) = gens[k].vset;
        gen_min_q(k) = gens[k].min_q;
        gen_max_q(k) = gens[k].max_q;
        gen_bus(k) = gens[k].bus;
    }
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    for (int k = 0; k < nb_gen; ++k) {
        if (gens[k].regulated_bus >= 0) grid.set_gen_regulated_bus(k, gens[k].regulated_bus);
    }
    if (std::isfinite(station_q_on_gen_bus)) {
        // a VSC station on bus 1 regulating bus 1 (the ordinary PV path, like a local
        // generator: it shares the bus' reactive residual -- see
        // LSGrid::_collect_q_residual_shares), its other end on bus 2 in fixed-Q mode
        Eigen::VectorXi bus1(1), bus2(1);
        bus1 << GEN_BUS;
        bus2 << GEN_BUS + 1;
        const std::vector<int> type1{0}, type2{0};
        const std::vector<bool> vreg1{true}, vreg2{false};
        RealVect loss1(1), loss2(1), vm1(1), vm2(1), q1(1), q2(1);
        RealVect minq1(1), maxq1(1), minq2(1), maxq2(1), pf1(1), pf2(1);
        RealVect p_set(1), r_ohm(1), vnom(1), droop_p0(1), droop_slope(1), pmax12(1), pmax21(1);
        loss1 << 0.; loss2 << 0.;
        vm1 << V_SET; vm2 << 1.;
        q1 << 0.; q2 << 0.;
        minq1 << -station_q_on_gen_bus; maxq1 << station_q_on_gen_bus;
        minq2 << -WIDE_Q; maxq2 << WIDE_Q;
        pf1 << 1.; pf2 << 1.;
        p_set << 1.;
        r_ohm << 1.;
        vnom << 100.;
        droop_p0 << 0.; droop_slope << 0.;
        pmax12 << 20.; pmax21 << 20.;
        const std::vector<int> conv_mode{0};
        const std::vector<bool> droop_enabled{false};
        grid.init_hvdc_lines(bus1, bus2, type1, type2, loss1, loss2, vreg1, vreg2, vm1, vm2,
                             q1, q2, minq1, maxq1, minq2, maxq2, pf1, pf2, conv_mode, p_set,
                             r_ohm, vnom, droop_enabled, droop_p0, droop_slope, pmax12, pmax21);
    }
    grid.tell_solver_need_reset();
    return grid;
}

// the slack generator, with limits wide enough never to be reported
GenSpec slack_gen() { return GenSpec{0, 1.02, 0., -WIDE_Q, WIDE_Q, -1}; }

// the same feeder, with a voltage-mode SVC holding the LOAD bus at V_SET and the susceptance
// range `[b_min_pu, b_max_pu]` (pu, base sn_mva). It is the only element regulating that bus
// -- LSGrid refuses an SVC sharing its regulated bus with another controller -- so the bus'
// whole reactive power is the SVC's own.
//
// A positive susceptance is CAPACITIVE and produces reactive power: a shunt `jb` consumes
// `-j.b.|V|^2`, ie injects `+b.|V|^2` in generator convention. So `[0, b]` is a machine that
// can only produce and `[-b, 0]` one that can only absorb -- which is what the asymmetric
// section of the SVC test below pins down, and the symmetric case cannot.
LSGrid make_svc_grid_asym(real_type b_min_pu, real_type b_max_pu);
LSGrid make_svc_grid(real_type b_pu) { return make_svc_grid_asym(-b_pu, b_pu); }

LSGrid make_svc_grid_asym(real_type b_min_pu, real_type b_max_pu)
{
    LSGrid grid = make_grid(std::vector<GenSpec>{slack_gen()});
    std::vector<int> mode{SvcContainer::RegulationMode::VOLTAGE};
    RealVect target_vm(1), q_set(1), slope(1), b_min(1), b_max(1);
    target_vm << V_SET;
    q_set << 0.;
    slope << 0.;
    b_min << b_min_pu;
    b_max << b_max_pu;
    Eigen::VectorXi reg_bus(1), svc_bus(1);
    reg_bus << SVC_BUS;
    svc_bus << SVC_BUS;
    grid.init_svcs(mode, target_vm, q_set, slope, b_min, b_max, reg_bus, svc_bus);
    grid.add_gen_slackbus(0, 1.);
    grid.tell_solver_need_reset();
    return grid;
}

// The same feeder plus an ANGLE-DROOP ("AC emulation") hvdc line from bus 1 to bus 3, in
// parallel with the 1-2-3 path: `p = p0_mw + k_mw_per_rad . (theta_1 - theta_3)`, with no
// converter or dc-line losses (lf = 0, r = 0) so that the flow leaving bus 1 into the hvdc
// is exactly that `p` and the reference below can be read straight off `res_p1_mw`.
//
// `nb_spare_bus` prepends buses that carry NOTHING: they are therefore outside the solved
// system, every element's grid bus id is shifted up by that many, and the solver ids of the
// live buses no longer equal their grid ids -- the labelling mistake this check could make
// silently (reading one bus' angle as another's).
// `reverse_ends` swaps which end of the line is side 1: the same physical flow (bus 1 ->
// bus 3, pulled by the load) is then described as leaving SIDE 2, and bounded by
// pmax_2to1 instead. Note `init_hvdc_lines` takes the droop slope in MW per DEGREE (it
// stores MW/rad, see HvdcLineContainer::init).
LSGrid make_droop_grid(real_type pmax_1to2, real_type pmax_2to1,
                       real_type p0_mw = 30., real_type k_mw_per_deg = 400.,
                       int nb_spare_bus = 0, bool reverse_ends = false)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    const int off = nb_spare_bus;  // grid bus id of what would otherwise be bus 0
    const int nb_bus = NB_BUS + off;
    grid.init_bus(static_cast<unsigned int>(nb_bus), 1, RealVect::Constant(nb_bus, 138.), 0, 0);
    const int n_line = NB_BUS - 1;
    Eigen::VectorXi from_id(n_line), to_id(n_line);
    for (int i = 0; i < n_line; ++i) { from_id(i) = off + i; to_id(i) = off + i + 1; }
    grid.init_powerlines(RealVect::Constant(n_line, 0.01), RealVect::Constant(n_line, 0.1),
                         CplxVect::Zero(n_line), from_id, to_id);

    RealVect load_p(1), load_q(1);
    load_p << LOAD_P;
    load_q << LOAD_Q;
    Eigen::VectorXi load_bus(1);
    load_bus << off + NB_BUS - 1;
    grid.init_loads(load_p, load_q, load_bus);

    RealVect gen_p(1), gen_v(1), gen_min_q(1), gen_max_q(1);
    Eigen::VectorXi gen_bus(1);
    gen_p << 0.;
    gen_v << 1.02;
    gen_min_q << -WIDE_Q;
    gen_max_q << WIDE_Q;
    gen_bus << off;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);

    Eigen::VectorXi bus1(1), bus2(1);
    const int hvdc_end_a = off + GEN_BUS;        // 1 (+ the spare offset)
    const int hvdc_end_b = off + NB_BUS - 1;     // 3, the load bus
    bus1 << (reverse_ends ? hvdc_end_b : hvdc_end_a);
    bus2 << (reverse_ends ? hvdc_end_a : hvdc_end_b);
    const std::vector<int> type1{0}, type2{0}, mode{0};       // VSC, SIDE_1_RECTIFIER
    const std::vector<bool> vreg1{false}, vreg2{false}, droop_on{true};
    const RealVect zero = RealVect::Zero(1);
    RealVect vm(1), q_lim(1), pf(1), p0(1), k(1), pmax12(1), pmax21(1);
    vm << 1.;
    q_lim << WIDE_Q;
    pf << 1.;
    p0 << p0_mw;
    k << k_mw_per_deg;
    pmax12 << pmax_1to2;
    pmax21 << pmax_2to1;
    grid.init_hvdc_lines(bus1, bus2, type1, type2, zero, zero, vreg1, vreg2, vm, vm,
                         zero, zero, -q_lim, q_lim, -q_lim, q_lim, pf, pf, mode, zero,
                         zero, zero, droop_on, p0, k, pmax12, pmax21);
    grid.tell_solver_need_reset();
    return grid;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

// every generator's reactive output as a single-shot ac_pf publishes it: the reference this
// whole file is written against
RealVect reference_gen_q(LSGrid & grid)
{
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    grid.ac_pf(flat_start(grid), 30, 1e-11);
    RealVect res(std::get<1>(grid.get_gen_res()));
    return res;
}

// ... summed over the generators of one bus: what this check reports
real_type reference_bus_q(const std::vector<GenSpec> & gens, int bus, bool meshed = false)
{
    LSGrid grid = make_grid(gens, meshed);
    const RealVect q = reference_gen_q(grid);
    real_type total = 0.;
    for (int k = 0; k < static_cast<int>(gens.size()); ++k) {
        if (gens[k].bus == bus) total += q(k);
    }
    return total;
}

// the violation reported on `bus_id`, or nullptr
const LimitViolation * find_bus(const std::vector<LimitViolation> & viols, int bus_id)
{
    for (std::size_t k = 0; k < viols.size(); ++k) {
        if (viols[k].element_type == ViolationElementType::BUS && viols[k].element_id == bus_id) {
            return &viols[k];
        }
    }
    return nullptr;
}

// one row whose injection is the grid's own state, with the check on. A batch is neither
// copyable nor movable (and it works on a private copy of the grid taken at construction),
// so it is built by the caller and configured here.
void setup_one_row(TimeSeries & ts)
{
    ts.set_compute_physical_violations(true);
    ts.set_physical_violation_tol_mva(0.);
    RealMat load_p(1, 1);
    load_p << LOAD_P;
    ts.modify_load_p(load_p);
    RealMat load_q(1, 1);
    load_q << LOAD_Q;
    ts.modify_load_q(load_q);
}

}  // namespace

TEST_CASE("every violation type says what kind of statement it is", "[batch][physical]")
{
    // the three categories are the point of the taxonomy: an operational limit the grid may
    // leave, a physical one it cannot, and the solver's own verdict -- which is not a limit
    // at all (a divergence does not say whether the state was feasible).
    CHECK(violation_category(LimitViolationType::LOW_VOLTAGE) == ViolationCategory::OPERATIONAL);
    CHECK(violation_category(LimitViolationType::HIGH_VOLTAGE) == ViolationCategory::OPERATIONAL);
    CHECK(violation_category(LimitViolationType::CURRENT) == ViolationCategory::OPERATIONAL);
    CHECK(violation_category(LimitViolationType::LOW_Q) == ViolationCategory::PHYSICAL);
    CHECK(violation_category(LimitViolationType::HIGH_Q) == ViolationCategory::PHYSICAL);
    CHECK(violation_category(LimitViolationType::NOT_SIMULATED) == ViolationCategory::SOLVER);
    CHECK(violation_category(LimitViolationType::DIVERGENCE) == ViolationCategory::SOLVER);
    // and a violation carries its own, derived from its type
    const LimitViolation viol{ViolationElementType::BUS, 1, 0, LimitViolationType::HIGH_Q,
                              50., 10., std::string()};
    CHECK(viol.category() == ViolationCategory::PHYSICAL);
}

TEST_CASE("a bus reports the reactive power its machines had to produce", "[batch][physical]")
{
    SECTION("within what the machines own: nothing reported")
    {
        std::vector<GenSpec> gens{slack_gen(), GenSpec{GEN_BUS, V_SET, 10., -WIDE_Q, WIDE_Q, -1}};
        LSGrid grid = make_grid(gens);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        ts.compute(flat_start(grid), 30, 1e-11);
        REQUIRE(ts.converged_mask()[0] == 1);
        REQUIRE(ts.get_physical_violations().size() == 1);
        CHECK(ts.get_physical_violations()[0].empty());
        CHECK(ts.get_physical_violations_n().empty());
    }

    SECTION("beyond it: reported, with the value ac_pf publishes and the SUMMED limit")
    {
        // +/- 10 MVAr on the machine holding bus 1 at 1.05 pu against a 80 MW / 60 MVAr
        // load two lines away: it takes far more than 10 MVAr to do it.
        std::vector<GenSpec> gens{slack_gen(), GenSpec{GEN_BUS, V_SET, 10., -10., 10., -1}};
        const real_type q_bus = reference_bus_q(gens, GEN_BUS);
        REQUIRE(q_bus > 10.);  // the case is only meaningful if the capability IS exceeded

        LSGrid grid = make_grid(gens);
        std::vector<std::string> sub_names{"sub0", "sub1", "sub2", "sub3"};
        grid.set_substation_names(sub_names);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        ts.compute(flat_start(grid), 30, 1e-11);

        REQUIRE(ts.converged_mask()[0] == 1);
        const std::vector<LimitViolation> & viols = ts.get_physical_violations()[0];
        REQUIRE(viols.size() == 1);
        CHECK(viols[0].element_type == ViolationElementType::BUS);
        CHECK(viols[0].element_id == GEN_BUS);
        CHECK(viols[0].side == 0);
        CHECK(viols[0].violation_type == LimitViolationType::HIGH_Q);
        CHECK(viols[0].category() == ViolationCategory::PHYSICAL);
        CHECK(viols[0].value == Approx(q_bus).margin(1e-6));
        CHECK(viols[0].limit == Approx(10.));
        CHECK(viols[0].name == "sub1");  // the bus' substation, as for a voltage violation
        // the base ("n") case solves that same grid, so it reports the same thing
        REQUIRE(ts.get_physical_violations_n().size() == 1);
        CHECK(ts.get_physical_violations_n()[0].value == Approx(q_bus).margin(1e-6));
    }

    SECTION("the tolerance is what keeps a bus resting on its capability quiet")
    {
        std::vector<GenSpec> gens{slack_gen(), GenSpec{GEN_BUS, V_SET, 10., -10., 10., -1}};
        const real_type q_bus = reference_bus_q(gens, GEN_BUS);

        LSGrid grid = make_grid(gens);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        // a tolerance wider than the overshoot hides it
        ts.set_physical_violation_tol_mva(q_bus - 10. + 1.);
        ts.compute(flat_start(grid), 30, 1e-11);
        CHECK(ts.get_physical_violations()[0].empty());
    }
}

TEST_CASE("a bus is fine while one of its machines is not", "[batch][physical]")
{
    // THE test that separates a bus-level check from a per-machine one. Two machines hold
    // bus 1; together they own more than the bus needs, so the bus is feasible and nothing
    // is reported -- even though the reactive power lightsim2grid's sharing convention
    // hands one of them is beyond that machine's own limit. Which machine "produces" which
    // share of a bus' reactive power is a convention (LSGrid::_split_q_residual_per_bus),
    // not something the solver decides, so it must not decide what is reported.
    std::vector<GenSpec> wide{slack_gen(),
                              GenSpec{GEN_BUS, V_SET, 10., -WIDE_Q, WIDE_Q, -1},
                              GenSpec{GEN_BUS, V_SET, 10., -WIDE_Q, WIDE_Q, -1}};
    const real_type q_bus = reference_bus_q(wide, GEN_BUS);
    REQUIRE(q_bus > 20.);

    // asymmetric ranges: one machine can only produce (0 .. 0.75 q_bus), the other is
    // symmetric and small. Their SUM covers q_bus with 10% to spare.
    const real_type max_a = 0.75 * q_bus;
    const real_type max_b = 0.35 * q_bus;
    std::vector<GenSpec> gens{slack_gen(),
                              GenSpec{GEN_BUS, V_SET, 10., 0., max_a, -1},
                              GenSpec{GEN_BUS, V_SET, 10., -max_b, max_b, -1}};
    REQUIRE(max_a + max_b > q_bus);

    // ... and a per-machine check WOULD fire: at least one machine's own published
    // reactive output is beyond its own limit. Read off ac_pf, so this is a statement
    // about the grid and not about the code under test.
    LSGrid ref_grid = make_grid(gens);
    const RealVect q_gen = reference_gen_q(ref_grid);
    CHECK(q_gen(1) + q_gen(2) == Approx(q_bus).margin(1e-6));  // the split does not move the total
    const bool a_machine_is_over = (q_gen(1) > max_a) || (q_gen(2) > max_b);
    REQUIRE(a_machine_is_over);

    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    setup_one_row(ts);
    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.converged_mask()[0] == 1);
    CHECK(ts.get_physical_violations()[0].empty());

    SECTION("... and the same bus IS reported once the two of them cannot cover it")
    {
        std::vector<GenSpec> tight{slack_gen(),
                                   GenSpec{GEN_BUS, V_SET, 10., 0., 0.3 * q_bus, -1},
                                   GenSpec{GEN_BUS, V_SET, 10., -0.2 * q_bus, 0.2 * q_bus, -1}};
        LSGrid grid2 = make_grid(tight);
        grid2.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts2(grid2);
        setup_one_row(ts2);
        ts2.compute(flat_start(grid2), 30, 1e-11);
        REQUIRE(ts2.converged_mask()[0] == 1);
        const LimitViolation * viol = find_bus(ts2.get_physical_violations()[0], GEN_BUS);
        REQUIRE(viol != nullptr);
        CHECK(viol->violation_type == LimitViolationType::HIGH_Q);
        CHECK(viol->value == Approx(q_bus).margin(1e-6));
        CHECK(viol->limit == Approx(0.5 * q_bus));  // the SUM of the two machines' max_q
    }
}

TEST_CASE("a bus held through remote regulation reads its reactive power off the controllers",
          "[batch][physical][vctrl]")
{
    // gen 1 stands on bus 1 and regulates bus 3: its reactive output is a Jacobian unknown
    // of the VoltageControl extension, and bus 1 keeps a reactive equation of its own (so
    // its mismatch is ~0). Reading the mismatch alone would report ~0 instead of the real
    // reactive power, so this fails outright if the controller term is missing.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{GEN_BUS, V_SET, 10., -10., 10., NB_BUS - 1}};
    const real_type q_bus = reference_bus_q(gens, GEN_BUS);
    REQUIRE(std::abs(q_bus) > 10.);

    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    setup_one_row(ts);
    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.converged_mask()[0] == 1);

    const LimitViolation * viol = find_bus(ts.get_physical_violations()[0], GEN_BUS);
    REQUIRE(viol != nullptr);
    CHECK(viol->value == Approx(q_bus).margin(1e-6));
    CHECK(viol->violation_type == (q_bus > 0. ? LimitViolationType::HIGH_Q
                                              : LimitViolationType::LOW_Q));
}

TEST_CASE("an hvdc converter station's reactive capability counts towards its bus'",
          "[batch][physical][vctrl]")
{
    // a station regulating the bus it stands on holds that bus exactly like a local
    // generator does, and its [min_q, max_q] is in the same currency (MVAr). So it belongs
    // in the bus' capability: the SAME solution is feasible or not depending only on how
    // much the station brings.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{GEN_BUS, V_SET, 10., -10., 10., -1}};

    // how much bus 1 needs, with the station present (the station changes the solution --
    // it injects active power and pins the bus too -- so the reference must include it)
    LSGrid ref_grid = make_grid(gens, /*meshed=*/false, /*station_q_on_gen_bus=*/WIDE_Q);
    ref_grid.change_algorithm(AlgorithmType::NR_SparseLU);
    ref_grid.ac_pf(flat_start(ref_grid), 30, 1e-11);
    const real_type q_gen = RealVect(std::get<1>(ref_grid.get_gen_res()))(1);
    const real_type q_station = ref_grid.get_dclines()[0].res_q1_mvar;
    const real_type q_bus = q_gen + q_station;
    REQUIRE(q_bus > 10.);  // more than the generator alone owns

    SECTION("a station that brings enough makes the bus feasible")
    {
        // generator 10 MVAr + station (q_bus - 10 + 5) MVAr > q_bus
        LSGrid grid = make_grid(gens, /*meshed=*/false, q_bus - 10. + 5.);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        ts.compute(flat_start(grid), 30, 1e-11);
        REQUIRE(ts.converged_mask()[0] == 1);
        CHECK(find_bus(ts.get_physical_violations()[0], GEN_BUS) == nullptr);
    }

    SECTION("a station that does not is reported, against the summed capability")
    {
        const real_type station_q = 0.25 * q_bus;
        LSGrid grid = make_grid(gens, /*meshed=*/false, station_q);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        ts.compute(flat_start(grid), 30, 1e-11);
        REQUIRE(ts.converged_mask()[0] == 1);
        const LimitViolation * viol = find_bus(ts.get_physical_violations()[0], GEN_BUS);
        REQUIRE(viol != nullptr);
        CHECK(viol->violation_type == LimitViolationType::HIGH_Q);
        // the value is the bus' whole reactive power (generator + station), and the limit
        // the sum of the two capabilities -- 10 MVAr of generator plus the station's
        CHECK(viol->value == Approx(q_bus).margin(1e-6));
        CHECK(viol->limit == Approx(10. + station_q));
    }
}

TEST_CASE("a voltage-mode SVC's susceptance range counts towards its bus'",
          "[batch][physical][vctrl]")
{
    // an SVC's capability is a SUSCEPTANCE range, so what it is worth in MVAr depends on
    // the solved voltage: q = b . |V|^2 . sn_mva. It holds bus 3 here (alone -- LSGrid
    // refuses an SVC sharing its regulated bus with any other controller), so the bus'
    // reactive power is the SVC's own and the check reduces to "is b enough at this V".
    const real_type b_pu = 2.0;  // +/- 2 pu of susceptance: plenty at this voltage
    LSGrid ref_grid = make_svc_grid(b_pu);
    ref_grid.change_algorithm(AlgorithmType::NR_SparseLU);
    ref_grid.ac_pf(flat_start(ref_grid), 30, 1e-11);
    const real_type q_svc = ref_grid.get_svcs()[0].res_q_mvar;
    const real_type v_svc = ref_grid.get_svcs()[0].res_v_kv / 138.;  // pu
    const real_type q_max = b_pu * v_svc * v_svc * 100.;             // b . |V|^2 . sn_mva
    REQUIRE(q_svc > 0.);

    SECTION("enough susceptance: nothing reported")
    {
        REQUIRE(q_svc < q_max);  // the case is only meaningful if b really does cover it
        LSGrid grid = make_svc_grid(b_pu);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        ts.compute(flat_start(grid), 30, 1e-11);
        REQUIRE(ts.converged_mask()[0] == 1);
        CHECK(find_bus(ts.get_physical_violations()[0], SVC_BUS) == nullptr);
    }

    SECTION("not enough: reported, against b . |V|^2 . sn_mva and not against b")
    {
        // half the susceptance the solution needs
        const real_type b_small = 0.5 * q_svc / (v_svc * v_svc * 100.);
        LSGrid grid = make_svc_grid(b_small);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        ts.compute(flat_start(grid), 30, 1e-11);
        REQUIRE(ts.converged_mask()[0] == 1);
        const LimitViolation * viol = find_bus(ts.get_physical_violations()[0], SVC_BUS);
        REQUIRE(viol != nullptr);
        CHECK(viol->violation_type == LimitViolationType::HIGH_Q);
        CHECK(viol->value == Approx(q_svc).margin(1e-6));
        // the reported limit is in MVAr, evaluated at the solved voltage
        CHECK(viol->limit == Approx(b_small * v_svc * v_svc * 100.).margin(1e-6));
        CHECK(viol->limit == Approx(0.5 * q_svc).margin(1e-6));
    }

    SECTION("an asymmetric range pins the sign: a positive susceptance PRODUCES reactive power")
    {
        // this bus needs reactive power produced (q_svc > 0), so a purely capacitive SVC
        // ([0, b]) covers it and a purely inductive one ([-b, 0]) cannot. Were the sign
        // convention inverted, the two verdicts would swap -- and a symmetric range, which
        // is what the sections above use, could never tell.
        LSGrid capacitive = make_svc_grid_asym(0., b_pu);
        capacitive.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts_cap(capacitive);
        setup_one_row(ts_cap);
        ts_cap.compute(flat_start(capacitive), 30, 1e-11);
        REQUIRE(ts_cap.converged_mask()[0] == 1);
        CHECK(find_bus(ts_cap.get_physical_violations()[0], SVC_BUS) == nullptr);

        LSGrid inductive = make_svc_grid_asym(-b_pu, 0.);
        inductive.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts_ind(inductive);
        setup_one_row(ts_ind);
        ts_ind.compute(flat_start(inductive), 30, 1e-11);
        REQUIRE(ts_ind.converged_mask()[0] == 1);
        const LimitViolation * viol = find_bus(ts_ind.get_physical_violations()[0], SVC_BUS);
        REQUIRE(viol != nullptr);
        CHECK(viol->violation_type == LimitViolationType::HIGH_Q);
        CHECK(viol->limit == Approx(0.).margin(1e-9));  // b_max = 0: it can produce nothing
        CHECK(viol->value == Approx(q_svc).margin(1e-6));
    }
}

TEST_CASE("the physical-limit checks are opt in, and say so when they are off",
          "[batch][physical]")
{
    std::vector<GenSpec> gens{slack_gen(), GenSpec{GEN_BUS, V_SET, 10., -10., 10., -1}};
    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    CHECK_FALSE(ts.get_compute_physical_violations());
    CHECK_THROWS_AS(ts.get_physical_violations(), std::runtime_error);
    CHECK_THROWS_AS(ts.get_physical_violations_n(), std::runtime_error);
    CHECK(ts.get_physical_violation_tol_mva() == Approx(1e-4));
    CHECK_THROWS_AS(ts.set_physical_violation_tol_mva(-1.), std::runtime_error);
}

TEST_CASE("the reactive half needs an algorithm that can feed it; the hvdc half does not",
          "[batch][physical]")
{
    // Every built-in AC family publishes its per-bus mismatch (NR, fast-decoupled AND
    // Gauss-Seidel -- see BaseAlgo::FILLS_BUS_MISMATCH and its overrides), so the only
    // algorithm this could refuse is a plugin that has not opted in. DC is NOT refused: a DC
    // powerflow has no reactive power at all, so the reactive half is not applicable rather
    // than missing, and the hvdc half -- which needs nothing but the bus angles -- still
    // runs.
    //
    // The algorithm a batch runs is its OWN: inherited from the grid at construction, and
    // changed afterwards through the BATCH, never through the grid (see
    // BaseBatchSolverSynch's constructor).
    LSGrid grid = make_droop_grid(/*pmax_1to2=*/5., /*pmax_2to1=*/500.);
    TimeSeries ts(grid);
    ts.change_algorithm(AlgorithmType::DC_SparseLU);
    ts.set_compute_physical_violations(true);
    ts.set_physical_violation_tol_mva(0.);
    RealMat load_p(1, 1);
    load_p << LOAD_P;
    ts.modify_load_p(load_p);
    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.converged_mask()[0] == 1);

    // the hvdc limit is reported, and against the DC solution's own angles: the reference is
    // a DC single shot of the same grid
    LSGrid ref = make_droop_grid(/*pmax_1to2=*/5., /*pmax_2to1=*/500.);
    ref.change_algorithm(AlgorithmType::DC_SparseLU);
    ref.dc_pf(flat_start(ref), 30, 1e-11);
    const real_type p_dc = -ref.get_dclines()[0].res_p1_mw;  // leaving bus 1 into the hvdc
    REQUIRE(p_dc > 5.);

    const std::vector<LimitViolation> & viols = ts.get_physical_violations()[0];
    REQUIRE(viols.size() == 1);
    CHECK(viols[0].element_type == ViolationElementType::HVDC);
    CHECK(viols[0].violation_type == LimitViolationType::HIGH_P);
    CHECK(viols[0].value == Approx(p_dc).margin(1e-6));
    // ... and nothing about reactive power, which a DC solve does not have
    for (std::size_t k = 0; k < viols.size(); ++k) {
        CHECK(viols[k].violation_type != LimitViolationType::LOW_Q);
        CHECK(viols[k].violation_type != LimitViolationType::HIGH_Q);
    }
}

TEST_CASE("an angle-droop hvdc beyond what its converters can transmit is reported",
          "[batch][physical][hvdc]")
{
    // `status_droop` is an INPUT of the solve (0 = linear), so nothing saturates the droop:
    // the flow is whatever the angle difference asks for, and a row can converge with a
    // converter transmitting power it does not have. That is the condition OpenLoadFlow's
    // HvdcAcEmulationLimits outer loop acts on.
    LSGrid ref = make_droop_grid(/*pmax_1to2=*/1000., /*pmax_2to1=*/1000.);
    ref.change_algorithm(AlgorithmType::NR_SparseLU);
    ref.ac_pf(flat_start(ref), 30, 1e-11);
    // no losses on this fixture, so the flow leaving bus 1 into the hvdc is p0 + k.dtheta
    const real_type p_flow = -ref.get_dclines()[0].res_p1_mw;
    CHECK(p_flow == Approx(ref.get_dclines()[0].res_p2_mw).margin(1e-9));
    REQUIRE(p_flow > 1.);  // it really does flow 1 -> 2

    SECTION("within the limits: nothing reported")
    {
        LSGrid grid = make_droop_grid(p_flow + 5., 1000.);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        ts.compute(flat_start(grid), 30, 1e-11);
        REQUIRE(ts.converged_mask()[0] == 1);
        CHECK(ts.get_physical_violations()[0].empty());
    }

    SECTION("beyond them: reported, with the flow ac_pf publishes and that direction's max")
    {
        const real_type pmax = p_flow - 5.;
        LSGrid grid = make_droop_grid(pmax, 1000.);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        std::vector<std::string> names{"dc_link"};
        grid.set_dcline_names(names);
        TimeSeries ts(grid);
        setup_one_row(ts);
        ts.compute(flat_start(grid), 30, 1e-11);
        REQUIRE(ts.converged_mask()[0] == 1);

        const std::vector<LimitViolation> & viols = ts.get_physical_violations()[0];
        REQUIRE(viols.size() == 1);
        CHECK(viols[0].element_type == ViolationElementType::HVDC);
        CHECK(viols[0].element_id == 0);
        CHECK(viols[0].side == 1);  // the flow leaves side 1
        CHECK(viols[0].violation_type == LimitViolationType::HIGH_P);
        CHECK(viols[0].category() == ViolationCategory::PHYSICAL);
        CHECK(viols[0].value == Approx(p_flow).margin(1e-6));
        CHECK(viols[0].limit == Approx(pmax));
        CHECK(viols[0].name == "dc_link");
        // the base ("n") case solves the same grid, so it reports the same thing
        REQUIRE(ts.get_physical_violations_n().size() == 1);
        CHECK(ts.get_physical_violations_n()[0].value == Approx(p_flow).margin(1e-6));
    }

    SECTION("each direction is judged against ITS OWN maximum")
    {
        // the flow goes 1 -> 2 here, so only pmax_1to2 can bound it: a tiny pmax_2to1 is
        // irrelevant, and swapping the two changes the verdict. A check comparing |p|
        // against whichever limit came to hand would pass one of these and fail the other.
        LSGrid wrong_way = make_droop_grid(/*pmax_1to2=*/1000., /*pmax_2to1=*/1.);
        wrong_way.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts_ok(wrong_way);
        setup_one_row(ts_ok);
        ts_ok.compute(flat_start(wrong_way), 30, 1e-11);
        REQUIRE(ts_ok.converged_mask()[0] == 1);
        CHECK(ts_ok.get_physical_violations()[0].empty());

        LSGrid right_way = make_droop_grid(/*pmax_1to2=*/1., /*pmax_2to1=*/1000.);
        right_way.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts_bad(right_way);
        setup_one_row(ts_bad);
        ts_bad.compute(flat_start(right_way), 30, 1e-11);
        REQUIRE(ts_bad.converged_mask()[0] == 1);
        REQUIRE(ts_bad.get_physical_violations()[0].size() == 1);
        CHECK(ts_bad.get_physical_violations()[0][0].side == 1);
        CHECK(ts_bad.get_physical_violations()[0][0].limit == Approx(1.));
    }
}

TEST_CASE("a droop flowing 2 -> 1 is reported on side 2, against pmax_2to1",
          "[batch][physical][hvdc]")
{
    // the same physical flow (bus 1 -> bus 3, pulled by the load), described from the other
    // end: side 1 now sits at the load bus, so the hvdc carries power OUT of its side 2 and
    // it is pmax_2to1 that bounds it.
    LSGrid ref = make_droop_grid(1000., 1000., 30., 400., 0, /*reverse_ends=*/true);
    ref.change_algorithm(AlgorithmType::NR_SparseLU);
    ref.ac_pf(flat_start(ref), 30, 1e-11);
    const real_type p_flow = ref.get_dclines()[0].res_p1_mw;  // positive: received at side 1
    REQUIRE(p_flow > 1.);

    const real_type pmax = p_flow - 5.;
    LSGrid grid = make_droop_grid(1000., pmax, 30., 400., 0, /*reverse_ends=*/true);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    setup_one_row(ts);
    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.converged_mask()[0] == 1);

    const std::vector<LimitViolation> & viols = ts.get_physical_violations()[0];
    REQUIRE(viols.size() == 1);
    CHECK(viols[0].side == 2);  // the flow leaves side 2
    CHECK(viols[0].violation_type == LimitViolationType::HIGH_P);
    CHECK(viols[0].value == Approx(p_flow).margin(1e-6));  // reported positive
    CHECK(viols[0].limit == Approx(pmax));

    // ... and the other direction's maximum has no say in it
    LSGrid other_way = make_droop_grid(/*pmax_1to2=*/1., 1000., 30., 400., 0, true);
    other_way.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts_ok(other_way);
    setup_one_row(ts_ok);
    ts_ok.compute(flat_start(other_way), 30, 1e-11);
    REQUIRE(ts_ok.converged_mask()[0] == 1);
    CHECK(ts_ok.get_physical_violations()[0].empty());
}

TEST_CASE("a droop the caller already saturated is not reported", "[batch][physical][hvdc]")
{
    // status_droop != 0 means someone has run the saturation logic between two solves (what
    // LSGrid::set_status_droop_hvdc is for): the solver then PINS the flow at the very limit
    // this would compare against, so there is nothing left to detect.
    LSGrid grid = make_droop_grid(/*pmax_1to2=*/5., /*pmax_2to1=*/1000.);
    grid.set_status_droop_hvdc(0, 1);  // saturated 1 -> 2
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    setup_one_row(ts);
    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.converged_mask()[0] == 1);
    CHECK(ts.get_physical_violations()[0].empty());
    CHECK(ts.get_physical_violations_n().empty());
}

TEST_CASE("the hvdc flow is read in the batch's own bus labelling", "[batch][physical][hvdc]")
{
    // Deactivated spare buses shift the solver ids away from the grid ones, so a check
    // reading theta through the wrong map takes another bus' angle -- a wrong number, not an
    // error. The reported flow must still be the one a single shot publishes.
    const int nb_spare = 3;
    LSGrid ref = make_droop_grid(1000., 1000., 30., 400., nb_spare);
    ref.change_algorithm(AlgorithmType::NR_SparseLU);
    ref.ac_pf(flat_start(ref), 30, 1e-11);
    const real_type p_flow = -ref.get_dclines()[0].res_p1_mw;
    REQUIRE(p_flow > 1.);

    LSGrid grid = make_droop_grid(p_flow - 5., 1000., 30., 400., nb_spare);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    setup_one_row(ts);
    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.converged_mask()[0] == 1);
    REQUIRE(ts.get_physical_violations()[0].size() == 1);
    CHECK(ts.get_physical_violations()[0][0].value == Approx(p_flow).margin(1e-6));
}

TEST_CASE("a second compute() that reuses the base case reports the same thing",
          "[batch][physical]")
{
    // `reuse_base_case` (on by default) skips the "n" solve of the second call, which
    // leaves the member algorithm holding the mismatch of the LAST ROW of the first call.
    // Deriving the base case's report from that would report a row as the base case; both
    // calls must report the very same thing.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{GEN_BUS, V_SET, 10., -10., 10., -1}};
    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    ts.set_compute_physical_violations(true);
    ts.set_physical_violation_tol_mva(0.);
    // two rows with DIFFERENT loads, so the last row's reactive power is not the base
    // case's and a stale read is visible
    RealMat load_p(2, 1);
    load_p << LOAD_P, 0.5 * LOAD_P;
    ts.modify_load_p(load_p);
    RealMat load_q(2, 1);
    load_q << LOAD_Q, 0.5 * LOAD_Q;
    ts.modify_load_q(load_q);

    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.get_compute_physical_violations());
    REQUIRE(ts.get_physical_violations_n().size() == 1);
    const real_type q_n_first = ts.get_physical_violations_n()[0].value;
    REQUIRE(ts.converged_mask()[1] == 1);
    REQUIRE(ts.get_physical_violations().size() == 2);
    REQUIRE(ts.get_physical_violations()[1].size() == 1);
    const real_type q_row1_first = ts.get_physical_violations()[1][0].value;
    REQUIRE(std::abs(q_row1_first - q_n_first) > 1.);  // the rows differ from the base case

    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.base_case_was_reused());
    REQUIRE(ts.get_physical_violations_n().size() == 1);
    CHECK(ts.get_physical_violations_n()[0].value == Approx(q_n_first));
    REQUIRE(ts.get_physical_violations()[1].size() == 1);
    CHECK(ts.get_physical_violations()[1][0].value == Approx(q_row1_first));
}

TEST_CASE("Gauss-Seidel reports the same reactive power as Newton-Raphson", "[batch][physical]")
{
    // the value is read off the ALGORITHM's mismatch, so it has to be right for every
    // family that publishes one -- not just for the Newton-Raphson the rest of this file
    // runs on.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{GEN_BUS, V_SET, 10., -10., 10., -1}};
    const real_type q_bus = reference_bus_q(gens, GEN_BUS);

    LSGrid grid = make_grid(gens);
    TimeSeries ts(grid);
    ts.change_algorithm(AlgorithmType::GaussSeidel);
    setup_one_row(ts);
    ts.compute(flat_start(grid), 10000, 1e-9);
    REQUIRE(ts.converged_mask()[0] == 1);
    const LimitViolation * viol = find_bus(ts.get_physical_violations()[0], GEN_BUS);
    REQUIRE(viol != nullptr);
    CHECK(viol->value == Approx(q_bus).margin(1e-4));
}

TEST_CASE("a contingency row reports its own reactive power, not the base case's",
          "[batch][physical][contingency]")
{
    // meshed: line 2 (bus2--bus3) can go without islanding the load, which now reaches bus
    // 3 through line 3 (bus0--bus3). The outage moves the reactive flows, so the base case
    // and the contingency row must report DIFFERENT values -- each matching the ac_pf of
    // the corresponding grid.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{GEN_BUS, V_SET, 10., -10., 10., -1}};
    const real_type q_n = reference_bus_q(gens, GEN_BUS, /*meshed=*/true);
    LSGrid ref_c = make_grid(gens, /*meshed=*/true);
    ref_c.deactivate_powerline(2);
    const real_type q_c = reference_gen_q(ref_c)(1);
    REQUIRE(q_n > 10.);
    REQUIRE(q_c > 10.);
    REQUIRE(std::abs(q_c - q_n) > 1e-3);  // the outage really does change it

    LSGrid grid = make_grid(gens, /*meshed=*/true);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    ContingencyAnalysis ca(grid);
    ca.set_compute_physical_violations(true);
    ca.set_physical_violation_tol_mva(0.);
    ca.add_n1(2);
    ca.compute(flat_start(grid), 30, 1e-11);

    REQUIRE(ca.converged_mask()[0] == 1);
    const LimitViolation * viol_n = find_bus(ca.get_physical_violations_n(), GEN_BUS);
    REQUIRE(viol_n != nullptr);
    CHECK(viol_n->value == Approx(q_n).margin(1e-6));
    CHECK(viol_n->limit == Approx(10.));

    REQUIRE(ca.get_physical_violations().size() == 1);
    const LimitViolation * viol_c = find_bus(ca.get_physical_violations()[0], GEN_BUS);
    REQUIRE(viol_c != nullptr);
    CHECK(viol_c->value == Approx(q_c).margin(1e-6));
}

TEST_CASE("a row that was never simulated reports nothing at all", "[batch][physical][contingency]")
{
    // line 0 (bus0--bus1) carries the whole radial feeder: taking it out islands everything
    // past bus 0, so the row is skipped before the solver. A skipped row must come back
    // empty rather than with a stale or zero-voltage answer.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{GEN_BUS, V_SET, 10., -10., 10., -1}};
    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    ContingencyAnalysis ca(grid);
    ca.set_compute_physical_violations(true);
    ca.add_n1(0);
    ca.compute(flat_start(grid), 30, 1e-11);

    CHECK(ca.converged_mask()[0] == 0);
    REQUIRE(ca.get_physical_violations().size() == 1);
    CHECK(ca.get_physical_violations()[0].empty());
    // ... while the base case, which did converge, still reports on its own
    CHECK(find_bus(ca.get_physical_violations_n(), GEN_BUS) != nullptr);
}

TEST_CASE("a row that disconnects a machine checks the bus against what is left",
          "[batch][physical][scenario_sweep]")
{
    // two machines hold bus 1 and together cover what it needs; row 1 disconnects one of
    // them, and the bus' capability shrinks with it. So the same bus is feasible in row 0
    // and not in row 1 -- and the reactive power reported in row 1 is what a grid with that
    // machine deactivated publishes.
    std::vector<GenSpec> wide{slack_gen(),
                              GenSpec{GEN_BUS, V_SET, 10., -WIDE_Q, WIDE_Q, -1},
                              GenSpec{GEN_BUS, V_SET, 10., -WIDE_Q, WIDE_Q, -1}};
    const real_type q_both = reference_bus_q(wide, GEN_BUS);
    // each machine owns 60% of what the bus needs: the two of them cover it, one does not
    const real_type max_each = 0.6 * q_both;
    std::vector<GenSpec> gens{slack_gen(),
                              GenSpec{GEN_BUS, V_SET, 10., -max_each, max_each, -1},
                              GenSpec{GEN_BUS, V_SET, 10., -max_each, max_each, -1}};
    LSGrid ref_alone = make_grid(gens);
    ref_alone.deactivate_gen(2);
    const real_type q_alone = reference_gen_q(ref_alone)(1);
    REQUIRE(q_alone > max_each);  // bus 1 alone cannot cover it

    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    ScenarioSweep sweep(grid);
    sweep.set_compute_physical_violations(true);
    sweep.set_physical_violation_tol_mva(0.);
    RealMat load_p(2, 1);
    load_p << LOAD_P, LOAD_P;
    sweep.modify_load_p(load_p);
    RealMat load_q(2, 1);
    load_q << LOAD_Q, LOAD_Q;
    sweep.modify_load_q(load_q);
    BoolMat gen_off(2, 3);
    gen_off << false, false, false,
               false, false, true;  // row 1 disconnects gen 2
    sweep.set_contingency_gens(gen_off);

    sweep.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(sweep.converged_mask()[0] == 1);
    REQUIRE(sweep.converged_mask()[1] == 1);

    // row 0: both machines on, together they cover the bus
    CHECK(sweep.get_physical_violations()[0].empty());
    // row 1: one machine left, and it cannot
    const LimitViolation * row1 = find_bus(sweep.get_physical_violations()[1], GEN_BUS);
    REQUIRE(row1 != nullptr);
    CHECK(row1->violation_type == LimitViolationType::HIGH_Q);
    CHECK(row1->value == Approx(q_alone).margin(1e-6));
    CHECK(row1->limit == Approx(max_each));  // only the machine that is still on
}
