// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LIMITVIOLATION_H
#define LIMITVIOLATION_H

#include "Utils.hpp"
#include <string>

namespace ls2g {

// the kind of element on which a limit was violated
enum class LS2G_API ViolationElementType : int {
    BUS = 0,
    LINE = 1,
    TRAFO = 2,
    GRID = 3,  // the whole grid / contingency, not a specific element (see LimitViolationType::DIVERGENCE)
    // an hvdc line, by its own id (see LimitViolationType::HIGH_P); also the release of one of
    // its VSC stations frozen at a reactive limit (LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q,
    // `side` the station's end, see LSGrid::set_hvdc_can_be_pv)
    HVDC = 4,
    // a generator, by its own id. Its ACTIVE power (LOW_P / HIGH_P) and, for a PQ machine
    // pinned at a reactive limit, the voltage of the bus it would regulate
    // (LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q), and, for a machine regulating a REMOTE
    // bus, the voltage of its own bus (LOW_VOLTAGE_REMOTE_CONTROL / HIGH_VOLTAGE_REMOTE_CONTROL,
    // see RemoteVoltageControlCheck.hpp). A REGULATING machine's reactive
    // violation is reported on the BUS instead, because how a bus' reactive power is
    // divided between its machines is a convention, while the active one is divided by the
    // participation factors the caller chose. See BusQCheck.hpp / GenPCheck.hpp /
    // GenPvReleaseCheck.hpp.
    GENERATOR = 5,
    // a storage unit, by its own id -- same statement as GENERATOR, since a storage unit
    // takes a share of the distributed slack under the same rule (SlackParticipation).
    // /!\ `value` and `limit` are in the GENERATOR convention (positive = injected into
    // the grid), like the unit's min_q_mvar / max_q_mvar and unlike its target_p_mw /
    // res_p_mw, which this library stores in the load convention.
    STORAGE = 6,
    // a static var compensator, by its own id: an idle SVC carrying a standby automaton
    // that the voltage of the bus it regulates would switch on
    // (LOW_VOLTAGE_SVC_STANDBY / HIGH_VOLTAGE_SVC_STANDBY, see SvcStandbyCheck.hpp), or a
    // fixed-Q SVC an outer loop froze at a reactive limit that would regulate again
    // (LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q, see GenPvReleaseCheck.hpp)
    SVC = 7
};

// the kind of limit that was violated
enum class LS2G_API LimitViolationType : int {
    LOW_VOLTAGE = 0,
    HIGH_VOLTAGE = 1,
    CURRENT = 2,
    NOT_SIMULATED = 3,  // a pre-check (graph connectivity) skipped the contingency: the solver was never invoked
    DIVERGENCE = 4,  // the solver was invoked (for the contingency, or the pre-contingency "n" case) but did not converge
    // The reactive power the machines holding ONE BUS' voltage had to produce left what
    // they are capable of -- the summed capability of its voltage-regulating generators,
    // hvdc converter stations and voltage-mode SVCs. Unlike the three above this is not a
    // limit one may choose to exceed: a machine simply cannot produce reactive power it
    // does not have, so the converged solution the solver returned is not a state the grid
    // can reach -- its voltage set-point could not actually be held. See
    // ViolationCategory::PHYSICAL and compute_physical_violations. Checked per bus, not per
    // machine, because which machine of a bus "produces" which share of its reactive power
    // is a modelling convention (see LSGrid::_split_q_residual_per_bus), while what the
    // bus as a whole can produce is not.
    LOW_Q = 5,
    HIGH_Q = 6,
    // An active power above what the element can deliver. On an HVDC: the active power of
    // an angle-droop ("AC emulation") line left what its converters can transmit, with
    // `side` saying which direction and therefore which limit (1 for 1 -> 2, against
    // pmax_1to2_mw; 2 for the other way, against pmax_2to1_mw). On a GENERATOR: the
    // distributed slack -- solved inside the Jacobian, by participation factors that know
    // nothing about limits -- asked a machine for more than its max_p_mw.
    // Physical for the same reason as LOW_Q / HIGH_Q and unlike CURRENT: a branch above
    // its thermal rating is a state the grid reaches and should not sit in, while a
    // converter beyond its maximum power is a state it does not reach -- its own control
    // saturates first, which is what `status_droop = +/-1` models. Reported (see
    // compute_physical_violations), never enforced: nothing clamps the droop.
    // On a STORAGE: the same as on a generator, in the generator convention (see above).
    HIGH_P = 7,
    // ... and the other way: below what it can deliver. Only a generator or a storage unit
    // carrying the distributed slack can be reported here (an hvdc line's two directions
    // are two HIGH_P with a different `side`, see above).
    LOW_P = 8,
    // A PQ GENERATOR an outer loop pinned at its MINIMUM reactive power (flagged as such by
    // the caller, see LSGrid::set_gen_can_be_pv) whose regulated bus sits BELOW the target
    // it would hold: it absorbs too much for that target, and OpenLoadFlow's ReactiveLimits
    // loop would switch it back to PV (its PQ -> PV direction, the mirror of LOW_Q /
    // HIGH_Q). Physical for the same reason: the converged solution assumes a control the
    // loop would not leave in place. `value` the regulated bus' voltage and `limit` the
    // target, both in kV; the "LOW" in the name says which side the excess is on (value
    // below limit), and downstream code relies on it. See GenPvReleaseCheck.hpp. Also on an
    // SVC flagged as frozen at a reactive limit (LSGrid::set_svc_can_be_pv).
    LOW_VOLTAGE_AT_MIN_Q = 9,
    // ... and the mirror: pinned at its MAXIMUM, regulated bus ABOVE the target (value
    // above limit).
    HIGH_VOLTAGE_AT_MAX_Q = 10,
    // A non-regulating SVC carrying a standby automaton (flagged by the caller, see
    // LSGrid::set_svc_standby) whose regulated bus sits BELOW the automaton's low voltage
    // threshold: OpenLoadFlow's MonitoringVoltageOuterLoop would switch it to voltage
    // control. Physical for the same reason as LOW_VOLTAGE_AT_MIN_Q: the converged
    // solution assumes a control the loop would not leave in place. `value` the regulated
    // bus' voltage and `limit` the threshold, both in kV (value below limit). See
    // SvcStandbyCheck.hpp.
    LOW_VOLTAGE_SVC_STANDBY = 11,
    // ... and the mirror: regulated bus ABOVE the high voltage threshold (value above limit).
    HIGH_VOLTAGE_SVC_STANDBY = 12,
    // A GENERATOR regulating a REMOTE bus whose OWN bus sits BELOW the lowest voltage the
    // caller deems realistic for a controller (LSGrid::set_remote_voltage_control_vm_range):
    // OpenLoadFlow's ReactiveLimits loop, in its "robust" remote voltage control mode,
    // switches such a controller to PQ at its target reactive power rather than let it hold
    // the remote target at that price. Physical for the same reason as LOW_VOLTAGE_AT_MIN_Q:
    // the converged solution assumes a control the loop would not leave in place. `value`
    // the generator's own bus voltage and `limit` the threshold, both in kV of that bus
    // (value below limit). See RemoteVoltageControlCheck.hpp.
    LOW_VOLTAGE_REMOTE_CONTROL = 13,
    // ... and the mirror: own bus ABOVE the highest realistic voltage (value above limit).
    HIGH_VOLTAGE_REMOTE_CONTROL = 14,
    // An angle-droop ("AC emulation") HVDC line an outer loop froze at its active power
    // limit (flagged by the caller, see LSGrid::set_hvdc_ac_emulation_frozen) whose droop
    // would now ask for LESS than that limit: OpenLoadFlow's AcHvdcAcEmulationLimits loop
    // would not saturate it, the line would follow its droop. Physical for the same reason
    // as LOW_VOLTAGE_AT_MIN_Q: the converged solution assumes a control the loop would not
    // leave in place. `side` the direction it is frozen in (1: 1 -> 2), `value` the flow its
    // droop asks for in that direction and `limit` the limit, both in MW (value below limit).
    // See HvdcPCheck.hpp.
    HVDC_AC_EMULATION_RELEASE = 15,
    // The active power the single slack bus absorbed is above OpenLoadFlow's
    // slackBusPMaxMismatch: its DistributedSlack loop would share it on the units that take
    // part in the slack (see outer_loop/DistributedSlackLoop.hpp). element_type GRID,
    // `value` the mismatch (MW, > 0: the units must inject more) and `limit` the threshold.
    SLACK_MISMATCH = 16
};

/**
 * What KIND of statement a violation is -- a property of its `LimitViolationType`, and the
 * first thing to read: the three kinds do not mean the same thing and must not be acted on
 * the same way.
 *
 * A pure function of the type (`violation_category`), not a field, so the two can never
 * disagree.
 */
enum class LS2G_API ViolationCategory : int {
    /// A limit chosen by an operator, which the grid CAN leave: a bus outside its voltage
    /// range, a branch above its current rating. The solution is a state the grid can
    /// reach -- it is just a state nobody wants to sit in, and things eventually break.
    /// LOW_VOLTAGE, HIGH_VOLTAGE, CURRENT.
    OPERATIONAL = 0,
    /// A limit of the equipment itself, which nothing can leave. A violation here says the
    /// converged solution is NOT physically realizable, whatever anyone decides: the
    /// control it assumes (a voltage set-point held by machines that would have to produce
    /// reactive power they do not have, an hvdc converter transmitting more than it can, a
    /// distributed slack asking a machine for power it does not have) cannot happen. It is a
    /// statement about the model's assumptions, not about how the grid is being operated.
    /// LOW_Q, HIGH_Q, LOW_P, HIGH_P, LOW_VOLTAGE_AT_MIN_Q, HIGH_VOLTAGE_AT_MAX_Q,
    /// LOW_VOLTAGE_SVC_STANDBY, HIGH_VOLTAGE_SVC_STANDBY, LOW_VOLTAGE_REMOTE_CONTROL,
    /// HIGH_VOLTAGE_REMOTE_CONTROL, HVDC_AC_EMULATION_RELEASE, SLACK_MISMATCH.
    PHYSICAL = 1,
    /// Not a limit at all: what the solver did. A divergence in particular says nothing
    /// about the grid -- the state may be perfectly feasible and the algorithm simply
    /// failed to find it, or there may be no solution; this does not distinguish the two.
    /// NOT_SIMULATED, DIVERGENCE.
    SOLVER = 2
};

/// The category of a violation type. Every type has exactly one; see ViolationCategory.
inline ViolationCategory violation_category(LimitViolationType violation_type) noexcept
{
    switch(violation_type){
        case LimitViolationType::LOW_VOLTAGE:
        case LimitViolationType::HIGH_VOLTAGE:
        case LimitViolationType::CURRENT:
            return ViolationCategory::OPERATIONAL;
        case LimitViolationType::LOW_Q:
        case LimitViolationType::HIGH_Q:
        case LimitViolationType::HIGH_P:
        case LimitViolationType::LOW_P:
        case LimitViolationType::LOW_VOLTAGE_AT_MIN_Q:
        case LimitViolationType::HIGH_VOLTAGE_AT_MAX_Q:
        case LimitViolationType::LOW_VOLTAGE_SVC_STANDBY:
        case LimitViolationType::HIGH_VOLTAGE_SVC_STANDBY:
        case LimitViolationType::LOW_VOLTAGE_REMOTE_CONTROL:
        case LimitViolationType::HIGH_VOLTAGE_REMOTE_CONTROL:
        case LimitViolationType::HVDC_AC_EMULATION_RELEASE:
        case LimitViolationType::SLACK_MISMATCH:
            return ViolationCategory::PHYSICAL;
        default:  // NOT_SIMULATED, DIVERGENCE
            return ViolationCategory::SOLVER;
    }
}

// Default relative tolerance of the OPERATIONAL checks (`violation_rel_tol`): a value
// is reported only when it is beyond its (threshold-tightened) limit by more than this
// fraction of that limit -- v < low_eff * (1 - rel_tol), v > high_eff * (1 + rel_tol),
// amps > threshold * limit * (1 + rel_tol). A bus a regulator holds exactly at its vmax
// comes out of a solve at vmax +/- a few ulps, and the last bit then decided whether it
// was reported (it differed between the one-off solve, the batch and gpusim2grid).
// 1e-9 is ~0.4 mV on a 400 kV bus: far below anything physical, far above rounding.
constexpr real_type DEFAULT_VIOLATION_REL_TOL = 1e-9;

// a single limit violation, as detected by ContingencyAnalysis (see compute_limit_violations)
// and by the physical-limit checks (see compute_physical_violations)
struct LS2G_API LimitViolation {
    ViolationElementType element_type;
    // grid-model bus id for BUS ; local (0-based, own type) line / trafo / hvdc line /
    // generator / storage / svc id otherwise ; unused (-1) for GRID
    int element_id;
    // 1 or 2 for LINE / TRAFO (the terminal) and for HVDC (the direction: 1 means the flow
    // leaves side 1, ie 1 -> 2) ; unused (0) for BUS / GENERATOR / STORAGE / SVC / GRID
    int side;
    LimitViolationType violation_type;
    // value reached: MVAr for LOW_Q / HIGH_Q (what the machines holding that bus had to
    // produce), MW for LOW_P / HIGH_P (an hvdc line's flow leaving `side`, always positive;
    // a generator's or a storage unit's own converged active power, signed and in the
    // generator convention), kV for LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q (the voltage
    // of the bus the pinned generator would regulate) and for LOW_VOLTAGE_SVC_STANDBY /
    // HIGH_VOLTAGE_SVC_STANDBY (the voltage of the bus the standby SVC regulates) and for
    // LOW_VOLTAGE_REMOTE_CONTROL / HIGH_VOLTAGE_REMOTE_CONTROL (the voltage of the remote
    // controller's own bus) ; unused
    // (NaN) for NOT_SIMULATED / DIVERGENCE
    real_type value;
    // limit that was violated. For LOW_Q / HIGH_Q the SUMMED capability of the machines
    // holding that bus, not one machine's: min_q_mvar / max_q_mvar for a generator or an
    // hvdc converter station, b_min / b_max at the solved voltage for a voltage-mode SVC.
    // For HIGH_P the pmax of the direction `side` names (HVDC) or the generator's / storage
    // unit's max_p_mw; for LOW_P its min_p_mw. For LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q
    // the generator's target voltage, in kV of the regulated bus; for LOW_VOLTAGE_SVC_STANDBY
    // / HIGH_VOLTAGE_SVC_STANDBY the automaton's low / high threshold, in kV of the
    // regulated bus; for LOW_VOLTAGE_REMOTE_CONTROL / HIGH_VOLTAGE_REMOTE_CONTROL the realistic
    // voltage bound, in kV of the controller's own bus. Unused (NaN) for NOT_SIMULATED / DIVERGENCE
    real_type limit;
    // element name: LINE / TRAFO / HVDC / GENERATOR / STORAGE / SVC (from
    // LSGrid::set_line_names / set_trafo_names / set_dcline_names / set_gen_names /
    // set_storage_names / set_svc_names) or, for BUS,
    // the name of the *substation* the violating bus belongs
    // to (from LSGrid::set_substation_names) -- there is no per-bus name in LSGrid, only
    // per-substation ones. Empty string if the grid never had names set for the relevant
    // kind, or for GRID.
    std::string name{};

    /// see ViolationCategory: what kind of statement this violation is
    ViolationCategory category() const noexcept { return violation_category(violation_type); }
};

} // namespace ls2g

#endif  // LIMITVIOLATION_H
