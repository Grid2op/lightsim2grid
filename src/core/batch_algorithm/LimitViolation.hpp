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
    HVDC = 4,  // an hvdc line, by its own id (see LimitViolationType::HIGH_P)
    // a generator, by its own id. Its ACTIVE power (LOW_P / HIGH_P) and, for a PQ machine
    // pinned at a reactive limit, the voltage of the bus it would regulate
    // (LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q). A REGULATING machine's reactive
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
    STORAGE = 6
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
    // below limit), and downstream code relies on it. See GenPvReleaseCheck.hpp.
    LOW_VOLTAGE_AT_MIN_Q = 9,
    // ... and the mirror: pinned at its MAXIMUM, regulated bus ABOVE the target (value
    // above limit).
    HIGH_VOLTAGE_AT_MAX_Q = 10
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
    /// LOW_Q, HIGH_Q, LOW_P, HIGH_P, LOW_VOLTAGE_AT_MIN_Q, HIGH_VOLTAGE_AT_MAX_Q.
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
            return ViolationCategory::PHYSICAL;
        default:  // NOT_SIMULATED, DIVERGENCE
            return ViolationCategory::SOLVER;
    }
}

// a single limit violation, as detected by ContingencyAnalysis (see compute_limit_violations)
// and by the physical-limit checks (see compute_physical_violations)
struct LS2G_API LimitViolation {
    ViolationElementType element_type;
    // grid-model bus id for BUS ; local (0-based, own type) line / trafo / hvdc line /
    // generator / storage id otherwise ; unused (-1) for GRID
    int element_id;
    // 1 or 2 for LINE / TRAFO (the terminal) and for HVDC (the direction: 1 means the flow
    // leaves side 1, ie 1 -> 2) ; unused (0) for BUS / GENERATOR / STORAGE / GRID
    int side;
    LimitViolationType violation_type;
    // value reached: MVAr for LOW_Q / HIGH_Q (what the machines holding that bus had to
    // produce), MW for LOW_P / HIGH_P (an hvdc line's flow leaving `side`, always positive;
    // a generator's or a storage unit's own converged active power, signed and in the
    // generator convention), kV for LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q (the voltage
    // of the bus the pinned generator would regulate) ; unused (NaN) for NOT_SIMULATED /
    // DIVERGENCE
    real_type value;
    // limit that was violated. For LOW_Q / HIGH_Q the SUMMED capability of the machines
    // holding that bus, not one machine's: min_q_mvar / max_q_mvar for a generator or an
    // hvdc converter station, b_min / b_max at the solved voltage for a voltage-mode SVC.
    // For HIGH_P the pmax of the direction `side` names (HVDC) or the generator's / storage
    // unit's max_p_mw; for LOW_P its min_p_mw. For LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q
    // the generator's target voltage, in kV of the regulated bus. Unused (NaN) for
    // NOT_SIMULATED / DIVERGENCE
    real_type limit;
    // element name: LINE / TRAFO / HVDC / GENERATOR / STORAGE (from LSGrid::set_line_names /
    // set_trafo_names / set_dcline_names / set_gen_names / set_storage_names) or, for BUS,
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
