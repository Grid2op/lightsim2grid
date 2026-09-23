// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LS2G_SLACK_REDISTRIBUTION_H
#define LS2G_SLACK_REDISTRIBUTION_H

#include <cmath>
#include <limits>
#include <vector>

#include "BaseConstants.hpp"
#include "TaggedIdVec.hpp"
#include "Utils.hpp"

namespace ls2g {

/**
 * OpenLoadFlow-style bounded redistribution of a KNOWN active-power imbalance
 * (OLF's `DistributedSlack` outer loop, `GenerationActivePowerDistributionStep`).
 *
 * The distributed slack of the Newton solve shares whatever imbalance is left by the
 * solution proportionally to fixed per-bus weights, with no limit. When a contingency
 * strands a big generator, that pushes the remaining machines past their `max_p`
 * (OLF stops each one at its bound and re-shares the excess). The pre-pass here does
 * what OLF's loop does, BEFORE the solve, on the part of the imbalance that is known
 * upfront (the setpoints of what the contingency took out): every participating unit
 * gets its share, clamped to `[min_p, max_p]`; a clamped unit leaves the pool and what
 * it could not take is shared again among the others. The caller then writes the new
 * setpoints and takes the saturated units out of the distributed slack, so the Newton
 * solve only shares what is left (the change in the losses) on the units that can
 * still move.
 *
 * Shared by `LSGrid::consider_only_main_component` / `LSGrid::redistribute_active_power`
 * and the batch sweeps (`BaseBatchSweep`, option `redistribute_slack`).
 *
 * Everything is in MW and in GENERATOR convention (a storage unit's target is in load
 * convention: the caller negates it on the way in and on the way out).
 */
namespace slack_redistribution {

enum class UnitKind : char { GENERATOR = 0, STORAGE = 1 };

struct Participant {
    UnitKind kind;
    int el_id;               // id in its own container
    int bus;                 // the caller's own bus numbering (informational)
    real_type injection_mw;  // current setpoint, generator convention
    real_type weight;        // raw (un-normalised) slack weight, > 0
    real_type min_p_mw;      // NaN: unbounded below
    real_type max_p_mw;      // NaN: unbounded above
};

/// What the redistribution did (exposed to python as `SlackRedistributionReport`).
struct Report {
    real_type mismatch_mw = 0.;        // what had to be shared (> 0: units inject more)
    int nb_participants = 0;
    int nb_saturated = 0;              // units that reached a bound (and left the slack)
    int nb_rounds = 0;                 // 0: nothing was shared
    real_type not_distributed_mw = 0.; // what no unit could take (all saturated)
    bool all_saturated = false;        // every unit hit its bound: none left the slack
};

constexpr real_type default_eps_mw = 1e-6;

/**
 * Share `mismatch_mw` (> 0: the units must inject more) on `units`, in their order
 * (the caller lists the generators by id then the storage units by id, so the result
 * is deterministic). Writes each unit's new injection in `new_injection_mw` and whether
 * it reached a bound in `saturated` (both resized to `units.size()`).
 *
 * When EVERY unit hits its bound, `saturated` is cleared (all zero) and
 * `all_saturated` is set: the caller keeps them all in the slack (a solve with no
 * slack is impossible, so a bound will have to give) and `not_distributed_mw` says
 * how much was left.
 */
inline Report distribute(const std::vector<Participant> & units,
                         real_type mismatch_mw,
                         real_type eps_mw,
                         std::vector<real_type> & new_injection_mw,
                         std::vector<char> & saturated)
{
    const std::size_t nb = units.size();
    Report report;
    report.mismatch_mw = mismatch_mw;
    report.nb_participants = static_cast<int>(nb);
    new_injection_mw.resize(nb);
    saturated.assign(nb, 0);
    for(std::size_t k = 0; k < nb; ++k) new_injection_mw[k] = units[k].injection_mw;
    if(nb == 0 || std::abs(mismatch_mw) <= eps_mw) return report;

    const real_type inf = std::numeric_limits<real_type>::infinity();
    std::vector<real_type> lo(nb), hi(nb);
    std::vector<char> active(nb, 1);
    for(std::size_t k = 0; k < nb; ++k){
        lo[k] = std::isfinite(units[k].min_p_mw) ? units[k].min_p_mw : -inf;
        hi[k] = std::isfinite(units[k].max_p_mw) ? units[k].max_p_mw : inf;
    }

    real_type remaining = mismatch_mw;
    std::size_t nb_active = nb;
    // each round either ends the loop or removes at least one unit: nb + 1 rounds is
    // a hard cap, never reached in practice
    while(nb_active > 0 && std::abs(remaining) > eps_mw && report.nb_rounds <= static_cast<int>(nb) + 1){
        ++report.nb_rounds;
        real_type factor_sum = 0.;
        for(std::size_t k = 0; k < nb; ++k) if(active[k]) factor_sum += units[k].weight;
        if(factor_sum <= 0.) break;
        real_type done = 0.;
        for(std::size_t k = 0; k < nb; ++k){
            if(!active[k]) continue;
            const real_type old = new_injection_mw[k];
            real_type cand = old + remaining * units[k].weight / factor_sum;
            if(remaining > 0. && cand >= hi[k]){
                // never DEcrease a unit already above its max (same rule for the min)
                cand = old > hi[k] ? old : hi[k];
                active[k] = 0;
                saturated[k] = 1;
                --nb_active;
            } else if(remaining < 0. && cand <= lo[k]){
                cand = old < lo[k] ? old : lo[k];
                active[k] = 0;
                saturated[k] = 1;
                --nb_active;
            }
            done += cand - old;
            new_injection_mw[k] = cand;
        }
        remaining -= done;
    }
    report.not_distributed_mw = remaining;
    for(std::size_t k = 0; k < nb; ++k) if(saturated[k]) ++report.nb_saturated;
    if(nb_active == 0){
        report.all_saturated = true;
        saturated.assign(nb, 0);
    }
    return report;
}

/**
 * Append the participating units of one family (a GeneratorContainer or a
 * StorageContainer) to `out`: connected, flagged slack, with a positive weight, whose
 * bus `keep_bus(bus_me)` accepts and that `is_off(el_id)` does not exclude. Their
 * injection is `target_sign * injection_of(el_id)` (`target_sign` is -1 for the
 * load-convention storage units), their limits the container's (NaN when unset).
 */
template<class Container, class KeepBus, class IsOff, class InjectionOf>
inline void append_participants(const Container & container,
                                UnitKind kind,
                                real_type target_sign,
                                KeepBus keep_bus,
                                IsOff is_off,
                                InjectionOf injection_of,
                                std::vector<Participant> & out)
{
    const int nb_el = container.nb();
    const std::vector<bool> & status = container.get_status();
    const GlobalBusIdVect & buses = container.get_bus_id();
    for(int el_id = 0; el_id < nb_el; ++el_id){
        if(!status[el_id]) continue;
        if(!container.is_slack(el_id)) continue;
        const real_type weight = container.get_slack_weight(el_id);
        if(weight <= BaseConstants::_tol_equal_float) continue;
        const int bus_me = buses(el_id).cast_int();
        if(bus_me == BaseConstants::_deactivated_bus_id) continue;
        if(!keep_bus(bus_me)) continue;
        if(is_off(el_id)) continue;
        Participant part;
        part.kind = kind;
        part.el_id = el_id;
        part.bus = bus_me;
        part.injection_mw = target_sign * injection_of(el_id);
        part.weight = weight;
        part.min_p_mw = container.get_min_p(el_id);
        part.max_p_mw = container.get_max_p(el_id);
        out.push_back(part);
    }
}

/**
 * `sign * Σ value_of(el_id)` over the connected elements of `container` whose bus
 * `bus_lost(bus_me)` flags: what a family loses (in generator convention when `sign`
 * is the family's own: +1 generators / static generators, -1 loads / storage units /
 * shunts) when those buses leave the solved grid.
 */
template<class Container, class BusLost, class ValueOf>
inline real_type sum_setpoints_if(const Container & container,
                                  real_type sign,
                                  BusLost bus_lost,
                                  ValueOf value_of)
{
    const int nb_el = container.nb();
    const std::vector<bool> & status = container.get_status();
    const GlobalBusIdVect & buses = container.get_bus_id();
    real_type res = 0.;
    for(int el_id = 0; el_id < nb_el; ++el_id){
        if(!status[el_id]) continue;
        const int bus_me = buses(el_id).cast_int();
        if(bus_me == BaseConstants::_deactivated_bus_id) continue;
        if(!bus_lost(bus_me)) continue;
        res += value_of(el_id);
    }
    return sign * res;
}

}  // namespace slack_redistribution
}  // namespace ls2g

#endif  // LS2G_SLACK_REDISTRIBUTION_H
