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
 * As OLF ("we don't want to change the generation sign"), a unit never crosses 0 MW: one
 * that injects stops at 0 when the mismatch pushes it down, one that draws (a charging
 * storage unit, a pumping machine, any unit with `min_p < 0` running below 0) stops at 0
 * when pushed up, whatever its `[min_p, max_p]` says. A storage unit whose whole range
 * straddles 0 is the common case: it is only ever moved on the side of 0 it started on.
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
    /// a participant of the Newton solve's distributed slack (the caller takes it out of
    /// it when it saturates), or only of this pre-pass (a unit flagged "can participate
    /// in the slack", see SlackParticipation::set_can_participate)
    bool in_slack = true;
    /// for a unit only flagged "can participate in the slack": how far BEYOND the limit it
    /// sits at it was in the reference solve, MW, > 0 above its upper limit and < 0 below
    /// its lower one (see SlackParticipation::set_can_participate_overshoot); 0 otherwise
    real_type overshoot_mw = 0.;
};

/// What the redistribution did (exposed to python as `SlackRedistributionReport`).
struct Report {
    real_type mismatch_mw = 0.;        // what had to be shared (> 0: units inject more)
    int nb_participants = 0;
    int nb_saturated = 0;              // units that reached a bound (and left the slack)
    int nb_rounds = 0;                 // 0: nothing was shared
    real_type not_distributed_mw = 0.; // what no unit could take (all saturated)
    bool all_saturated = false;        // every unit of the slack hit its bound: none left it
};

constexpr real_type default_eps_mw = 1e-6;

/**
 * Share `mismatch_mw` (> 0: the units must inject more) on `units`, in their order
 * (the caller lists the generators by id then the storage units by id, so the result
 * is deterministic). Writes each unit's new injection in `new_injection_mw` and whether
 * it reached a bound in `saturated` (both resized to `units.size()`).
 *
 * When every unit of the Newton solve's distributed slack (`in_slack`) hits its bound,
 * `saturated` is cleared (all zero) and `all_saturated` is set: the caller keeps them
 * all in the slack (a solve with no slack is impossible, so a bound will have to give)
 * and `not_distributed_mw` says how much no unit could take. That holds even when a
 * unit only flagged "can participate" still has room: it takes the rest of the
 * mismatch here, but the solve cannot distribute on it.
 */
inline Report distribute_with_overshoot(const std::vector<Participant> & units,
                                        real_type mismatch_mw,
                                        real_type eps_mw,
                                        std::vector<real_type> & new_injection_mw,
                                        std::vector<char> & saturated);

inline Report distribute(const std::vector<Participant> & units,
                         real_type mismatch_mw,
                         real_type eps_mw,
                         std::vector<real_type> & new_injection_mw,
                         std::vector<char> & saturated)
{
    // a unit that sat beyond its limit in the reference solve: the exact rule (see
    // distribute_with_overshoot); without one, this one, which gives the same answer
    for(const Participant & unit : units){
        if(unit.overshoot_mw != 0.){
            return distribute_with_overshoot(units, mismatch_mw, eps_mw, new_injection_mw, saturated);
        }
    }
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
        // the sign of the injection is kept (OLF's GenerationActivePowerDistributionStep):
        // 0 MW is a bound on the side the unit is not on
        if(units[k].injection_mw < 0.) hi[k] = std::min(hi[k], 0.);
        else lo[k] = std::max(lo[k], 0.);
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
    // what is left of the solve's slack: the units still active that are in it (a unit only
    // flagged "can participate" is not, so it cannot keep the slack from being emptied)
    bool any_in_slack = false;
    bool slack_left = false;
    for(std::size_t k = 0; k < nb; ++k){
        if(!units[k].in_slack) continue;
        any_in_slack = true;
        if(active[k]){
            slack_left = true;
            break;
        }
    }
    if(nb_active == 0 || (any_in_slack && !slack_left)){
        report.all_saturated = true;
        saturated.assign(nb, 0);
    }
    return report;
}

/**
 * `distribute`, when some unit carries an overshoot: OpenLoadFlow's rule written as what it
 * is, `p_k = clamp(v_k + delta * w_k / W, lo_k, hi_k)` with ONE common shift `delta` (MW of the
 * whole mismatch) chosen so that the units take `mismatch_mw` between them. `v_k` is where the
 * unit would be without its bounds: its injection, plus its overshoot beyond the limit it sits
 * at -- so a unit OpenLoadFlow capped well beyond its max_p stays there until the shift has
 * used that up, as when OpenLoadFlow shares the slack from the raw set-points -- and `w_k / W`
 * its normalised weight. The bounds are the round-based algorithm's (0 MW on the side the unit
 * is not on), widened to its injection when it already sits beyond one, so that a unit above
 * its max_p is never pulled down to it by a positive mismatch. The total is monotone in the
 * shift: bisection.
 *
 * Same outputs and same "all saturated" convention as `distribute`; `nb_rounds` counts the
 * bisection steps.
 */
inline Report distribute_with_overshoot(const std::vector<Participant> & units,
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
    std::vector<real_type> lo(nb), hi(nb), virt(nb), share(nb);
    real_type weight_sum = 0.;
    for(std::size_t k = 0; k < nb; ++k) weight_sum += units[k].weight;
    if(weight_sum <= 0.) return report;
    for(std::size_t k = 0; k < nb; ++k){
        const real_type inj = units[k].injection_mw;
        const real_type over = units[k].overshoot_mw;
        real_type l = std::isfinite(units[k].min_p_mw) ? units[k].min_p_mw : -inf;
        real_type h = std::isfinite(units[k].max_p_mw) ? units[k].max_p_mw : inf;
        // the side of 0 MW the unit is on: that of its injection, and at 0 MW the one its
        // overshoot gives -- a drawing unit capped at 0 MW sits beyond it from below (> 0)
        const bool drawing = inj < 0. || (inj <= eps_mw && over > 0.);
        if(drawing) h = std::min(h, 0.);
        else l = std::max(l, 0.);
        // the overshoot is beyond the limit its sign says the unit sits at
        real_type offset = 0.;
        if(over > 0. && std::isfinite(h) && inj >= h - eps_mw) offset = over;
        else if(over < 0. && std::isfinite(l) && inj <= l + eps_mw) offset = over;
        lo[k] = std::min(l, inj);
        hi[k] = std::max(h, inj);
        virt[k] = inj + offset;
        share[k] = units[k].weight / weight_sum;
    }
    auto taken = [&](real_type delta){
        real_type res = 0.;
        for(std::size_t k = 0; k < nb; ++k){
            const real_type p = std::min(hi[k], std::max(lo[k], virt[k] + delta * share[k]));
            res += p - units[k].injection_mw;
        }
        return res;
    };
    // what the units can take at most in the direction asked
    real_type capacity = 0.;
    for(std::size_t k = 0; k < nb; ++k){
        capacity += (mismatch_mw > 0.) ? hi[k] - units[k].injection_mw : lo[k] - units[k].injection_mw;
    }
    const bool all_saturated = std::isfinite(capacity) && std::abs(capacity) <= std::abs(mismatch_mw) + eps_mw;
    real_type delta = 0.;
    if(all_saturated){
        // every unit ends at its bound: as far as the shift can push them
        for(std::size_t k = 0; k < nb; ++k){
            new_injection_mw[k] = (mismatch_mw > 0.) ? hi[k] : lo[k];
        }
        report.not_distributed_mw = mismatch_mw - capacity;
        report.nb_rounds = 1;
        // as `distribute`: counted, then left out of the "saturated" mask (none leaves the slack)
        report.nb_saturated = static_cast<int>(nb);
        report.all_saturated = true;
        return report;
    }
    // bracket then bisect: the total is continuous and non decreasing in the shift
    real_type a = 0., b = mismatch_mw;
    int guard = 0;
    while(std::abs(taken(b)) < std::abs(mismatch_mw) && guard < 200){ a = b; b *= 2.; ++guard; }
    for(int it = 0; it < 200; ++it){
        const real_type mid = 0.5 * (a + b);
        const real_type t = taken(mid);
        ++report.nb_rounds;
        if(std::abs(t - mismatch_mw) <= 1e-3 * eps_mw) { a = b = mid; break; }
        if((t < mismatch_mw) == (mismatch_mw > 0.)) a = mid; else b = mid;
    }
    delta = 0.5 * (a + b);
    real_type done = 0.;
    for(std::size_t k = 0; k < nb; ++k){
        const real_type target = virt[k] + delta * share[k];
        const real_type p = std::min(hi[k], std::max(lo[k], target));
        new_injection_mw[k] = p;
        done += p - units[k].injection_mw;
        if((mismatch_mw > 0. && target >= hi[k]) || (mismatch_mw < 0. && target <= lo[k])){
            saturated[k] = 1;
            ++report.nb_saturated;
        }
    }
    report.not_distributed_mw = mismatch_mw - done;
    // as `distribute`: a unit only flagged "can participate" still moving cannot keep the
    // solve's slack from being emptied, every unit of it at its bound keeps them all in it
    bool any_in_slack = false;
    bool slack_left = false;
    for(std::size_t k = 0; k < nb; ++k){
        if(!units[k].in_slack) continue;
        any_in_slack = true;
        if(!saturated[k]){
            slack_left = true;
            break;
        }
    }
    if(any_in_slack && !slack_left){
        report.all_saturated = true;
        saturated.assign(nb, 0);
    }
    return report;
}

/**
 * Append the participating units of one family (a GeneratorContainer or a
 * StorageContainer) to `out`: connected, flagged slack with a positive weight -- or
 * flagged "can participate in the slack" (an outer loop only left it out because it sat
 * at an active limit, see SlackParticipation::set_can_participate), with that weight
 * and `in_slack` false -- whose bus `keep_bus(bus_me)` accepts and that `is_off(el_id)`
 * does not exclude. Their injection is `target_sign * injection_of(el_id)`
 * (`target_sign` is -1 for the load-convention storage units), their limits the
 * container's (NaN when unset): a unit sitting at a limit only ever moves away from it.
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
        const bool in_slack = container.is_slack(el_id)
                              && container.get_slack_weight(el_id) > BaseConstants::_tol_equal_float;
        const real_type weight = in_slack ? container.get_slack_weight(el_id)
                                          : container.get_can_participate_slack_weight(el_id);
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
        part.in_slack = in_slack;
        part.overshoot_mw = in_slack ? 0. : container.get_can_participate_slack_overshoot(el_id);
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

/**
 * Σ of the active setpoints (generator convention) of the converter stations of the
 * active HVDC lines of `hvdc_lines` whose bus `bus_lost(bus_me)` flags: what the solved
 * grid loses when those converters leave it. A line keeps the converter that stays in
 * the main component injecting its scheduled power (see
 * `HvdcLineContainer::_disconnect_if_not_in_main_component`), so only the stranded
 * station(s) count -- a rectifier draws from the grid it leaves (a negative setpoint,
 * the balance loses a consumption), an inverter feeds it (a positive one). A station
 * already open (the far end of a cross-border link, outside the solved grid since the
 * start: its side keeps its bus id) was never in the balance and does not count.
 */
template<class HvdcContainer, class BusLost>
inline real_type sum_hvdc_station_setpoints_if(const HvdcContainer & hvdc_lines,
                                               BusLost bus_lost)
{
    const int nb_el = hvdc_lines.nb();
    const std::vector<bool> & status = hvdc_lines.get_status_global();
    const GlobalBusIdVect & bus_1 = hvdc_lines.get_bus_id_side_1();
    const GlobalBusIdVect & bus_2 = hvdc_lines.get_bus_id_side_2();
    const auto & side_1 = hvdc_lines.get_stations_side_1();
    const auto & side_2 = hvdc_lines.get_stations_side_2();
    const std::vector<bool> & status_1 = side_1.get_status();
    const std::vector<bool> & status_2 = side_2.get_status();
    real_type res = 0.;
    for(int el_id = 0; el_id < nb_el; ++el_id){
        if(!status[el_id]) continue;
        const int b1 = bus_1(el_id).cast_int();
        const int b2 = bus_2(el_id).cast_int();
        if(status_1[el_id] && b1 != BaseConstants::_deactivated_bus_id && bus_lost(b1)) res += side_1.get_target_p(el_id);
        if(status_2[el_id] && b2 != BaseConstants::_deactivated_bus_id && bus_lost(b2)) res += side_2.get_target_p(el_id);
    }
    return res;
}

}  // namespace slack_redistribution
}  // namespace ls2g

#endif  // LS2G_SLACK_REDISTRIBUTION_H
