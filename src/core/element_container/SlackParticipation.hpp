// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LS2G_SLACK_PARTICIPATION_H
#define LS2G_SLACK_PARTICIPATION_H

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Eigen/Core"

#include "BaseConstants.hpp"
#include "SlackRedistribution.hpp"
#include "TaggedIdVec.hpp"
#include "Utils.hpp"

// feature test for code built against this header (eg gpusim2grid): the per-element
// "can participate in the slack" weight (SlackParticipation::set_can_participate)
#define LS2G_HAS_CAN_PARTICIPATE_SLACK 1
// ... and how far beyond its limit each such element was (SlackParticipation::set_can_participate_overshoot)
#define LS2G_HAS_CAN_PARTICIPATE_SLACK_OVERSHOOT 1

namespace ls2g {

/**
 * Which elements of ONE container take part in the distributed slack, and with what
 * (un-normalised) weight.
 *
 * Generators and storage units both carry one (OpenLoadFlow distributes the slack on
 * batteries exactly as on generators), and LSGrid adds the families up: the slack bus
 * set is the union of theirs, a bus' weight is the sum of every participant on it,
 * and the share a bus absorbed is split back onto its participants by raw weight
 * whatever their family. The rule for who participates -- connected, flagged, with a
 * non-zero weight, and not taken out by the caller -- lives here once, so that the
 * weights a solve used and the split of its result can never disagree, and neither
 * can two families.
 *
 * It holds data and rules only: the element status, bus and result vectors stay the
 * container's, and are passed in.
 *
 * TODO (see the CHANGELOG's [TODO]): participating in the distributed slack is a purely
 * ACTIVE statement -- this element takes a share of the power imbalance -- and says
 * nothing about voltage, unlike being the REFERENCE slack (theta known, |V| known, P and
 * Q unknown). Today only a voltage source can be given a share, which does not follow: a
 * load could take one. Separating the two roles is what would allow it.
 */
class SlackParticipation
{
    public:
        void reset(std::size_t nb_el){
            slackbus_.assign(nb_el, false);
            weight_.assign(nb_el, 0.);
            can_participate_weight_.assign(nb_el, 0.);
            can_participate_overshoot_mw_.assign(nb_el, 0.);
        }

        /**
         * "Can participate in the slack": every element that can be in the distributed
         * slack, with the weight it has there (on the same scale as the slack weights, 0
         * for an element that cannot). A slack participant carries it on its own (`add`
         * sets it, `remove` clears it); what this call adds are the elements a caller
         * knows an outer loop only left out of the slack because they sat at an active
         * limit in the reference solve (OpenLoadFlow caps a unit at max_p when the
         * mismatch it distributes is positive, at min_p when it is negative -- so that
         * unit takes no share of a mismatch of that sign, and a full one of the other).
         * The slack participants keep their flag, whatever `flags` says of them.
         *
         * A flagged element out of the slack is not read by the Newton solve (its
         * distributed slack has no bounds, so it would push such a unit past its limit).
         * The bounded redistribution pre-pass (SlackRedistribution.hpp) counts it as a
         * participant, within its [min_p, max_p] -- which lets it move away from the limit
         * it sits at, never across it -- and it goes back into the slack as soon as it is
         * connected with a set-point off its limits (`rejoin_if_able`). A flagged element
         * that is also a slack participant takes part with its slack weight.
         */
        void set_can_participate(const std::vector<bool> & flags,
                                 const Eigen::Ref<const RealVect> & weights,
                                 const char * fun_name){
            const std::size_t nb_el = slackbus_.size();
            if(flags.size() != nb_el || static_cast<std::size_t>(weights.size()) != nb_el){
                std::ostringstream exc_;
                exc_ << fun_name << ": expected " << nb_el << " flags and weights, got "
                     << flags.size() << " and " << weights.size() << ".";
                throw std::runtime_error(exc_.str());
            }
            std::vector<real_type> res(nb_el, 0.);
            for(std::size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!flags[el_id]) continue;
                const real_type w = weights(static_cast<Eigen::Index>(el_id));
                if(!std::isfinite(w) || w <= 0.){
                    std::ostringstream exc_;
                    exc_ << fun_name << ": the element with id " << el_id
                         << " is flagged but its weight is not a finite, positive number (got " << w << ").";
                    throw std::runtime_error(exc_.str());
                }
                res[el_id] = w;
            }
            for(std::size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!flags[el_id] && slackbus_[el_id]) res[el_id] = weight_[el_id];
            }
            // nothing a powerflow reads: no AlgoControl flag to raise
            can_participate_weight_ = res;
        }
        [[nodiscard]] bool can_participate(int el_id) const {
            return can_participate_weight_[el_id] > 0.;
        }
        [[nodiscard]] real_type can_participate_weight(int el_id) const {return can_participate_weight_[el_id];}
        [[nodiscard]] const std::vector<real_type> & can_participate_weights() const {return can_participate_weight_;}
        /// restore a serialized state (sizes are checked by the caller). A slack participant
        /// of a state saved before it carried its own flag gets it.
        void set_can_participate_weights(const std::vector<real_type> & weights){
            can_participate_weight_ = weights;
            for(std::size_t el_id = 0; el_id < slackbus_.size(); ++el_id){
                if(slackbus_[el_id] && !(can_participate_weight_[el_id] > 0.)) can_participate_weight_[el_id] = weight_[el_id];
            }
        }

        /**
         * How far BEYOND the limit it sits at an element flagged "can participate in the
         * slack" was in the reference solve, in MW (>= 0, 0 by default): OpenLoadFlow shares
         * the slack from the raw set-points, `p = clamp(raw + lambda * weight)`, so a unit it
         * capped at max_p had `raw + lambda * weight` above max_p by that much (below min_p
         * for a unit capped at min_p). A later imbalance of the other sign only moves it once
         * the common shift of the distribution has used that up -- before, OpenLoadFlow
         * still caps it. 0 means it leaves its limit at once (an element whose reference
         * solve sat exactly at it).
         *
         * Only the bounded redistribution pre-pass (SlackRedistribution.hpp) reads it, and
         * only for an element flagged "can participate" (see set_can_participate).
         */
        void set_can_participate_overshoot(const Eigen::Ref<const RealVect> & overshoot_mw,
                                           const char * fun_name){
            const std::size_t nb_el = slackbus_.size();
            if(static_cast<std::size_t>(overshoot_mw.size()) != nb_el){
                std::ostringstream exc_;
                exc_ << fun_name << ": expected " << nb_el << " values, got " << overshoot_mw.size() << ".";
                throw std::runtime_error(exc_.str());
            }
            std::vector<real_type> res(nb_el, 0.);
            for(std::size_t el_id = 0; el_id < nb_el; ++el_id){
                const real_type o = overshoot_mw(static_cast<Eigen::Index>(el_id));
                if(!std::isfinite(o) || o < 0.){
                    std::ostringstream exc_;
                    exc_ << fun_name << ": the element with id " << el_id
                         << " has an overshoot that is not a finite, non negative number (got " << o << ").";
                    throw std::runtime_error(exc_.str());
                }
                res[el_id] = o;
            }
            // nothing a powerflow reads: no AlgoControl flag to raise
            can_participate_overshoot_mw_ = res;
        }
        [[nodiscard]] real_type can_participate_overshoot(int el_id) const {return can_participate_overshoot_mw_[el_id];}
        [[nodiscard]] const std::vector<real_type> & can_participate_overshoots() const {return can_participate_overshoot_mw_;}
        /// restore a serialized state (sizes are checked by the caller)
        void set_can_participate_overshoots(const std::vector<real_type> & overshoot_mw){
            can_participate_overshoot_mw_ = overshoot_mw;
        }

        [[nodiscard]] bool is_slack(int el_id) const {return slackbus_[el_id];}
        [[nodiscard]] real_type weight(int el_id) const {return weight_[el_id];}
        /// a non-zero weight: an element carrying one is not "pseudo off" (see GeneratorContainer)
        [[nodiscard]] bool has_weight(int el_id) const {
            return std::abs(weight_[el_id]) >= BaseConstants::_tol_equal_float;
        }
        [[nodiscard]] const std::vector<bool> & flags() const {return slackbus_;}
        [[nodiscard]] const std::vector<real_type> & weights() const {return weight_;}

        /// restore a serialized state (sizes are checked by the caller)
        void set(const std::vector<bool> & flags, const std::vector<real_type> & weights){
            slackbus_ = flags;
            weight_ = weights;
        }

        /**
         * Make `el_id` a participant with the given (> 0) weight, or update its weight.
         *
         * Why `tell_slack_participate_changed()` and nothing about voltage: taking the
         * slack role is an ACTIVE-power role, but the pv/pq split reads the slack set
         * directly -- `fillpv` skips a bus that is in `slack_bus_id_solver` ("slack bus
         * is not PV"), and the PQ loop right after it skips it too -- so a bus joining
         * or leaving the slack set moves the split, whatever it does to voltage. That
         * flag is what says so, and it is measured, not assumed: see the `[pv_pq]` cases
         * in test_cache_reuse.cpp, which reach a state where it is the only term raised
         * and fail if it is dropped.
         *
         * The voltage side needs no flag of its own for the same reason. It is real but
         * secondary -- `GeneratorContainer::is_pseudo_off()` answers false for a slack
         * generator whatever its active power, so a zero-P slack generator IS a voltage
         * controller where an ordinary one would not be -- and the voltage-control plan
         * is built around the split this same flag rebuilds. See
         * AlgoControl::need_recompute_voltage_control.
         */
        void add(int el_id, real_type weight, DualAlgoControl & solver_control, const char * fun_name){
            if(weight <= 0.){
                std::ostringstream exc_;
                exc_ << fun_name << " Cannot assign a negative (<=0) weight to the slack bus.";
                throw std::runtime_error(exc_.str());
            }
            if(!slackbus_[el_id]){ solver_control.tell_slack_participate_changed(); }
            slackbus_[el_id] = true;
            if(std::abs(weight_[el_id] - weight) > BaseConstants::_tol_equal_float){
                solver_control.tell_slack_weight_changed();
                weight_[el_id] = weight;
            }
            // it can be in the slack: it can come back to it (see leave / rejoin_if_able)
            can_participate_weight_[el_id] = weight;
        }
        /// take `el_id` out of the slack for good: it can no longer participate either
        void remove(int el_id, DualAlgoControl & solver_control){
            _take_out(el_id, solver_control);
            can_participate_weight_[el_id] = 0.;
        }
        /**
         * Take `el_id` out of the slack while it sits at the active limit it saturated
         * at (LSGrid::redistribute_active_power), keeping it flagged "can participate"
         * with its weight -- what a unit the bake capped there looks like -- so that
         * `rejoin_if_able` puts it back once it is moved off that limit. (A unit that is
         * merely disconnected does not leave the slack: see `append_slack_buses`.)
         */
        void leave(int el_id, DualAlgoControl & solver_control){
            if(slackbus_[el_id] && weight_[el_id] > 0.) can_participate_weight_[el_id] = weight_[el_id];
            _take_out(el_id, solver_control);
        }
        /**
         * Put back into the slack an element that can participate in it and is out of it
         * only because it sat at an active limit (saturated by the pre-pass, or capped by
         * the bake): `connected` (the container's status, on a bus) and its injection (MW,
         * generator convention) strictly inside its [min_p, max_p] (NaN: no limit on that
         * side). A unit ON a limit, up to the pre-pass tolerance, stays out: that is where
         * the pre-pass leaves a saturated one. Returns whether it rejoined.
         */
        bool rejoin_if_able(int el_id, bool connected, real_type injection_mw,
                            real_type min_p_mw, real_type max_p_mw,
                            DualAlgoControl & solver_control, const char * fun_name){
            if(slackbus_[el_id] || !connected) return false;
            const real_type w = can_participate_weight_[el_id];
            if(!(w > 0.)) return false;
            const real_type eps = slack_redistribution::default_eps_mw;
            if(std::isfinite(max_p_mw) && injection_mw >= max_p_mw - eps) return false;
            if(std::isfinite(min_p_mw) && injection_mw <= min_p_mw + eps) return false;
            add(el_id, w, solver_control, fun_name);
            return true;
        }
        void remove_all(){
            DualAlgoControl unused_solver_control;
            const int nb_el = static_cast<int>(slackbus_.size());
            for(int el_id = 0; el_id < nb_el; ++el_id) remove(el_id, unused_solver_control);
        }

        /// does `el_id` take a share of the slack right now (`status` is the container's)?
        [[nodiscard]] bool participates(int el_id, const std::vector<bool> & status) const {
            return status[el_id] && slackbus_[el_id] && has_weight(el_id);
        }

        /**
         * Add every participant's raw weight to its solver bus in `res` (sized by the
         * number of solver buses). `off`, when non-null, is a nb()-sized mask of
         * elements to leave out on top of the participation rule -- what a batch sweep
         * uses for a row whose contingency disconnects a participant.
         */
        void accumulate_raw(RealVect & res,
                            const std::vector<bool> & status,
                            const GlobalBusIdVect & bus_id,
                            const SolverBusIdVect & id_grid_to_solver,
                            const std::vector<bool> * off,
                            const char * element_name) const
        {
            const int nb_el = static_cast<int>(slackbus_.size());
            for(int el_id = 0; el_id < nb_el; ++el_id){
                if(!participates(el_id, status)) continue;
                if(off != nullptr && (*off)[el_id]) continue;
                const int bus_me = bus_id(el_id).cast_int();
                const int bus_solver = bus_me == BaseConstants::_deactivated_bus_id ?
                                       BaseConstants::_deactivated_bus_id :
                                       id_grid_to_solver[bus_me].cast_int();
                if(bus_solver == BaseConstants::_deactivated_bus_id){
                    // TODO DEBUG MODE: only check in debug mode
                    std::ostringstream exc_;
                    exc_ << "LSGrid::get_slack_weights_solver: the " << element_name << " with id " << el_id
                         << " is connected to a disconnected bus while being connected to the grid.";
                    throw std::runtime_error(exc_.str());
                }
                res.coeffRef(bus_solver) += weight_[el_id];
            }
        }

        /// append the grid buses of the flagged elements (connected or not) that are not in
        /// `buses` yet. A disconnected participant stays flagged and takes no share until it
        /// is reconnected (see `participates`); its bus stays a slack bus as long as it is in
        /// the grid (LSGrid::_slack_bus_id_me drops the ones that are not).
        void append_slack_buses(std::vector<int> & buses, const GlobalBusIdVect & bus_id) const {
            const int nb_el = static_cast<int>(slackbus_.size());
            for(int el_id = 0; el_id < nb_el; ++el_id){
                if(!slackbus_[el_id]) continue;
                const int bus_me = bus_id(el_id).cast_int();
                if(bus_me == BaseConstants::_deactivated_bus_id) continue;
                bool already_there = false;
                for(int b : buses) if(b == bus_me) { already_there = true; break; }
                if(!already_there) buses.push_back(bus_me);
            }
        }

        /**
         * Split the active power each slack bus absorbed (`node_mismatch`, MW, per solver
         * bus) onto its participants, proportionally to their raw weight over
         * `bus_raw_total` -- the raw weight of EVERY participant of that bus, all
         * families included. `sign` is +1 for a container in generator convention, -1
         * for one in load convention.
         */
        void split(RealVect & res_p,
                   real_type sign,
                   const Eigen::Ref<const RealVect> & node_mismatch,
                   const Eigen::Ref<const RealVect> & bus_raw_total,
                   const std::vector<bool> & status,
                   const GlobalBusIdVect & bus_id,
                   const SolverBusIdVect & id_grid_to_solver,
                   const char * fun_name) const
        {
            if(bus_raw_total.size() != node_mismatch.size()){
                // TODO DEBUG MODE: perform this check only in debug mode
                std::ostringstream exc_;
                exc_ << fun_name << ": Impossible to set the active value of the slack participants: no known slack "
                     << "(the slack weights of this solve were never computed).";
                throw std::runtime_error(exc_.str());
            }
            const int nb_el = static_cast<int>(slackbus_.size());
            for(int el_id = 0; el_id < nb_el; ++el_id){
                if(!participates(el_id, status)) continue;
                const int bus_solver = id_grid_to_solver[bus_id(el_id).cast_int()].cast_int();
                // TODO DEBUG MODE: check bus_solver >= 0 and bus_raw_total[bus_solver] > 0
                res_p(el_id) += sign * node_mismatch(bus_solver) * weight_[el_id] / bus_raw_total(bus_solver);
            }
        }

        /// a flagged element must carry a finite, strictly positive weight
        void check_weights(const char * element_name) const {
            const int nb_el = static_cast<int>(slackbus_.size());
            for(int el_id = 0; el_id < nb_el; ++el_id){
                if(!slackbus_[el_id]) continue;
                const real_type w = weight_[el_id];
                if((!std::isfinite(w)) || (w <= BaseConstants::_tol_equal_float)){
                    std::ostringstream exc_;
                    exc_ << "LSGrid::check_grid: " << element_name << " id " << el_id
                         << " is flagged as a slack but has a non-positive or non-finite slack weight ("
                         << w << ").";
                    throw std::runtime_error(exc_.str());
                }
            }
        }

        /// is any element flagged, and is any flagged element connected?
        void summary(const std::vector<bool> & status, bool & any_flagged, bool & any_connected) const {
            const int nb_el = static_cast<int>(slackbus_.size());
            for(int el_id = 0; el_id < nb_el; ++el_id){
                if(!slackbus_[el_id]) continue;
                any_flagged = true;
                if(status[el_id]) any_connected = true;
            }
        }

    private:
        void _take_out(int el_id, DualAlgoControl & solver_control){
            if(slackbus_[el_id]){ solver_control.tell_slack_participate_changed(); }
            if(std::abs(weight_[el_id]) > BaseConstants::_tol_equal_float){ solver_control.tell_slack_weight_changed(); }
            slackbus_[el_id] = false;
            weight_[el_id] = 0.;
        }

        std::vector<bool> slackbus_;     // is this element flagged a slack participant
        std::vector<real_type> weight_;  // its raw weight (does not sum to 1)
        // the weight it has when it can be in the slack, 0 if it cannot: its slack weight
        // for a participant, what the pre-pass uses for one out of it (see set_can_participate)
        std::vector<real_type> can_participate_weight_;
        // how far beyond its limit it was in the reference solve, MW (see set_can_participate_overshoot)
        std::vector<real_type> can_participate_overshoot_mw_;
};

} // namespace ls2g

#endif // LS2G_SLACK_PARTICIPATION_H
