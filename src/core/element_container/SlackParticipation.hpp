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
#include "TaggedIdVec.hpp"
#include "Utils.hpp"

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
        }
        void remove(int el_id, DualAlgoControl & solver_control){
            if(slackbus_[el_id]){ solver_control.tell_slack_participate_changed(); }
            if(std::abs(weight_[el_id]) > BaseConstants::_tol_equal_float){ solver_control.tell_slack_weight_changed(); }
            slackbus_[el_id] = false;
            weight_[el_id] = 0.;
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

        /// append the grid buses of the flagged elements (connected or not) that are not in `buses` yet
        void append_slack_buses(std::vector<int> & buses, const GlobalBusIdVect & bus_id) const {
            const int nb_el = static_cast<int>(slackbus_.size());
            for(int el_id = 0; el_id < nb_el; ++el_id){
                if(!slackbus_[el_id]) continue;
                const int bus_me = bus_id(el_id).cast_int();
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
        std::vector<bool> slackbus_;     // is this element flagged a slack participant
        std::vector<real_type> weight_;  // its raw weight (does not sum to 1)
};

} // namespace ls2g

#endif // LS2G_SLACK_PARTICIPATION_H
