// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef VOLTAGE_SOURCE_CONTAINER_H
#define VOLTAGE_SOURCE_CONTAINER_H

#include <cmath>
#include <sstream>

#include "Utils.hpp"
#include "OneSideContainer_PQ.hpp"

namespace ls2g {

/**
 * A one-sided element that can regulate a voltage: a generator, an HVDC
 * converter station, an SVC -- and, later, a STATCOM.
 *
 * On top of the setpoints of OneSideContainer_PQ it owns the three facts that
 * make an element a voltage source: whether it regulates at all
 * (`voltage_regulator_on_`), what magnitude it holds (`target_vm_pu_`) and at
 * which bus (`regulated_bus_id_`; its own bus for local control, which is what
 * a bus change follows). Everything the solver asks a voltage source is written
 * once here and read by every leaf: the voltage initialisation (`set_vm`), the
 * classical PV path (`_fillpv`), the reactive output that needs no solve
 * (`set_q`, `_compute_res_pq`), the controller predicates the voltage-control
 * plan is built from, and the flags a status, bus or setpoint change raises.
 *
 * It is a CRTP base rather than a virtual one for the reason spelled out on
 * GenericContainer: the loops here run over every element at every solve
 * (`set_vm`) or every rebuild (`_fillpv`, `set_q`), and the one thing a leaf
 * adds to them is a per-element predicate. Reaching it through `Leaf`
 * statically keeps that predicate inlined; a virtual call per element would be
 * the first one on those paths. The leaf provides:
 *
 *   bool _treated_as_off(int el_id) const;
 *       an element that is connected and regulating but must be treated as
 *       switched off anyway -- a generator at 0 MW when `turnedoff_no_pv` is set.
 *       False for a station or an SVC.
 *   bool _set_vm_skips(int el_id) const;
 *       a regulating element `set_vm` must NOT pin at its setpoint -- a sloped
 *       SVC, whose regulated bus does not sit at the setpoint. False otherwise.
 *   static constexpr bool set_vm_throws_on_unresolved;
 *       whether `set_vm` throws when the regulated bus is not in the solved
 *       system (a generator: yes, the grid is inconsistent) or skips it (an SVC
 *       may point at a bus outside the main component, and that is legitimate).
 *   static real_type _vm_scale(real_type target_vm, real_type current_vm);
 *       the factor `set_vm` multiplies V by to bring |V| from `current_vm` to
 *       `target_vm`. Mathematically `target_vm / current_vm` for everyone; the
 *       generators (and stations) have always computed it as
 *       `(1 / current_vm) * target_vm`, the SVCs as `target_vm / current_vm`, and
 *       the two differ in the last bit. Each leaf keeps its own so that every
 *       solve starts from exactly the voltage it started from before -- the
 *       A/B harness of benchmarks/cache_profiling compares answers to the ulp.
 *   static const char * _element_name();
 *       the noun of the error messages ("generator", ...).
 *
 * and, being a leaf, its own `init` / state, `_fillSbus` (the injection differs
 * per element), and whatever it adds to the `_on_xxx` flags (a generator's
 * slack terms, an SVC's `one_el_changed_bus`).
 */
template<class Leaf>
class VoltageSourceContainer : public OneSideContainer_PQ
{
    public:
        VoltageSourceContainer() noexcept = default;
        ~VoltageSourceContainer() noexcept override = default;

        // ---- what regulates what -----------------------------------------------
        bool get_voltage_regulator_on(int el_id) const {return voltage_regulator_on_[el_id];}
        real_type get_target_vm_pu(int el_id) const {return target_vm_pu_(el_id);}
        /// grid bus id whose voltage this element regulates (== its own bus for local control)
        int get_regulated_bus_id(int el_id) const {return regulated_bus_id_(el_id);}
        bool regulates_remote(int el_id) const {
            return regulated_bus_id_(el_id) != bus_id_(el_id).cast_int();
        }
        /// an ACTIVE voltage regulator: connected, regulating, and not treated as off
        bool is_voltage_controller(int el_id) const {
            if(!status_[el_id]) return false;
            if(!voltage_regulator_on_[el_id]) return false;
            if(leaf()._treated_as_off(el_id)) return false;
            return true;
        }
        /// ... pinning the magnitude of its OWN bus (the classical PV path). A bus
        /// with such an element is Vm-fixed; a bus without one is PQ-for-voltage
        /// (free Vm + Q equation), even when it also carries the slack role.
        bool is_local_voltage_controller(int el_id) const {
            return is_voltage_controller(el_id) && !regulates_remote(el_id);
        }
        /// ... regulating ANOTHER bus: a controller of a VoltageControl group
        /// instead of the PV path (its own bus stays PQ)
        bool is_remote_voltage_controller(int el_id) const {
            return is_voltage_controller(el_id) && regulates_remote(el_id);
        }

        /**
         * Is this element one of those that share the reactive residual of their bus?
         * The complement of what set_q() publishes, minus the ones the algorithm solved
         * for (`solved_by_algo`, built once per solve by LSGrid from the controller
         * list). Asked per element by LSGrid::compute_results, which then groups them
         * by bus.
         */
        bool takes_q_residual_share(int el_id, const std::vector<bool> & solved_by_algo) const
        {
            if(!is_voltage_controller(el_id)) return false;
            // the algorithm solved this one: LSGrid writes it back from the controller list
            if(_is_solved_by_algo(solved_by_algo, el_id)) return false;
            return true;
        }

        /// write the converged reactive output (MVAr) of an element the algorithm
        /// solved for, supplied by the VoltageControl extension (LSGrid::compute_results)
        void set_voltage_control_q(int el_id, real_type q_mvar) {res_q_(el_id) = q_mvar;}

        // ---- setpoints ------------------------------------------------------------
        void set_regulated_bus(int el_id, int bus_id, DualAlgoControl & solver_control){
            // el_id indexes regulated_bus_id_ with an unchecked Eigen operator() below
            // (OOB write for an out-of-range / negative id). bus_id itself is validated
            // by the caller (LSGrid::set_gen_regulated_bus) against the grid bus count.
            _check_in_range(el_id, regulated_bus_id_, "set_regulated_bus");
            if(regulated_bus_id_(el_id) != bus_id){
                regulated_bus_id_(el_id) = bus_id;
                // Only an element that actually regulates reads this field: it decides
                // whether it pins its own bus (PV) or joins a control group at another
                // one, which is the pv/pq split and, through it, the whole
                // voltage-control plan. For a non-regulating element the value is
                // stored and nothing else changes.
                //
                // No tell_recompute_sbus(): fillSbus never reads regulated_bus_id_, and
                // what it does read of a regulating element -- target_p only, the
                // reactive being free -- is the same either way.
                if(voltage_regulator_on_[el_id]) solver_control.tell_pv_changed();
            }
        }

        /// the voltage setpoint, in pu; refuses a disconnected element
        void change_v(int el_id, real_type new_v_pu, DualAlgoControl & solver_control)
        {
            _check_in_range(el_id, status_, "change_v");
            if(!status_[el_id])
            {
                std::ostringstream exc_;
                exc_ << "VoltageSourceContainer::change_v: Impossible to change the voltage setpoint of a disconnected "
                     << Leaf::_element_name() << " (check id " << el_id << ")";
                throw std::runtime_error(exc_.str());
            }
            change_v_nothrow(el_id, new_v_pu, solver_control);
        }
        void change_v_nothrow(int el_id, real_type new_v_pu, DualAlgoControl & solver_control)
        {
            _check_in_range(el_id, status_, "change_v");
            if (std::abs(target_vm_pu_(el_id) - new_v_pu) > _tol_equal_float)
            {
                solver_control.tell_v_changed();
                // A setpoint is the one input of the voltage-control plan the pv/pq split does
                // not read: moving it changes no bus' class at all (the regulated bus of a
                // group stays PQ), so nothing else would notice. Only for an element that
                // regulates -- target_vm_pu_ of one that does not is never read.
                if(voltage_regulator_on_[el_id]) solver_control.ac_algo_controler().tell_voltage_control_changed();
                target_vm_pu_(el_id) = new_v_pu;
            }
        }

        // ---- the solver's questions -----------------------------------------------
        /**
        this functions makes sure that the voltage magnitude of every connected bus is properly used to initialize
        the ac powerflow
        **/
        void set_vm(Eigen::Ref<CplxVect> V, const SolverBusIdVect & id_grid_to_solver) const
        {
            _set_vm_impl(V, id_grid_to_solver, target_vm_pu_);
        }

        /**
        same as set_vm(V, id_grid_to_solver) above, but reads the per-element target
        vm_pu from `target_vm_pu_row` instead of the member `target_vm_pu_` -- used by
        BaseBatchSweep to re-seed |V| with a per-step (per-scenario) generator voltage
        setpoint (see modify_gen_v) rather than the grid's own, fixed, target. Same
        "last writer wins" behaviour as the other overload when several elements
        regulate the same bus (the loop below applies whichever one it visits last).
        **/
        void set_vm(Eigen::Ref<CplxVect> V, const SolverBusIdVect & id_grid_to_solver,
                    const Eigen::Ref<const RealVect> & target_vm_pu_row) const
        {
            if(target_vm_pu_row.size() != static_cast<Eigen::Index>(nb())){
                std::ostringstream exc_;
                exc_ << "VoltageSourceContainer::set_vm: got a target_vm_pu vector of size " << target_vm_pu_row.size()
                     << ", expected " << nb() << " (one entry per " << Leaf::_element_name() << ").";
                throw std::runtime_error(exc_.str());
            }
            _set_vm_impl(V, id_grid_to_solver, target_vm_pu_row);
        }

        /**
         * Publish the reactive output of the elements whose value is known without
         * looking at the powerflow (disconnected, non-regulating, treated as off).
         * A voltage-regulating one is left untouched: its value is either written back
         * by LSGrid from the algorithm's controller list, or it is a share of its bus'
         * reactive residual -- which depends on the other elements of that bus and so
         * is LSGrid's to compute. See takes_q_residual_share.
         */
        void set_q(bool ac)
        {
            const int nb_el = nb();
            if(!ac){
                // do not consider Q values in dc mode
                for(int el_id = 0; el_id < nb_el; ++el_id) res_q_(el_id) = 0.;
                return;
            }
            for(int el_id = 0; el_id < nb_el; ++el_id)
            {
                if(!status_[el_id]){
                    res_q_(el_id) = 0.;  // disconnected
                    continue;
                }
                if (!voltage_regulator_on_[el_id]){
                    // purposedly not pv, so output MVAr = input MVAr (just like a load)
                    res_q_(el_id) = target_q_mvar_(el_id);
                    continue;
                }
                if (leaf()._treated_as_off(el_id)) {
                    // it's as if the element were turned off
                    res_q_(el_id) = 0.;
                    continue;
                }
                // a voltage-regulating element: left for LSGrid, see the note above
            }
        }

    protected:
        const Leaf & leaf() const { return static_cast<const Leaf &>(*this); }

        // shared body of both set_vm() overloads -- `target_vm` is either the member
        // target_vm_pu_ or a caller-supplied per-element vector
        void _set_vm_impl(Eigen::Ref<CplxVect> V, const SolverBusIdVect & id_grid_to_solver,
                          const Eigen::Ref<const RealVect> & target_vm) const
        {
            const int nb_el = nb();
            for(int el_id = 0; el_id < nb_el; ++el_id){
                //  i don't do anything if the element is disconnected
                if(!status_[el_id]) continue;
                if (!voltage_regulator_on_[el_id]) continue;  // purposedly not pv
                if (leaf()._treated_as_off(el_id)) continue;  // in this case turned off elements are not pv
                if (leaf()._set_vm_skips(el_id)) continue;

                // a remote-regulating element sets the magnitude of the REGULATED bus (init quality)
                const int target_grid_bus = regulated_bus_id_(el_id);
                if(target_grid_bus == _deactivated_bus_id){
                    if(!Leaf::set_vm_throws_on_unresolved) continue;
                    std::ostringstream exc_;
                    exc_ << "VoltageSourceContainer::set_vm: " << Leaf::_element_name() << " with id " << el_id
                         << " regulates a disconnected bus while being connected to the grid.";
                    throw std::runtime_error(exc_.str());
                }
                const SolverBusId bus_id_solver = id_grid_to_solver[target_grid_bus];
                if(bus_id_solver.cast_int() == _deactivated_bus_id){
                    if(!Leaf::set_vm_throws_on_unresolved) continue;
                    std::ostringstream exc_;
                    exc_ << "VoltageSourceContainer::set_vm: " << Leaf::_element_name() << " with id " << el_id
                         << " is connected to a disconnected bus while being connected to the grid.";
                    throw std::runtime_error(exc_.str());
                }
                _scale_v_to_target(V, bus_id_solver, target_vm(el_id));
            }
        }

        // scale V at `bus` so that its magnitude is `target_vm` (a zero magnitude is
        // first set to 1., otherwise it would stay 0. whatever the scaling)
        static void _scale_v_to_target(Eigen::Ref<CplxVect> V, SolverBusId bus, real_type target_vm)
        {
            real_type tmp = std::abs(V(bus.cast_int()));
            if(std::abs(tmp) < _tol_equal_float)
            {
                V(bus.cast_int()) = 1.0;
                tmp = 1.0;
            }
            V(bus.cast_int()) *= Leaf::_vm_scale(target_vm, tmp);
        }

        // ---- the hooks ------------------------------------------------------------
        // The classical PV path: every element pinning its OWN bus. (An SVC never
        // goes through it -- it is always a group controller -- and overrides this
        // with nothing.)
        void _fillpv(std::vector<int> & bus_pv,
                     std::vector<bool> & has_bus_been_added,
                     const SolverBusIdVect & slack_bus_id_solver,
                     const SolverBusIdVect & id_grid_to_solver) const override
        {
            const int nb_el = nb();
            for(int el_id = 0; el_id < nb_el; ++el_id){
                if(!is_local_voltage_controller(el_id)) continue;
                const SolverBusId bus_id_solver = _solver_bus(el_id, bus_id_(el_id), id_grid_to_solver, "VoltageSourceContainer::fillpv");
                if(is_in_vect(bus_id_solver.cast_int(), slack_bus_id_solver.to_int_vector())) continue;  // slack bus is not PV
                if(has_bus_been_added[bus_id_solver.cast_int()]) continue; // i already added this bus
                bus_pv.push_back(bus_id_solver.cast_int());
                has_bus_been_added[bus_id_solver.cast_int()] = true;  // don't add it a second time
            }
        }

        // the reactive output a non-regulating element publishes is its setpoint; a
        // regulating one is filled later (set_q, then LSGrid's write-back or residual share)
        void _compute_res_pq(const Eigen::Ref<const RealVect> & /*Va*/,
                             const Eigen::Ref<const RealVect> & /*Vm*/,
                             const Eigen::Ref<const CplxVect> & /*V*/,
                             const SolverBusIdVect & /*id_grid_to_solver*/,
                             const Eigen::Ref<const RealVect> & /*bus_vn_kv*/,
                             real_type /*sn_mva*/,
                             bool ac) override
        {
            set_osc_pq_res_p();
            if(!ac){
                set_osc_pq_res_q(ac);  // nothing special to do here
                return;
            }
            const int nb_el = nb();
            for(int el_id = 0; el_id < nb_el; ++el_id)
            {
                if(!status_[el_id]){
                    res_q_[el_id] = 0.;  // turned off: no q
                    continue;
                }
                if(voltage_regulator_on_[el_id]) continue;  // filled later by set_q
                res_q_(el_id) = target_q_mvar_(el_id);
            }
        }

        // the one-side index checks, plus the range of `regulated_bus_id_`: a grid bus
        // id the powerflow uses *as an index* without re-checking it
        void _check_valid(int nb_bus,
                          int nb_sub,
                          const SubstationContainer & substations,
                          std::vector<int> & all_pos_topo_vect) const override
        {
            check_valid_osc(nb_bus, nb_sub, substations, all_pos_topo_vect, Leaf::_element_name());
            const int nb_el = nb();
            if(regulated_bus_id_.size() == 0) return;
            for(int el_id = 0; el_id < nb_el; ++el_id)
            {
                // -1 is legal: this element regulates no bus (disconnected / no remote target)
                const int reg = regulated_bus_id_(el_id);
                if(reg == _deactivated_bus_id) continue;
                if((reg < 0) || (reg >= nb_bus))
                {
                    std::ostringstream exc_;
                    exc_ << "LSGrid::check_grid: " << Leaf::_element_name() << " id " << el_id << " regulates bus id "
                         << reg << " which is out of range [0, " << nb_bus << ").";
                    throw std::out_of_range(exc_.str());
                }
            }
        }

        // an element in Sbus whose regulation, when on, pins a bus: its status moves
        // the injections always, and the pv/pq split (and the voltage-control plan
        // built on it) when it regulates
        void _on_deactivate(int el_id, DualAlgoControl & solver_control) override {
            solver_control.tell_recompute_sbus();
            if(voltage_regulator_on_[el_id]) solver_control.tell_pv_changed();
        }
        void _on_reactivate(int el_id, DualAlgoControl & solver_control) override {
            solver_control.tell_recompute_sbus();
            if(voltage_regulator_on_[el_id]) solver_control.tell_pv_changed();
        }
        void _on_change_bus(int el_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control) override {
            // keep a LOCAL regulator local across a bus change: its regulated bus follows
            // its own bus (bus_id_ is still the OLD bus here, reassigned by the caller after).
            // A REMOTE regulator keeps its independent target bus.
            // TODO: a REMOTE regulator's target bus is whatever was resolved at import time
            // (e.g. by `init_from_pypowsybl`). If the *regulated element* itself changes bus
            // here, we have no way to know it (we only store the resolved bus id), so the
            // regulated bus stays frozen and desynchronises from the source grid. Tracking the
            // regulated element id (not just the bus) would let us follow such a move.
            if(regulated_bus_id_(el_id) == bus_id_(el_id).cast_int()){
                regulated_bus_id_(el_id) = new_bus_id.cast_int();
            }
            solver_control.tell_recompute_sbus();
            solver_control.tell_one_el_changed_bus();
            if(voltage_regulator_on_[el_id]) solver_control.tell_pv_changed();
        }

    protected:
        // input data
        std::vector<bool> voltage_regulator_on_;
        RealVect target_vm_pu_;
        // grid bus id whose voltage is regulated (defaults to own bus = local control)
        Eigen::VectorXi regulated_bus_id_;
};

} // namespace ls2g

#endif  // VOLTAGE_SOURCE_CONTAINER_H
