// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GeneratorContainer.hpp"
#include "BinaryArchive.hpp"

#include <iostream>
#include <sstream>
#include <cmath>  // for std::isfinite (check_valid)

namespace ls2g {

void GeneratorContainer::init(const Eigen::Ref<const RealVect> & generators_p,
                              const Eigen::Ref<const RealVect> & generators_v,
                              const Eigen::Ref<const RealVect> & generators_min_q,
                              const Eigen::Ref<const RealVect> & generators_max_q,
                              const Eigen::Ref<const Eigen::VectorXi> & generators_bus_id)
{
    const auto generators_q = RealVect::Zero(generators_p.size());
    const auto voltage_regulator_on = std::vector<bool>(generators_p.size(), true);
    init_full(
        generators_p,
        generators_v,
        generators_q,
        voltage_regulator_on,
        generators_min_q,
        generators_max_q,
        generators_bus_id);
}

void GeneratorContainer::init_full(const Eigen::Ref<const RealVect> & generators_p,
                                   const Eigen::Ref<const RealVect> & generators_v,
                                   const Eigen::Ref<const RealVect> & generators_q,
                                   const std::vector<bool> & voltage_regulator_on,
                                   const Eigen::Ref<const RealVect> & generators_min_q,
                                   const Eigen::Ref<const RealVect> & generators_max_q,
                                   const Eigen::Ref<const Eigen::VectorXi> & generators_bus_id
                                   )
{
    init_osc_pq(generators_p, generators_q, generators_bus_id, "generators");

    // check the sizes
    int size = nb();
    check_size(generators_v, size, "generators_v");
    check_size(generators_q, size, "generators_q");
    check_size(generators_min_q, size, "generators_min_q");
    check_size(generators_max_q, size, "generators_max_q");
    check_size(voltage_regulator_on, size, "voltage_regulator_on");

    // fill the data
    target_vm_pu_ = generators_v;
    min_q_ = generators_min_q;
    max_q_ = generators_max_q;
    for(int gen_id = 0; gen_id < size; ++gen_id){
        if (min_q_(gen_id) > max_q_(gen_id))
        {
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::init: Impossible to initialize generator min_q being above max_q for generator ";
            exc_ << gen_id;
            throw std::runtime_error(exc_.str());
        }
    }
    gen_slackbus_ = std::vector<bool>(generators_p.size(), false);
    gen_slack_weight_ = std::vector<real_type>(generators_p.size(), 0.);
    turnedoff_gen_pv_ = true;
    voltage_regulator_on_ = voltage_regulator_on;
    // local control by default: the regulated bus is the generator's own bus
    regulated_bus_id_ = generators_bus_id;
    reset_results();
}

GeneratorContainer::StateRes GeneratorContainer::get_state() const  // osc : one side container
{
     std::vector<real_type> vm_pu(target_vm_pu_.begin(), target_vm_pu_.end());
     std::vector<real_type> min_q(min_q_.begin(), min_q_.end());
     std::vector<real_type> max_q(max_q_.begin(), max_q_.end());
     std::vector<int> regulated_bus(regulated_bus_id_.begin(), regulated_bus_id_.end());
     GeneratorContainer::StateRes res(get_osc_pq_state(),  // osc : one side container
                                      turnedoff_gen_pv_,
                                      voltage_regulator_on_,
                                      vm_pu,
                                      min_q,
                                      max_q,
                                      gen_slackbus_,
                                      gen_slack_weight_,
                                      regulated_bus);
     return res;
}

void GeneratorContainer::set_state(GeneratorContainer::StateRes & my_state)
{
    set_osc_pq_state(std::get<0>(my_state));
    turnedoff_gen_pv_ = std::get<1>(my_state);

    // the generators themelves
    std::vector<bool> & voltage_regulator_on = std::get<2>(my_state);
    std::vector<real_type> & vm_pu = std::get<3>(my_state);
    std::vector<real_type> & min_q = std::get<4>(my_state);
    std::vector<real_type> & max_q = std::get<5>(my_state);
    std::vector<bool> & slack_bus = std::get<6>(my_state);
    std::vector<real_type> & slack_weight = std::get<7>(my_state);
    std::vector<int> & regulated_bus = std::get<8>(my_state);

    // check sizes
    const auto size = nb();
    check_size(voltage_regulator_on, size, "voltage_regulator_on");
    check_size(vm_pu, size, "vm_pu");
    check_size(min_q, size, "min_q");
    check_size(max_q, size, "max_q");
    check_size(slack_bus, size, "slack_bus");
    check_size(slack_weight, size, "slack_weight");
    check_size(regulated_bus, size, "regulated_bus");

    // assign data
    voltage_regulator_on_ = voltage_regulator_on;
    target_vm_pu_ = RealVect::Map(vm_pu.data(), vm_pu.size());
    min_q_ = RealVect::Map(min_q.data(), min_q.size());
    max_q_ = RealVect::Map(max_q.data(), max_q.size());
    gen_slackbus_ = slack_bus;
    gen_slack_weight_ = slack_weight;
    regulated_bus_id_ = Eigen::VectorXi::Map(regulated_bus.data(), regulated_bus.size());
    reset_results();
}

void GeneratorContainer::check_valid(int nb_bus,
                                     int nb_sub,
                                     const SubstationContainer & substations,
                                     std::vector<int> & all_pos_topo_vect) const
{
    // one-side (bus / subid / pos_topo_vect) + (p, q) finiteness checks
    check_valid_osc_pq(nb_bus, nb_sub, substations, all_pos_topo_vect, "generator");

    // generator-specific electrical inputs must be finite
    check_all_finite(target_vm_pu_, "generator target_vm_pu");
    check_all_finite(min_q_, "generator min_q");
    check_all_finite(max_q_, "generator max_q");

    // slack coherence + remote-regulated bus id range
    const int nb_gen = nb();
    const bool has_slack_info = !gen_slackbus_.empty();
    const bool has_reg_info = regulated_bus_id_.size() > 0;
    bool any_slack = false;
    bool any_connected_slack = false;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(has_reg_info)
        {
            // regulated bus may be -1 (no remote regulation / disconnected)
            const int reg = regulated_bus_id_(gen_id);
            if((reg != _deactivated_bus_id) && ((reg < 0) || (reg >= nb_bus)))
            {
                std::ostringstream exc_;
                exc_ << "LSGrid::check_grid: generator id " << gen_id << " regulates bus id "
                     << reg << " which is out of range [0, " << nb_bus << ").";
                throw std::out_of_range(exc_.str());
            }
        }
        if(has_slack_info && gen_slackbus_[gen_id])
        {
            any_slack = true;
            const real_type w = gen_slack_weight_[gen_id];
            if((!std::isfinite(w)) || (w <= _tol_equal_float))
            {
                std::ostringstream exc_;
                exc_ << "LSGrid::check_grid: generator id " << gen_id
                     << " is flagged as a slack but has a non-positive or non-finite slack weight ("
                     << w << ").";
                throw std::runtime_error(exc_.str());
            }
            if(status_[gen_id]) any_connected_slack = true;
        }
    }
    // if a slack is declared at all, at least one slack generator must be connected
    // (the powerflow cannot solve otherwise). We do NOT require a slack to exist:
    // that stays the solver's responsibility, exactly as before.
    if(any_slack && !any_connected_slack)
    {
        throw std::runtime_error("LSGrid::check_grid: at least one generator is flagged as a "
                                 "slack, but none of the slack generators is connected.");
    }
}

RealVect GeneratorContainer::get_slack_weights_solver(
    size_t nb_bus_solver,
    const SolverBusIdVect & id_grid_to_solver){
    const int nb_gen = nb();
    GlobalBusId bus_id_me;
    SolverBusId bus_id_solver;
    RealVect res = RealVect::Zero(nb_bus_solver);
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the load is disconnected
        if(!status_[gen_id]) continue;
        if(!gen_slackbus_[gen_id]) continue;
        if(abs(gen_slack_weight_[gen_id]) < _tol_equal_float) continue;

        bus_id_me = bus_id_(gen_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::get_slack_weights_solver: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::get_slack_weights_solver: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        if(gen_slackbus_[gen_id]) res.coeffRef(bus_id_solver.cast_int()) += gen_slack_weight_[gen_id];
    }
    bus_slack_weight_ = res;
    real_type sum_res = res.sum();
    res /= sum_res;
    return res;
}

void GeneratorContainer::fillSbus(Eigen::Ref<CplxVect> Sbus, const SolverBusIdVect & id_grid_to_solver, bool /*ac*/) const {
    const int nb_gen = nb();
    GlobalBusId bus_id_me;
    SolverBusId bus_id_solver;
    cplx_type tmp;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the gen is disconnected
        if(!status_[gen_id]) continue;

        // a pv gen that is "pseudo off" (if the flag is set) is turned off, so disconnected
        if ((!turnedoff_gen_pv_) && is_pseudo_off(gen_id) && voltage_regulator_on_[gen_id]) continue;  

        bus_id_me = bus_id_(gen_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::get_slack_weights_solver: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::fillSbus: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        tmp = {target_p_mw_(gen_id), 0.};
        if(!voltage_regulator_on_[gen_id]){
            // gen is pq if voltage regulaton is off
            tmp += my_i * target_q_mvar_(gen_id);
        }
        Sbus.coeffRef(bus_id_solver.cast_int()) += tmp;
    }
}

void GeneratorContainer::fillpv(std::vector<int> & bus_pv,
                                std::vector<bool> & has_bus_been_added,
                                const SolverBusIdVect & slack_bus_id_solver,
                                const SolverBusIdVect & id_grid_to_solver) const
{
    const int nb_gen = nb();
    GlobalBusId bus_id_me;
    SolverBusId bus_id_solver;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;

        // gen is purposedly not pv
        if (!voltage_regulator_on_[gen_id]) continue;

        // a remote-regulating gen does NOT go through the PV path: it is a
        // controller of a VoltageControl group instead (its own bus stays PQ)
        if (regulates_remote(gen_id)) continue;

        // in this case turned off generators are not pv
        // except the slack that can have a target of 0MW but is still "on"
        // no matter what
        bool gen_pseudo_off = is_pseudo_off(gen_id);
        // if (gen_slack_weight_[gen_id] != 0.) gen_pseudo_off = false;  // useless: slack is not PV anyway
        if ((!turnedoff_gen_pv_) && gen_pseudo_off) continue;  
        bus_id_me = bus_id_(gen_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::get_slack_weights_solver: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::fillpv: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }

        if(is_in_vect(bus_id_solver.cast_int(), slack_bus_id_solver.to_int_vector())) continue;  // slack bus is not PV
        if(has_bus_been_added[bus_id_solver.cast_int()]) continue; // i already added this bus
        bus_pv.push_back(bus_id_solver.cast_int());
        has_bus_been_added[bus_id_solver.cast_int()] = true;  // don't add it a second time
    }
}

void GeneratorContainer::get_vm_for_dc(Eigen::Ref<RealVect> Vm){
    const int nb_gen = nb();
    GlobalBusId bus_id_me;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;

        if (!voltage_regulator_on_[gen_id]) continue;  // gen is purposedly not pv
        if ((!turnedoff_gen_pv_) && is_pseudo_off(gen_id)) continue;  // in this case turned off generators are not pv

        // a remote-regulating gen sets the magnitude of the REGULATED bus
        const int target_bus = regulated_bus_id_(gen_id);
        real_type tmp = target_vm_pu_(gen_id);
        if(abs(tmp) > _tol_equal_float) Vm(target_bus) = tmp;
    }
}

void GeneratorContainer::_change_p(int gen_id, real_type new_p, bool /*my_status*/, DualAlgoControl & solver_control)
{
    if (abs(target_p_mw_(gen_id) - new_p) > _tol_equal_float) {
        solver_control.ac_algo_controler().tell_recompute_sbus(); solver_control.dc_algo_controler().tell_recompute_sbus();
    }
    if(!turnedoff_gen_pv_){
        // if turned off generators (including these with p==0)
        // are not pv, if we change the active generation, it changes
        // the list of pv buses, so I need to refactorize the solver
        // on the other hand, if all generators are pv then I do not need to refactorize in this case

        if (gen_slackbus_[gen_id]) return;  // slack is not pseudo off
        if ((abs(gen_slack_weight_[gen_id]) >= _tol_equal_float)) return;  // slack is not pseudo off

        bool pseudo_off_before = abs(target_p_mw_(gen_id)) < _tol_equal_float;
        bool pseudo_off_now = abs(new_p) < _tol_equal_float;
        if((pseudo_off_before && !pseudo_off_now) || 
           (!pseudo_off_before && pseudo_off_now)){
            solver_control.ac_algo_controler().tell_pv_changed(); solver_control.dc_algo_controler().tell_pv_changed();
           }
    }
}

bool GeneratorContainer::_deactivate(int el_id, DualAlgoControl & solver_control) {
    if(!status_[el_id]) return false;  // nothing to do if it was already deactivated
    solver_control.ac_algo_controler().tell_recompute_sbus(); solver_control.dc_algo_controler().tell_recompute_sbus();
    if(voltage_regulator_on_[el_id]){ solver_control.ac_algo_controler().tell_pv_changed(); solver_control.dc_algo_controler().tell_pv_changed(); }
    if(!turnedoff_gen_pv_){ solver_control.ac_algo_controler().tell_pv_changed(); solver_control.dc_algo_controler().tell_pv_changed(); }
    if(gen_slackbus_[el_id]){ solver_control.ac_algo_controler().tell_slack_participate_changed(); solver_control.dc_algo_controler().tell_slack_participate_changed(); }
    return true;
};

bool GeneratorContainer::_reactivate(int el_id, DualAlgoControl & solver_control) {
    if(status_[el_id]) return false;  // nothing to do if gen already connected
    solver_control.ac_algo_controler().tell_recompute_sbus(); solver_control.dc_algo_controler().tell_recompute_sbus();
    if(voltage_regulator_on_[el_id]){ solver_control.ac_algo_controler().tell_pv_changed(); solver_control.dc_algo_controler().tell_pv_changed(); }
    if(!turnedoff_gen_pv_){ solver_control.ac_algo_controler().tell_pv_changed(); solver_control.dc_algo_controler().tell_pv_changed(); }
    if(gen_slackbus_[el_id]){ solver_control.ac_algo_controler().tell_slack_participate_changed(); solver_control.dc_algo_controler().tell_slack_participate_changed(); }
    return true;
};

void GeneratorContainer::change_v(int gen_id, real_type new_v_pu, DualAlgoControl & solver_control)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        // TODO DEBUG MODE only this in debug mode
        std::ostringstream exc_;
        exc_ << "GeneratorContainer::change_p: Impossible to change the voltage setpoint of a disconnected generator (check gen. id ";
        exc_ << gen_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    change_v_nothrow(gen_id, new_v_pu, solver_control);
}

void GeneratorContainer::change_v_nothrow(int gen_id, real_type new_v_pu, DualAlgoControl & solver_control)
{
    [[maybe_unused]] bool my_status = status_.at(gen_id); // and this check that gen_id is not out of bound [[maybe_unused]] 
    if (abs(target_vm_pu_(gen_id) - new_v_pu) > _tol_equal_float)
    {
        solver_control.ac_algo_controler().tell_v_changed(); solver_control.dc_algo_controler().tell_v_changed();
        target_vm_pu_(gen_id) = new_v_pu;
    }
}

bool GeneratorContainer::_change_bus(int el_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control, int /*nb_bus*/) {
    // el_id is validated (and the proper IndexError raised) by `_generic_change_bus`,
    // which the caller runs *after* this function. Bail out here on an out-of-range
    // id so the `regulated_bus_id_` write below never touches memory out of bounds.
    if(el_id < 0 || el_id >= nb()) return false;
    if(bus_id_(el_id) == new_bus_id) return false;  // nothing to do if the bus did not changed
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
    solver_control.ac_algo_controler().tell_recompute_sbus(); solver_control.dc_algo_controler().tell_recompute_sbus();
    solver_control.ac_algo_controler().tell_one_el_changed_bus(); solver_control.dc_algo_controler().tell_one_el_changed_bus();
    if(voltage_regulator_on_[el_id]) { solver_control.ac_algo_controler().tell_pv_changed(); solver_control.dc_algo_controler().tell_pv_changed(); }
    if(gen_slackbus_[el_id]) { solver_control.ac_algo_controler().tell_slack_participate_changed(); solver_control.dc_algo_controler().tell_slack_participate_changed(); }
    return true;
// bool GeneratorContainer::_change_bus(int el_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control, int nb_bus) {
//     if(bus_id_(el_id) == new_bus_id) return false;  // nothing to do if the bus did not changed
//     solver_control.ac_algo_controler().tell_recompute_sbus(); solver_control.dc_algo_controler().tell_recompute_sbus();
//     solver_control.ac_algo_controler().tell_one_el_changed_bus(); solver_control.dc_algo_controler().tell_one_el_changed_bus();
//     if(voltage_regulator_on_[el_id]){ solver_control.ac_algo_controler().tell_pv_changed(); solver_control.dc_algo_controler().tell_pv_changed(); }
//     if(gen_slackbus_[el_id]){ solver_control.ac_algo_controler().tell_slack_participate_changed(); solver_control.dc_algo_controler().tell_slack_participate_changed(); }
    // return true;
};

void GeneratorContainer::set_vm(Eigen::Ref<CplxVect> V, const SolverBusIdVect & id_grid_to_solver) const
{
    const int nb_gen = nb();
    SolverBusId bus_id_solver;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;
        
        if (!voltage_regulator_on_[gen_id]) continue;  // gen is purposedly not pv

        
        bool pseudo_off = is_pseudo_off(gen_id);
        if ((!turnedoff_gen_pv_) && pseudo_off) continue;  // in this case turned off generators are not pv

        // a remote-regulating gen sets the magnitude of the REGULATED bus (init quality)
        const int target_grid_bus = regulated_bus_id_(gen_id);
        if(target_grid_bus == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::set_vm: Generator with id ";
            exc_ << gen_id;
            exc_ << " regulates a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver = id_grid_to_solver[target_grid_bus];
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::set_vm: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        // scale the input V such that abs(V) = Vm for this generator
        real_type tmp = std::abs(V(bus_id_solver.cast_int()));
        if(abs(tmp) < _tol_equal_float)
        {
            // if it was 0. i force it to 1. (otherwise the rest of the computation would make it O. still)
            V(bus_id_solver.cast_int()) = 1.0;
            tmp = 1.0;
        }
        tmp = 1.0 / tmp;
        tmp *= target_vm_pu_(gen_id);
        V(bus_id_solver.cast_int()) *= tmp;
    }
}

GlobalBusIdVect GeneratorContainer::get_slack_bus_id() const{
    std::vector<int> tmp;
    tmp.reserve(gen_slackbus_.size());
    GlobalBusIdVect res;
    const int nb_gen = nb();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(gen_slackbus_[gen_id]){
            const GlobalBusId my_bus = bus_id_(gen_id);
            // do not add twice the same "slack bus"
            if(!is_in_vect(my_bus.cast_int(), tmp)) tmp.push_back(my_bus.cast_int());
        }
    }
    if(tmp.empty()) throw std::runtime_error("GeneratorContainer::get_slack_bus_id: no generator are tagged slack bus for this grid.");
    res = GlobalBusIdVect(tmp);
    return res;
}

void GeneratorContainer::set_p_slack(const Eigen::Ref<const RealVect>& node_mismatch,
                                     const SolverBusIdVect & id_grid_to_solver)
{
    if(bus_slack_weight_.size() == 0){
        // TODO DEBUG MODE: perform this check only in debug mode
        throw std::runtime_error("Generator::set_p_slack: Impossible to set the active value of generators for the slack bus: no known slack (you should haved called Generator::get_slack_weights_solver first)");
    }
    const auto nb_gen = nb();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(!status_[gen_id]) continue;  // nothing to do if gen is disconnected
        if(!gen_slackbus_[gen_id]) continue;  // nothing to do if it's not a slack
        if(abs(gen_slack_weight_[gen_id]) < _tol_equal_float) continue; // nothing to do if no weights are associated to it
        const GlobalBusId bus_id_me = bus_id_(gen_id);
        const SolverBusId bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
        // TODO DEBUG MODE: check bus_id_solver >= 0
        // TODO DEBUG MODE: check bus_slack_weight_[bus_id_solver] > 0
        const real_type total_contrib_slack = bus_slack_weight_(bus_id_solver.cast_int());
        const real_type my_contrib_slack = gen_slack_weight_[gen_id];
        res_p_(gen_id) += node_mismatch(bus_id_solver.cast_int()) * my_contrib_slack / total_contrib_slack;
    }
}

void GeneratorContainer::init_q_vector(int /*nb_bus*/,
                                       Eigen::Ref<Eigen::VectorXi> total_gen_per_bus,
                                       Eigen::Ref<RealVect> total_q_min_per_bus,
                                       Eigen::Ref<RealVect> total_q_max_per_bus) const
{
    const int nb_gen = nb();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(!status_[gen_id]) continue;

        if (!voltage_regulator_on_[gen_id]) continue;  // gen is purposedly not pv
        if ((!turnedoff_gen_pv_) && is_pseudo_off(gen_id)) continue;  // in this case "turned off" generators are not pv
        // a remote-regulating gen's Q comes from the VoltageControl extension, not
        // from the per-bus reactive redistribution: exclude it from the bus totals
        if (regulates_remote(gen_id)) continue;

        const GlobalBusId bus_id = bus_id_(gen_id);
        total_q_min_per_bus(bus_id.cast_int()) += min_q_(gen_id);
        total_q_max_per_bus(bus_id.cast_int()) += max_q_(gen_id);
        total_gen_per_bus(bus_id.cast_int()) += 1;
    }
}

void GeneratorContainer::set_q(
    const Eigen::Ref<const RealVect> & reactive_mismatch,
    const SolverBusIdVect & id_grid_to_solver,
    bool ac,
    const Eigen::Ref<const Eigen::VectorXi> & total_gen_per_bus,
    const Eigen::Ref<const RealVect> & total_q_min_per_bus,
    const Eigen::Ref<const RealVect> & total_q_max_per_bus)
{
    const int nb_gen = nb();
    if(!ac){
        // do not consider Q values in dc mode
        for(int gen_id = 0; gen_id < nb_gen; ++gen_id) res_q_(gen_id) = 0.;
        return;
    }
    
    real_type eps_q = 1e-8;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(!status_[gen_id]){
            // set at 0 for disconnected generators
            res_q_(gen_id) = 0.;
            continue;  
        }
        real_type real_q = 0.;
        if (!voltage_regulator_on_[gen_id]){
            // gen is purposedly not pv, so output MVAr = input MVAr (just like a load)
            res_q_(gen_id) = target_q_mvar_(gen_id);
            continue;
        } 
        if ((!turnedoff_gen_pv_) && is_pseudo_off(gen_id)) {
            // in this case turned off generators are not pv
            // it's as if the generator were turned off
            res_q_(gen_id) = 0.;
            continue;
        }
        if (regulates_remote(gen_id)) {
            // a remote-regulating gen's reactive output is set by the VoltageControl
            // write-back (LSGrid::compute_results), not by the per-bus redistribution
            continue;
        }

        const GlobalBusId bus_id = bus_id_(gen_id);
        const SolverBusId bus_solver = id_grid_to_solver[bus_id.cast_int()];
        // TODO DEBUG MODE: check that the bus is correct!
        real_type q_to_absorb = reactive_mismatch[bus_solver.cast_int()];
        real_type max_q_me = max_q_(gen_id);
        real_type min_q_me = min_q_(gen_id);
        real_type max_q_bus = total_q_max_per_bus(bus_id.cast_int());
        real_type min_q_bus = total_q_min_per_bus(bus_id.cast_int());
        real_type nb_gen_with_me = static_cast<real_type>(total_gen_per_bus(bus_id.cast_int()));
        if(nb_gen_with_me == 1.){
            real_q = q_to_absorb;
        }else{
            real_type ratio = (max_q_me - min_q_me + eps_q) / (max_q_bus - min_q_bus + nb_gen_with_me * eps_q) ;
            real_q = q_to_absorb * ratio ;
        }
        res_q_(gen_id) = real_q;
    }
}

void GeneratorContainer::update_slack_weights(
    const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & could_be_slack,
    DualAlgoControl & solver_control)
{
    const int nb_gen = nb();
    std::vector<int> gen_slack_id;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(could_be_slack(gen_id)) gen_slack_id.push_back(gen_id);
    }
    Eigen::Ref<const IntVect> gen_slack_id_ref = IntVect::Map(gen_slack_id.data(), gen_slack_id.size());
    update_slack_weights_by_id(
        gen_slack_id_ref,
        solver_control);
}

void GeneratorContainer::update_slack_weights_by_id(
    const Eigen::Ref<const IntVect> & gen_slack_id,
    DualAlgoControl & solver_control)
{
    // TODO speed: the solver_control will always tell that the slacks changed
    // even if it's not the case.
    // Because the 
    int nb_gen = nb();
    std::vector<bool> maybe_slack_bus(nb_gen, false);

    // validate every caller-supplied id before it is used to index status_ /
    // maybe_slack_bus / target_p_mw_ below (raw operator[] / Eigen operator() are
    // unchecked, and a negative id would wrap to a huge size_t -> OOB write).
    for(int gen_id : gen_slack_id) _check_in_range(gen_id, status_, "update_slack_weights_by_id");

    // find which generators can be slack
    real_type total_target_p = 0.;
    for(int gen_id : gen_slack_id)
    {
        if(status_[gen_id])
        {
            maybe_slack_bus[gen_id] = true;
            total_target_p += abs(target_p_mw_(gen_id));
        }
    }

    // assign the slack to the generators
    if(abs(total_target_p) < _tol_equal_float){
        // all gen to the slacks produces 0.
        // slacks weights are equal for all generators
        real_type slack_weight = 1. / static_cast<real_type>(gen_slack_id.size());
        for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
            if(maybe_slack_bus[gen_id]) 
                add_slackbus(gen_id, slack_weight, solver_control);
            else remove_slackbus(gen_id, solver_control);
        }
    }else{
        // slack weights prop to abs(target_p)
        for(int gen_id : gen_slack_id)
        {
            if(maybe_slack_bus[gen_id] && (abs(target_p_mw_[gen_id]) > _tol_equal_float)) 
                add_slackbus(gen_id, abs(target_p_mw_[gen_id]), solver_control);
            else remove_slackbus(gen_id, solver_control);
        }
    }
}

void GeneratorContainer::save_binary(const std::string & path, bool atomic) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
}

GeneratorContainer GeneratorContainer::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<GeneratorContainer>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

} // namespace ls2g
