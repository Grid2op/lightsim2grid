// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GeneratorContainer.hpp"
#include "BinaryArchive.hpp"

#include <limits>
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
    slack_.reset(static_cast<std::size_t>(generators_p.size()));
    turnedoff_gen_pv_ = true;
    voltage_regulator_on_ = voltage_regulator_on;
    // local control by default: the regulated bus is the generator's own bus
    regulated_bus_id_ = generators_bus_id;
    // no reactive sharing key: the reactive range decides (see has_reactive_key)
    reactive_key_ = RealVect::Constant(size, std::numeric_limits<real_type>::quiet_NaN());
    // the active power limits are optional and belong to the machines this call replaces:
    // a fresh container has none (see set_p_limits). Without this, re-initialising a grid
    // with a different number of generators left limits sized for the previous ones.
    p_min_mw_ = RealVect();
    p_max_mw_ = RealVect();
    // same for the "can be PV" hint: it describes the machines this call replaces
    can_be_pv_ = std::vector<bool>(size, false);
    reset_results();
}

GeneratorContainer::StateRes GeneratorContainer::get_state() const  // osc : one side container
{
     std::vector<real_type> vm_pu(target_vm_pu_.begin(), target_vm_pu_.end());
     std::vector<real_type> min_q(min_q_.begin(), min_q_.end());
     std::vector<real_type> max_q(max_q_.begin(), max_q_.end());
     std::vector<int> regulated_bus(regulated_bus_id_.begin(), regulated_bus_id_.end());
     // optional, and stays empty when unset -- see set_p_limits
     std::vector<real_type> p_min(p_min_mw_.begin(), p_min_mw_.end());
     std::vector<real_type> p_max(p_max_mw_.begin(), p_max_mw_.end());
     std::vector<real_type> reactive_key(reactive_key_.begin(), reactive_key_.end());
     GeneratorContainer::StateRes res(get_osc_pq_state(),  // osc : one side container
                                      turnedoff_gen_pv_,
                                      voltage_regulator_on_,
                                      vm_pu,
                                      min_q,
                                      max_q,
                                      slack_.flags(),
                                      slack_.weights(),
                                      regulated_bus,
                                      p_min,
                                      p_max,
                                      reactive_key,
                                      can_be_pv_);
     return res;
}

void GeneratorContainer::set_state(GeneratorContainer::StateRes & my_state)
{
    set_osc_pq_state(std::get<StateResIdx::OSC_PQ_STATE>(my_state));
    turnedoff_gen_pv_ = std::get<StateResIdx::TURNEDOFF_GEN_PV>(my_state);

    // the generators themelves
    std::vector<bool> & voltage_regulator_on = std::get<StateResIdx::VREG_ON>(my_state);
    std::vector<real_type> & vm_pu = std::get<StateResIdx::TARGET_VM_PU>(my_state);
    std::vector<real_type> & min_q = std::get<StateResIdx::MIN_Q>(my_state);
    std::vector<real_type> & max_q = std::get<StateResIdx::MAX_Q>(my_state);
    std::vector<bool> & slack_bus = std::get<StateResIdx::GEN_SLACKBUS>(my_state);
    std::vector<real_type> & slack_weight = std::get<StateResIdx::GEN_SLACK_WEIGHT>(my_state);
    std::vector<int> & regulated_bus = std::get<StateResIdx::REGULATED_BUS_ID>(my_state);
    std::vector<real_type> & p_min = std::get<StateResIdx::P_MIN_MW>(my_state);
    std::vector<real_type> & p_max = std::get<StateResIdx::P_MAX_MW>(my_state);
    std::vector<real_type> & reactive_key = std::get<StateResIdx::REACTIVE_KEY>(my_state);
    std::vector<bool> & can_be_pv = std::get<StateResIdx::CAN_BE_PV>(my_state);

    // check sizes
    const auto size = nb();
    check_size(voltage_regulator_on, size, "voltage_regulator_on");
    check_size(vm_pu, size, "vm_pu");
    check_size(min_q, size, "min_q");
    check_size(max_q, size, "max_q");
    check_size(slack_bus, size, "slack_bus");
    check_size(slack_weight, size, "slack_weight");
    check_size(regulated_bus, size, "regulated_bus");
    // the active power limits are OPTIONAL: a grid that was never given any carries two
    // empty vectors, and so does a file written before they existed. Only a non-empty one
    // has to match the number of generators.
    if(!p_min.empty() || !p_max.empty()){
        check_size(p_min, size, "p_min_mw");
        check_size(p_max, size, "p_max_mw");
    }
    check_size(reactive_key, size, "reactive_key");
    check_size(can_be_pv, size, "can_be_pv");

    // assign data
    voltage_regulator_on_ = voltage_regulator_on;
    target_vm_pu_ = RealVect::Map(vm_pu.data(), vm_pu.size());
    min_q_ = RealVect::Map(min_q.data(), min_q.size());
    max_q_ = RealVect::Map(max_q.data(), max_q.size());
    slack_.set(slack_bus, slack_weight);
    regulated_bus_id_ = Eigen::VectorXi::Map(regulated_bus.data(), regulated_bus.size());
    p_min_mw_ = p_min.empty() ? RealVect() : RealVect::Map(p_min.data(), p_min.size());
    p_max_mw_ = p_max.empty() ? RealVect() : RealVect::Map(p_max.data(), p_max.size());
    reactive_key_ = RealVect::Map(reactive_key.data(), reactive_key.size());
    can_be_pv_ = can_be_pv;
    reset_results();
}

void GeneratorContainer::_check_valid(int nb_bus,
                                     int nb_sub,
                                     const SubstationContainer & substations,
                                     std::vector<int> & all_pos_topo_vect) const
{
    // one-side index checks (bus / subid / pos_topo_vect) + the regulated bus range
    VoltageSourceContainer<GeneratorContainer>::_check_valid(nb_bus, nb_sub, substations, all_pos_topo_vect);

    // a generator flagged as a slack participant carries a usable weight. Whether ANY
    // participant is connected is a grid-wide question -- storage units take part in
    // the slack too -- answered by LSGrid::check_grid.
    slack_.check_weights(_element_name());
}

void GeneratorContainer::_fillSbus(Eigen::Ref<CplxVect> Sbus, const SolverBusIdVect & id_grid_to_solver, bool /*ac*/) const {
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
#ifndef NDEBUG
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::get_slack_weights_solver: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
#endif
        bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
#ifndef NDEBUG
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::fillSbus: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
#endif
        tmp = {target_p_mw_(gen_id), 0.};
        if(!voltage_regulator_on_[gen_id]){
            // gen is pq if voltage regulaton is off
            tmp += my_i * target_q_mvar_(gen_id);
        }
        Sbus.coeffRef(bus_id_solver.cast_int()) += tmp;
    }
}

void GeneratorContainer::_on_change_p(int gen_id, real_type new_p, DualAlgoControl & solver_control)
{
    if (abs(target_p_mw_(gen_id) - new_p) > _tol_equal_float) {
        solver_control.tell_recompute_sbus();
    }
    if(!turnedoff_gen_pv_){
        // if turned off generators (including these with p==0)
        // are not pv, if we change the active generation, it changes
        // the list of pv buses, so I need to refactorize the solver
        // on the other hand, if all generators are pv then I do not need to refactorize in this case

        if (slack_.is_slack(gen_id)) return;  // slack is not pseudo off
        if (slack_.has_weight(gen_id)) return;  // slack is not pseudo off

        bool pseudo_off_before = abs(target_p_mw_(gen_id)) < _tol_equal_float;
        bool pseudo_off_now = abs(new_p) < _tol_equal_float;
        if((pseudo_off_before && !pseudo_off_now) ||
           (!pseudo_off_before && pseudo_off_now)){
            // (crossing p == 0 also makes this generator start or stop being a voltage
            // controller -- is_remote_voltage_controller gates on is_pseudo_off -- which
            // the pv/pq split, and so the voltage-control plan built around it, already
            // follows from the flag above.)
            solver_control.tell_pv_changed();
           }
    }
}

void GeneratorContainer::_on_deactivate(int el_id, DualAlgoControl & solver_control) {
    VoltageSourceContainer<GeneratorContainer>::_on_deactivate(el_id, solver_control);
    if(!turnedoff_gen_pv_){ solver_control.tell_pv_changed(); }
    if(slack_.is_slack(el_id)){ solver_control.tell_slack_participate_changed(); }
}

void GeneratorContainer::_on_reactivate(int el_id, DualAlgoControl & solver_control) {
    VoltageSourceContainer<GeneratorContainer>::_on_reactivate(el_id, solver_control);
    if(!turnedoff_gen_pv_){ solver_control.tell_pv_changed(); }
    if(slack_.is_slack(el_id)){ solver_control.tell_slack_participate_changed(); }
}

void GeneratorContainer::_on_change_bus(int el_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control) {
    VoltageSourceContainer<GeneratorContainer>::_on_change_bus(el_id, new_bus_id, solver_control);
    if(slack_.is_slack(el_id)) { solver_control.tell_slack_participate_changed(); }
}

void GeneratorContainer::update_slack_weights(
    const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & could_be_slack,
    DualAlgoControl & solver_control)
{
    const int nb_gen = nb();
    // `could_be_slack` comes from python and is indexed by generator id below with
    // an unchecked Eigen operator(): a shorter array would be read out of bounds
    // (release wheels are -O3 -DNDEBUG, so Eigen's own assert is gone).
    if(could_be_slack.rows() != nb_gen){
        std::ostringstream exc_;
        exc_ << "GeneratorContainer::update_slack_weights: 'could_be_slack' has "
             << could_be_slack.rows() << " elements but this grid has " << nb_gen
             << " generators. It is indexed by generator id, so both must match.";
        throw std::runtime_error(exc_.str());
    }
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
