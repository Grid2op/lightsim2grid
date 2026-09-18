// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "StorageContainer.hpp"
#include "BinaryArchive.hpp"

#include <limits>
#include <sstream>
#include <iostream>

namespace ls2g {

void StorageContainer::init(const Eigen::Ref<const RealVect> & storage_p_mw,
                            const Eigen::Ref<const RealVect> & storage_q_mvar,
                            const Eigen::Ref<const Eigen::VectorXi> & storage_bus_id)
{
    // no voltage regulation: the voltage side is filled with "never read" values --
    // an unbounded reactive range, a 1 pu target, the own bus as regulated bus
    const Eigen::Index n = storage_p_mw.size();
    const real_type unbounded = std::numeric_limits<real_type>::infinity();
    init_full(storage_p_mw,
              storage_q_mvar,
              std::vector<bool>(static_cast<std::size_t>(n), false),
              RealVect::Constant(n, 1.),
              RealVect::Constant(n, -unbounded),
              RealVect::Constant(n, unbounded),
              storage_bus_id);
}

void StorageContainer::init_full(const Eigen::Ref<const RealVect> & storage_p_mw,
                                 const Eigen::Ref<const RealVect> & storage_q_mvar,
                                 const std::vector<bool> & voltage_regulator_on,
                                 const Eigen::Ref<const RealVect> & storage_target_vm_pu,
                                 const Eigen::Ref<const RealVect> & storage_min_q,
                                 const Eigen::Ref<const RealVect> & storage_max_q,
                                 const Eigen::Ref<const Eigen::VectorXi> & storage_bus_id)
{
    init_osc_pq(storage_p_mw, storage_q_mvar, storage_bus_id, "storages");

    const int size = nb();
    check_size(voltage_regulator_on, size, "voltage_regulator_on");
    check_size(storage_target_vm_pu, size, "storage_target_vm_pu");
    check_size(storage_min_q, size, "storage_min_q");
    check_size(storage_max_q, size, "storage_max_q");

    target_vm_pu_ = storage_target_vm_pu;
    min_q_ = storage_min_q;
    max_q_ = storage_max_q;
    for(int storage_id = 0; storage_id < size; ++storage_id){
        if (min_q_(storage_id) > max_q_(storage_id))
        {
            std::ostringstream exc_;
            exc_ << "StorageContainer::init: Impossible to initialize storage min_q being above max_q for storage ";
            exc_ << storage_id;
            throw std::runtime_error(exc_.str());
        }
    }
    voltage_regulator_on_ = voltage_regulator_on;
    // local control only: the regulated bus is the unit's own bus (see _check_valid)
    regulated_bus_id_ = storage_bus_id;
    // no unit takes part in the distributed slack until told so (LSGrid::add_storage_slackbus)
    slack_.reset(static_cast<std::size_t>(size));
    reset_results();
}

StorageContainer::StateRes StorageContainer::get_state() const
{
    std::vector<real_type> vm_pu(target_vm_pu_.begin(), target_vm_pu_.end());
    std::vector<real_type> min_q(min_q_.begin(), min_q_.end());
    std::vector<real_type> max_q(max_q_.begin(), max_q_.end());
    std::vector<int> regulated_bus(regulated_bus_id_.begin(), regulated_bus_id_.end());
    StorageContainer::StateRes res(get_osc_pq_state(),  // osc : one side container
                                   voltage_regulator_on_,
                                   vm_pu,
                                   min_q,
                                   max_q,
                                   regulated_bus,
                                   slack_.flags(),
                                   slack_.weights());
    return res;
}

void StorageContainer::set_state(StorageContainer::StateRes & my_state)
{
    set_osc_pq_state(std::get<StateResIdx::OSC_PQ_STATE>(my_state));  // osc : one side container

    std::vector<bool> & voltage_regulator_on = std::get<StateResIdx::VREG_ON>(my_state);
    std::vector<real_type> & vm_pu = std::get<StateResIdx::TARGET_VM_PU>(my_state);
    std::vector<real_type> & min_q = std::get<StateResIdx::MIN_Q>(my_state);
    std::vector<real_type> & max_q = std::get<StateResIdx::MAX_Q>(my_state);
    std::vector<int> & regulated_bus = std::get<StateResIdx::REGULATED_BUS_ID>(my_state);
    std::vector<bool> & slack_bus = std::get<StateResIdx::SLACKBUS>(my_state);
    std::vector<real_type> & slack_weight = std::get<StateResIdx::SLACK_WEIGHT>(my_state);

    const auto size = nb();
    check_size(voltage_regulator_on, size, "voltage_regulator_on");
    check_size(vm_pu, size, "vm_pu");
    check_size(min_q, size, "min_q");
    check_size(max_q, size, "max_q");
    check_size(regulated_bus, size, "regulated_bus");
    check_size(slack_bus, size, "slack_bus");
    check_size(slack_weight, size, "slack_weight");

    voltage_regulator_on_ = voltage_regulator_on;
    target_vm_pu_ = RealVect::Map(vm_pu.data(), vm_pu.size());
    min_q_ = RealVect::Map(min_q.data(), min_q.size());
    max_q_ = RealVect::Map(max_q.data(), max_q.size());
    regulated_bus_id_ = Eigen::VectorXi::Map(regulated_bus.data(), regulated_bus.size());
    slack_.set(slack_bus, slack_weight);
    reset_results();
}

void StorageContainer::save_binary(const std::string & path, bool atomic) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
}

StorageContainer StorageContainer::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<StorageContainer>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

void StorageContainer::_fillSbus(Eigen::Ref<CplxVect> Sbus, const SolverBusIdVect & id_grid_to_solver, bool /*ac*/) const
{
    const int nb_storage = nb();
    for(int storage_id = 0; storage_id < nb_storage; ++storage_id){
        if(!status_[storage_id]) continue;
        const SolverBusId bus_id_solver = _solver_bus(storage_id, bus_id_(storage_id), id_grid_to_solver,
                                                      "StorageContainer::fillSbus");
        // load convention: drawn from the grid
        cplx_type tmp = {-target_p_mw_(storage_id), 0.};
        if(!voltage_regulator_on_[storage_id]){
            // a unit that does not regulate its bus is PQ: its reactive setpoint is drawn too
            tmp -= my_i * target_q_mvar_(storage_id);
        }
        Sbus.coeffRef(bus_id_solver.cast_int()) += tmp;
    }
}

void StorageContainer::_check_valid(int nb_bus,
                                    int nb_sub,
                                    const SubstationContainer & substations,
                                    std::vector<int> & all_pos_topo_vect) const
{
    // one-side index checks (bus / subid / pos_topo_vect) + the regulated bus range
    VoltageSourceContainer<StorageContainer>::_check_valid(nb_bus, nb_sub, substations, all_pos_topo_vect);

    // local regulation only: a storage unit is not enrolled in the VoltageControl plan
    const int nb_storage = nb();
    for(int storage_id = 0; storage_id < nb_storage; ++storage_id){
        if(!is_remote_voltage_controller(storage_id)) continue;
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: storage id " << storage_id << " regulates the voltage of bus "
             << regulated_bus_id_(storage_id) << " while being connected to bus " << bus_id_(storage_id).cast_int()
             << ": remote voltage regulation is not supported for storage units (only the unit's own bus).";
        throw std::runtime_error(exc_.str());
    }

    // a unit flagged as a slack participant carries a usable weight (whether ANY
    // participant is connected is a grid-wide question, see LSGrid::check_grid)
    slack_.check_weights(_element_name());
}

void StorageContainer::_on_deactivate(int storage_id, DualAlgoControl & solver_control) {
    VoltageSourceContainer<StorageContainer>::_on_deactivate(storage_id, solver_control);
    if(slack_.is_slack(storage_id)){ solver_control.tell_slack_participate_changed(); }
}

void StorageContainer::_on_reactivate(int storage_id, DualAlgoControl & solver_control) {
    VoltageSourceContainer<StorageContainer>::_on_reactivate(storage_id, solver_control);
    if(slack_.is_slack(storage_id)){ solver_control.tell_slack_participate_changed(); }
}

void StorageContainer::_on_change_bus(int storage_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control) {
    VoltageSourceContainer<StorageContainer>::_on_change_bus(storage_id, new_bus_id, solver_control);
    if(slack_.is_slack(storage_id)){ solver_control.tell_slack_participate_changed(); }
}

} // namespace ls2g
