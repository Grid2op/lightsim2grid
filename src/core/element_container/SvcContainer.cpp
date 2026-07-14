// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SvcContainer.hpp"
#include "BinaryArchive.hpp"

#include <iostream>
#include <sstream>

namespace ls2g {

void SvcContainer::init(const std::vector<int> & regulation_mode,
                        const RealVect & target_vm_pu,
                        const RealVect & q_setpoint_mvar,
                        const RealVect & slope_pu,
                        const RealVect & b_min,
                        const RealVect & b_max,
                        const Eigen::VectorXi & regulated_bus_id,
                        const Eigen::VectorXi & bus_id)
{
    const int size = static_cast<int>(bus_id.size());
    // an SVC has NO active power
    const RealVect p_zero = RealVect::Zero(size);
    init_osc_pq(p_zero, q_setpoint_mvar, bus_id, "svc");

    check_size(target_vm_pu, size, "svc_target_vm_pu");
    check_size(q_setpoint_mvar, size, "svc_q_setpoint");
    check_size(slope_pu, size, "svc_slope_pu");
    check_size(b_min, size, "svc_b_min");
    check_size(b_max, size, "svc_b_max");
    check_size(regulation_mode, size, "svc_regulation_mode");
    if(static_cast<int>(regulated_bus_id.size()) != size)
        throw std::runtime_error("SvcContainer::init: regulated_bus_id has a wrong size.");

    regulation_mode_ = IntVect::Map(regulation_mode.data(), regulation_mode.size());
    target_vm_pu_ = target_vm_pu;
    slope_pu_ = slope_pu;
    b_min_ = b_min;
    b_max_ = b_max;
    regulated_bus_id_ = regulated_bus_id;
    reset_results();
}

SvcContainer::StateRes SvcContainer::get_state() const
{
    std::vector<int> mode(regulation_mode_.begin(), regulation_mode_.end());
    std::vector<real_type> vm_pu(target_vm_pu_.begin(), target_vm_pu_.end());
    std::vector<real_type> slope(slope_pu_.begin(), slope_pu_.end());
    std::vector<real_type> bmin(b_min_.begin(), b_min_.end());
    std::vector<real_type> bmax(b_max_.begin(), b_max_.end());
    std::vector<int> regulated_bus(regulated_bus_id_.begin(), regulated_bus_id_.end());
    SvcContainer::StateRes res(get_osc_pq_state(), mode, vm_pu, slope, bmin, bmax, regulated_bus);
    return res;
}

void SvcContainer::set_state(SvcContainer::StateRes & my_state)
{
    set_osc_pq_state(std::get<StateResIdx::OSC_PQ_STATE>(my_state));
    std::vector<int> & mode = std::get<StateResIdx::REGULATION_MODE>(my_state);
    std::vector<real_type> & vm_pu = std::get<StateResIdx::TARGET_VM_PU>(my_state);
    std::vector<real_type> & slope = std::get<StateResIdx::SLOPE_PU>(my_state);
    std::vector<real_type> & bmin = std::get<StateResIdx::B_MIN>(my_state);
    std::vector<real_type> & bmax = std::get<StateResIdx::B_MAX>(my_state);
    std::vector<int> & regulated_bus = std::get<StateResIdx::REGULATED_BUS_ID>(my_state);

    const auto size = nb();
    check_size(mode, size, "regulation_mode");
    check_size(vm_pu, size, "target_vm_pu");
    check_size(slope, size, "slope_pu");
    check_size(bmin, size, "b_min");
    check_size(bmax, size, "b_max");
    check_size(regulated_bus, size, "regulated_bus");

    regulation_mode_ = IntVect::Map(mode.data(), mode.size());
    target_vm_pu_ = RealVect::Map(vm_pu.data(), vm_pu.size());
    slope_pu_ = RealVect::Map(slope.data(), slope.size());
    b_min_ = RealVect::Map(bmin.data(), bmin.size());
    b_max_ = RealVect::Map(bmax.data(), bmax.size());
    regulated_bus_id_ = Eigen::VectorXi::Map(regulated_bus.data(), regulated_bus.size());
    reset_results();
}

void SvcContainer::fillSbus(CplxVect & Sbus, const SolverBusIdVect & id_grid_to_solver, bool /*ac*/) const
{
    const int nb_svc = nb();
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(!status_[svc_id]) continue;
        // only the REACTIVE_POWER mode injects into Sbus (P = 0). VOLTAGE mode is
        // solved by the VoltageControl extension; OFF stamps nothing.
        if(regulation_mode_(svc_id) != RegulationMode::REACTIVE_POWER) continue;

        const GlobalBusId bus_id_me = bus_id_(svc_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "SvcContainer::fillSbus: Svc with id " << svc_id
                 << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        const SolverBusId bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "SvcContainer::fillSbus: Svc with id " << svc_id
                 << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        // P = 0, Q = setpoint (generator injection convention)
        Sbus.coeffRef(bus_id_solver.cast_int()) += my_i * target_q_mvar_(svc_id);
    }
}

void SvcContainer::get_vm_for_dc(RealVect & Vm)
{
    const int nb_svc = nb();
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(!status_[svc_id]) continue;
        if(regulation_mode_(svc_id) != RegulationMode::VOLTAGE) continue;
        const int target_bus = regulated_bus_id_(svc_id);
        const real_type tmp = target_vm_pu_(svc_id);
        if(abs(tmp) > _tol_equal_float) Vm(target_bus) = tmp;
    }
}

void SvcContainer::set_vm(CplxVect & V, const SolverBusIdVect & id_grid_to_solver) const
{
    const int nb_svc = nb();
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(!status_[svc_id]) continue;
        if(regulation_mode_(svc_id) != RegulationMode::VOLTAGE) continue;
        // a sloped SVC does NOT hold its regulated bus exactly at the setpoint
        // (Vm = v_set - s.Q): forcing it here would corrupt an already-solved V
        // (e.g. in check_solution). The flat init is good enough for the NR.
        if(has_slope(svc_id)) continue;

        const int target_grid_bus = regulated_bus_id_(svc_id);
        if(target_grid_bus == _deactivated_bus_id) continue;
        const SolverBusId bus_id_solver = id_grid_to_solver[target_grid_bus];
        if(bus_id_solver.cast_int() == _deactivated_bus_id) continue;
        real_type tmp = std::abs(V(bus_id_solver.cast_int()));
        if(abs(tmp) < _tol_equal_float){
            V(bus_id_solver.cast_int()) = 1.0;
            tmp = 1.0;
        }
        tmp = target_vm_pu_(svc_id) / tmp;
        V(bus_id_solver.cast_int()) *= tmp;
    }
}

void SvcContainer::_compute_results(
    const Eigen::Ref<const RealVect> & /*Va*/,
    const Eigen::Ref<const RealVect> & /*Vm*/,
    const Eigen::Ref<const CplxVect> & /*V*/,
    const SolverBusIdVect & /*id_grid_to_solver*/,
    const Eigen::Ref<const RealVect> & /*bus_vn_kv*/,
    real_type /*sn_mva*/,
    bool ac)
{
    const int nb_svc = nb();
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        res_p_(svc_id) = 0.;  // an SVC has no active power
        if(!ac){
            res_q_(svc_id) = 0.;  // no reactive result in DC
            continue;
        }
        if(!status_[svc_id]){
            res_q_(svc_id) = 0.;
            continue;
        }
        if(regulation_mode_(svc_id) == RegulationMode::REACTIVE_POWER){
            res_q_(svc_id) = target_q_mvar_(svc_id);
        } else {
            // VOLTAGE mode: set by the VoltageControl write-back (LSGrid::compute_results);
            // OFF: nothing. Initialise to 0 here (write-back overrides VOLTAGE).
            res_q_(svc_id) = 0.;
        }
    }
}

bool SvcContainer::_deactivate(int svc_id, DualAlgoControl & solver_control)
{
    if(status_[svc_id]){
        solver_control.ac_algo_controler().tell_recompute_sbus();
        solver_control.ac_algo_controler().tell_one_el_changed_bus();
        solver_control.dc_algo_controler().tell_recompute_sbus();
        solver_control.dc_algo_controler().tell_one_el_changed_bus();
        if(regulation_mode_(svc_id) == RegulationMode::VOLTAGE){
            solver_control.ac_algo_controler().tell_pv_changed();
            solver_control.dc_algo_controler().tell_pv_changed();
        }
        return true;
    }
    return false;
}

bool SvcContainer::_reactivate(int svc_id, DualAlgoControl & solver_control)
{
    if(!status_[svc_id]){
        solver_control.ac_algo_controler().tell_recompute_sbus();
        solver_control.ac_algo_controler().tell_one_el_changed_bus();
        solver_control.dc_algo_controler().tell_recompute_sbus();
        solver_control.dc_algo_controler().tell_one_el_changed_bus();
        if(regulation_mode_(svc_id) == RegulationMode::VOLTAGE){
            solver_control.ac_algo_controler().tell_pv_changed();
            solver_control.dc_algo_controler().tell_pv_changed();
        }
        return true;
    }
    return false;
}

bool SvcContainer::_change_bus(int svc_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control, int /*nb_bus*/)
{
    // el_id is validated (and the proper IndexError raised) by `_generic_change_bus`,
    // which the caller runs *after* this function. Bail out here on an out-of-range
    // id so the `regulated_bus_id_` write below never touches memory out of bounds.
    if(svc_id < 0 || svc_id >= nb()) return false;
    if(bus_id_(svc_id) == new_bus_id) return false;
    // a LOCAL voltage controller's regulated bus follows its own bus
    // TODO: a REMOTE controller's regulated bus is whatever was resolved at import time
    // (e.g. by `init_from_pypowsybl`); if the regulated *element* changes bus we cannot
    // tell (only the resolved bus id is stored), so it stays frozen and desynchronises
    // from the source grid. Tracking the regulated element id would let us follow it.
    if(regulated_bus_id_(svc_id) == bus_id_(svc_id).cast_int()) regulated_bus_id_(svc_id) = new_bus_id.cast_int();
    solver_control.ac_algo_controler().tell_recompute_sbus();
    solver_control.ac_algo_controler().tell_one_el_changed_bus();
    solver_control.dc_algo_controler().tell_recompute_sbus();
    solver_control.dc_algo_controler().tell_one_el_changed_bus();
    if(regulation_mode_(svc_id) == RegulationMode::VOLTAGE){
        solver_control.ac_algo_controler().tell_pv_changed();
        solver_control.dc_algo_controler().tell_pv_changed();
    }
    return true;
}

void SvcContainer::save_binary(const std::string & path, bool atomic) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
}

SvcContainer SvcContainer::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<SvcContainer>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

} // namespace ls2g
