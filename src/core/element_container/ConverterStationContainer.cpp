// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ConverterStationContainer.hpp"

#include <cmath>
#include <sstream>

namespace ls2g {

void ConverterStationContainer::init(const std::vector<int> & type,
                                     const Eigen::Ref<const RealVect> & loss_factor,
                                     const std::vector<bool> & voltage_regulator_on,
                                     const Eigen::Ref<const RealVect> & target_vm_pu,
                                     const Eigen::Ref<const RealVect> & q_setpoint_mvar,
                                     const Eigen::Ref<const RealVect> & min_q,
                                     const Eigen::Ref<const RealVect> & max_q,
                                     const Eigen::Ref<const RealVect> & power_factor,
                                     const Eigen::Ref<const Eigen::VectorXi> & bus_id)
{
    // active power is derived (owned by the HvdcLineContainer), it starts at 0
    const auto p_init = RealVect::Zero(bus_id.size());
    init_osc_pq(p_init, q_setpoint_mvar, bus_id, "converter_stations");

    const int size = nb();
    check_size(type, size, "type");
    check_size(loss_factor, size, "loss_factor");
    check_size(voltage_regulator_on, size, "voltage_regulator_on");
    check_size(target_vm_pu, size, "target_vm_pu");
    check_size(min_q, size, "min_q");
    check_size(max_q, size, "max_q");
    check_size(power_factor, size, "power_factor");

    type_ = IntVect(size);
    for(int station_id = 0; station_id < size; ++station_id){
        if((type[station_id] != ConverterType::VSC) && (type[station_id] != ConverterType::LCC)){
            std::ostringstream exc_;
            exc_ << "ConverterStationContainer::init: unknown converter type ";
            exc_ << type[station_id];
            exc_ << " for station ";
            exc_ << station_id;
            exc_ << " (0 = VSC, 1 = LCC).";
            throw std::runtime_error(exc_.str());
        }
        type_(station_id) = type[station_id];
        if((type[station_id] == ConverterType::LCC) && voltage_regulator_on[station_id]){
            std::ostringstream exc_;
            exc_ << "ConverterStationContainer::init: LCC station ";
            exc_ << station_id;
            exc_ << " cannot regulate voltage.";
            throw std::runtime_error(exc_.str());
        }
        if(min_q(station_id) > max_q(station_id)){
            std::ostringstream exc_;
            exc_ << "ConverterStationContainer::init: Impossible to initialize converter station with min_q being above max_q for station ";
            exc_ << station_id;
            throw std::runtime_error(exc_.str());
        }
    }
    loss_factor_ = loss_factor;
    voltage_regulator_on_ = voltage_regulator_on;
    target_vm_pu_ = target_vm_pu;
    min_q_ = min_q;
    max_q_ = max_q;
    power_factor_ = power_factor;
    reset_results();
}

ConverterStationContainer::StateRes ConverterStationContainer::get_state() const
{
    std::vector<int> type(type_.begin(), type_.end());
    std::vector<real_type> loss_factor(loss_factor_.begin(), loss_factor_.end());
    std::vector<real_type> vm_pu(target_vm_pu_.begin(), target_vm_pu_.end());
    std::vector<real_type> min_q(min_q_.begin(), min_q_.end());
    std::vector<real_type> max_q(max_q_.begin(), max_q_.end());
    std::vector<real_type> power_factor(power_factor_.begin(), power_factor_.end());
    ConverterStationContainer::StateRes res(get_osc_pq_state(),
                                            type,
                                            loss_factor,
                                            voltage_regulator_on_,
                                            vm_pu,
                                            min_q,
                                            max_q,
                                            power_factor);
    return res;
}

void ConverterStationContainer::set_state(ConverterStationContainer::StateRes & my_state)
{
    set_osc_pq_state(std::get<StateResIdx::OSC_PQ_STATE>(my_state));
    std::vector<int> & type = std::get<StateResIdx::TYPE>(my_state);
    std::vector<real_type> & loss_factor = std::get<StateResIdx::LOSS_FACTOR>(my_state);
    std::vector<bool> & voltage_regulator_on = std::get<StateResIdx::VREG_ON>(my_state);
    std::vector<real_type> & vm_pu = std::get<StateResIdx::TARGET_VM_PU>(my_state);
    std::vector<real_type> & min_q = std::get<StateResIdx::MIN_Q>(my_state);
    std::vector<real_type> & max_q = std::get<StateResIdx::MAX_Q>(my_state);
    std::vector<real_type> & power_factor = std::get<StateResIdx::POWER_FACTOR>(my_state);

    const auto size = nb();
    check_size(type, size, "type");
    check_size(loss_factor, size, "loss_factor");
    check_size(voltage_regulator_on, size, "voltage_regulator_on");
    check_size(vm_pu, size, "vm_pu");
    check_size(min_q, size, "min_q");
    check_size(max_q, size, "max_q");
    check_size(power_factor, size, "power_factor");

    type_ = IntVect::Map(type.data(), type.size());
    loss_factor_ = RealVect::Map(loss_factor.data(), loss_factor.size());
    voltage_regulator_on_ = voltage_regulator_on;
    target_vm_pu_ = RealVect::Map(vm_pu.data(), vm_pu.size());
    min_q_ = RealVect::Map(min_q.data(), min_q.size());
    max_q_ = RealVect::Map(max_q.data(), max_q.size());
    power_factor_ = RealVect::Map(power_factor.data(), power_factor.size());
    reset_results();
}

void ConverterStationContainer::set_station_p(int station_id, real_type p_mw, DualAlgoControl & solver_control)
{
    [[maybe_unused]] bool my_status = status_.at(station_id);  // also checks station_id is in range
    this->_change_p(station_id, p_mw, my_status, solver_control);
    if (abs(target_p_mw_(station_id) - p_mw) > _tol_equal_float) {
        target_p_mw_(station_id) = p_mw;
    }
    if(is_lcc(station_id)){
        // LCC always consumes Q = |P| * tan(acos(power_factor)) (generator sign convention)
        const real_type new_q = -abs(p_mw) * std::tan(std::acos(power_factor_(station_id)));
        if (abs(target_q_mvar_(station_id) - new_q) > _tol_equal_float) {
            solver_control.ac_algo_controler().tell_recompute_sbus();
            solver_control.dc_algo_controler().tell_recompute_sbus();
            target_q_mvar_(station_id) = new_q;
        }
    }
}

void ConverterStationContainer::change_v(int station_id, real_type new_v_pu, DualAlgoControl & solver_control)
{
    bool my_status = status_.at(station_id);
    if(!my_status)
    {
        std::ostringstream exc_;
        exc_ << "ConverterStationContainer::change_v: Impossible to change the voltage setpoint of a disconnected converter station (check station id ";
        exc_ << station_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    if (abs(target_vm_pu_(station_id) - new_v_pu) > _tol_equal_float)
    {
        solver_control.ac_algo_controler().tell_v_changed();
        solver_control.dc_algo_controler().tell_v_changed();
        target_vm_pu_(station_id) = new_v_pu;
    }
}

void ConverterStationContainer::fillSbus_station(Eigen::Ref<CplxVect> Sbus,
                                                 const SolverBusIdVect & id_grid_to_solver,
                                                 bool /*ac*/,
                                                 const std::vector<bool> & skip_p) const
{
    const int nb_station = nb();
    GlobalBusId bus_id_me;
    SolverBusId bus_id_solver;
    cplx_type tmp;
    for(int station_id = 0; station_id < nb_station; ++station_id){
        //  i don't do anything if the station is disconnected
        if(!status_[station_id]) continue;

        bus_id_me = bus_id_(station_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "ConverterStationContainer::fillSbus_station: Converter station with id ";
            exc_ << station_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "ConverterStationContainer::fillSbus_station: Converter station with id ";
            exc_ << station_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        tmp = {skip_p[station_id] ? 0. : target_p_mw_(station_id), 0.};
        if(!voltage_regulator_on_[station_id]){
            // station is pq if voltage regulation is off (VSC q setpoint or LCC consumption)
            tmp += my_i * target_q_mvar_(station_id);
        }
        Sbus.coeffRef(bus_id_solver.cast_int()) += tmp;
    }
}

void ConverterStationContainer::fillpv(std::vector<int> & bus_pv,
                                       std::vector<bool> & has_bus_been_added,
                                       const SolverBusIdVect & slack_bus_id_solver,
                                       const SolverBusIdVect & id_grid_to_solver) const
{
    const int nb_station = nb();
    GlobalBusId bus_id_me;
    SolverBusId bus_id_solver;
    for(int station_id = 0; station_id < nb_station; ++station_id){
        if(!status_[station_id]) continue;
        if (!voltage_regulator_on_[station_id]) continue;  // station is purposedly not pv

        bus_id_me = bus_id_(station_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "ConverterStationContainer::fillpv: Converter station with id ";
            exc_ << station_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "ConverterStationContainer::fillpv: Converter station with id ";
            exc_ << station_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }

        if(is_in_vect(bus_id_solver.cast_int(), slack_bus_id_solver.to_int_vector())) continue;  // slack bus is not PV
        if(has_bus_been_added[bus_id_solver.cast_int()]) continue; // i already added this bus
        bus_pv.push_back(bus_id_solver.cast_int());
        has_bus_been_added[bus_id_solver.cast_int()] = true;  // don't add it a second time
    }
}

void ConverterStationContainer::init_q_vector(int /*nb_bus*/,
                                              Eigen::Ref<Eigen::VectorXi> total_gen_per_bus,
                                              Eigen::Ref<RealVect> total_q_min_per_bus,
                                              Eigen::Ref<RealVect> total_q_max_per_bus) const
{
    const int nb_station = nb();
    for(int station_id = 0; station_id < nb_station; ++station_id)
    {
        if(!status_[station_id]) continue;
        if (!voltage_regulator_on_[station_id]) continue;  // station is purposedly not pv

        const GlobalBusId bus_id = bus_id_(station_id);
        total_q_min_per_bus(bus_id.cast_int()) += min_q_(station_id);
        total_q_max_per_bus(bus_id.cast_int()) += max_q_(station_id);
        total_gen_per_bus(bus_id.cast_int()) += 1;
    }
}

void ConverterStationContainer::set_q(const Eigen::Ref<const RealVect> & reactive_mismatch,
                                      const SolverBusIdVect & id_grid_to_solver,
                                      bool ac,
                                      const Eigen::Ref<const Eigen::VectorXi> & total_gen_per_bus,
                                      const Eigen::Ref<const RealVect> & total_q_min_per_bus,
                                      const Eigen::Ref<const RealVect> & total_q_max_per_bus)
{
    const int nb_station = nb();
    if(!ac){
        // do not consider Q values in dc mode
        for(int station_id = 0; station_id < nb_station; ++station_id) res_q_(station_id) = 0.;
        return;
    }

    real_type eps_q = 1e-8;
    for(int station_id = 0; station_id < nb_station; ++station_id)
    {
        if(!status_[station_id]){
            // set at 0 for disconnected stations
            res_q_(station_id) = 0.;
            continue;
        }
        if (!voltage_regulator_on_[station_id]){
            // station is purposedly not pv, so output MVAr = input MVAr (just like a load)
            res_q_(station_id) = target_q_mvar_(station_id);
            continue;
        }
        const GlobalBusId bus_id = bus_id_(station_id);
        const SolverBusId bus_solver = id_grid_to_solver[bus_id.cast_int()];
        // TODO DEBUG MODE: check that the bus is correct!
        real_type q_to_absorb = reactive_mismatch[bus_solver.cast_int()];
        real_type max_q_me = max_q_(station_id);
        real_type min_q_me = min_q_(station_id);
        real_type max_q_bus = total_q_max_per_bus(bus_id.cast_int());
        real_type min_q_bus = total_q_min_per_bus(bus_id.cast_int());
        real_type nb_gen_with_me = static_cast<real_type>(total_gen_per_bus(bus_id.cast_int()));
        real_type real_q;
        if(nb_gen_with_me == 1.){
            real_q = q_to_absorb;
        }else{
            real_type ratio = (max_q_me - min_q_me + eps_q) / (max_q_bus - min_q_bus + nb_gen_with_me * eps_q) ;
            real_q = q_to_absorb * ratio ;
        }
        res_q_(station_id) = real_q;
    }
}

void ConverterStationContainer::get_vm_for_dc(Eigen::Ref<RealVect> Vm)
{
    const int nb_station = nb();
    GlobalBusId bus_id_me;
    for(int station_id = 0; station_id < nb_station; ++station_id){
        if(!status_[station_id]) continue;
        if (!voltage_regulator_on_[station_id]) continue;  // station is purposedly not pv

        bus_id_me = bus_id_(station_id);
        real_type tmp = target_vm_pu_(station_id);
        if(abs(tmp) > _tol_equal_float) Vm(bus_id_me.cast_int()) = tmp;
    }
}

void ConverterStationContainer::set_vm(Eigen::Ref<CplxVect> V, const SolverBusIdVect & id_grid_to_solver) const
{
    const int nb_station = nb();
    GlobalBusId bus_id_me;
    SolverBusId bus_id_solver;
    for(int station_id = 0; station_id < nb_station; ++station_id){
        if(!status_[station_id]) continue;
        if (!voltage_regulator_on_[station_id]) continue;  // station is purposedly not pv

        bus_id_me = bus_id_(station_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "ConverterStationContainer::set_vm: Converter station with id ";
            exc_ << station_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "ConverterStationContainer::set_vm: Converter station with id ";
            exc_ << station_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        // scale the input V such that abs(V) = Vm for this station
        real_type tmp = std::abs(V(bus_id_solver.cast_int()));
        if(abs(tmp) < _tol_equal_float)
        {
            // if it was 0. i force it to 1. (otherwise the rest of the computation would make it O. still)
            V(bus_id_solver.cast_int()) = 1.0;
            tmp = 1.0;
        }
        tmp = 1.0 / tmp;
        tmp *= target_vm_pu_(station_id);
        V(bus_id_solver.cast_int()) *= tmp;
    }
}

void ConverterStationContainer::_compute_results(
    const Eigen::Ref<const RealVect> & /*Va*/,
    const Eigen::Ref<const RealVect> & /*Vm*/,
    const Eigen::Ref<const CplxVect> & /*V*/,
    const SolverBusIdVect & /*id_grid_to_solver*/,
    const Eigen::Ref<const RealVect> & /*bus_vn_kv*/,
    real_type /*sn_mva*/,
    bool ac)
{
    set_osc_pq_res_p();
    if(ac){
        const int nb_station = nb();
        for(int station_id = 0; station_id < nb_station; ++station_id)
        {
            if(!status_[station_id]){
                // turned off station does not have q
                res_q_[station_id] = 0.;
                continue;
            }
            if(voltage_regulator_on_[station_id]) continue;  // filled later by set_q
            res_q_(station_id) = target_q_mvar_(station_id);
        }
    }else{
        // nothing special to do here
        set_osc_pq_res_q(ac);
    }
}

void ConverterStationContainer::_change_p(int station_id, real_type new_p, bool /*my_status*/, DualAlgoControl & solver_control)
{
    if (abs(target_p_mw_(station_id) - new_p) > _tol_equal_float) {
        solver_control.ac_algo_controler().tell_recompute_sbus();
        solver_control.dc_algo_controler().tell_recompute_sbus();
    }
    // turned off stations (p == 0) are not pv: if the active power changes,
    // the list of pv buses may change (legacy dcline behaviour, cf `turnedoff_no_pv`)
    bool pseudo_off_before = abs(target_p_mw_(station_id)) < _tol_equal_float;
    bool pseudo_off_now = abs(new_p) < _tol_equal_float;
    if((pseudo_off_before && !pseudo_off_now) ||
       (!pseudo_off_before && pseudo_off_now)){
        solver_control.ac_algo_controler().tell_pv_changed();
        solver_control.dc_algo_controler().tell_pv_changed();
    }
}

bool ConverterStationContainer::_deactivate(int el_id, DualAlgoControl & solver_control) {
    if(!status_[el_id]) return false;  // nothing to do if it was already deactivated
    solver_control.ac_algo_controler().tell_recompute_sbus();
    solver_control.dc_algo_controler().tell_recompute_sbus();
    solver_control.ac_algo_controler().tell_pv_changed();
    solver_control.dc_algo_controler().tell_pv_changed();
    return true;
}

bool ConverterStationContainer::_reactivate(int el_id, DualAlgoControl & solver_control) {
    if(status_[el_id]) return false;  // nothing to do if station already connected
    solver_control.ac_algo_controler().tell_recompute_sbus();
    solver_control.dc_algo_controler().tell_recompute_sbus();
    solver_control.ac_algo_controler().tell_pv_changed();
    solver_control.dc_algo_controler().tell_pv_changed();
    return true;
}

bool ConverterStationContainer::_change_bus(int el_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control, int /*nb_bus*/) {
    if(bus_id_(el_id) == new_bus_id) return false;  // nothing to do if the bus did not changed
    solver_control.ac_algo_controler().tell_recompute_sbus();
    solver_control.dc_algo_controler().tell_recompute_sbus();
    solver_control.ac_algo_controler().tell_one_el_changed_bus();
    solver_control.dc_algo_controler().tell_one_el_changed_bus();
    if(voltage_regulator_on_[el_id]) {
        solver_control.ac_algo_controler().tell_pv_changed();
        solver_control.dc_algo_controler().tell_pv_changed();
    }
    return true;
}

} // namespace ls2g
