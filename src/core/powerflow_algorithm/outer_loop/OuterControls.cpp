// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "OuterControls.hpp"

namespace ls2g {

namespace {

template<class Control>
Control * find_control(const std::map<int, Control *> & of, int key)
{
    const auto it = of.find(key);
    return it == of.end() ? nullptr : it->second;
}

}  // namespace

// out-of-class definition: the C++14 build odr-uses it (vector::assign takes a const reference)
constexpr int HvdcRegimeControl::KEEP;

BusVoltageControl * OuterControls::reserve_bus_voltage(int bus)
{
    BusVoltageControl * known = bus_voltage(bus);
    if (known != nullptr) return known;
    bus_voltage_store_.push_back(BusVoltageControl(bus));
    BusVoltageControl * res = &bus_voltage_store_.back();
    bus_voltage_of_[bus] = res;
    return res;
}

BusVoltageControl * OuterControls::bus_voltage(int bus)
{
    const auto it = bus_voltage_of_.find(bus);
    return it == bus_voltage_of_.end() ? nullptr : it->second;
}

const BusVoltageControl * OuterControls::bus_voltage(int bus) const
{
    const auto it = bus_voltage_of_.find(bus);
    return it == bus_voltage_of_.end() ? nullptr : it->second;
}

VoltageControllerHold * OuterControls::reserve_controller_hold(int controller)
{
    VoltageControllerHold * known = controller_hold(controller);
    if (known != nullptr) return known;
    controller_hold_store_.push_back(VoltageControllerHold(controller));
    VoltageControllerHold * res = &controller_hold_store_.back();
    controller_hold_of_[controller] = res;
    return res;
}

VoltageControllerHold * OuterControls::controller_hold(int controller)
{
    const auto it = controller_hold_of_.find(controller);
    return it == controller_hold_of_.end() ? nullptr : it->second;
}

const VoltageControllerHold * OuterControls::controller_hold(int controller) const
{
    const auto it = controller_hold_of_.find(controller);
    return it == controller_hold_of_.end() ? nullptr : it->second;
}

StandbySvcControl * OuterControls::reserve_standby_svc(int svc)
{
    StandbySvcControl * known = standby_svc(svc);
    if (known != nullptr) return known;
    standby_svc_store_.push_back(StandbySvcControl(svc));
    StandbySvcControl * res = &standby_svc_store_.back();
    standby_svc_of_[svc] = res;
    return res;
}

HvdcRegimeControl * OuterControls::reserve_hvdc_regime(int line)
{
    HvdcRegimeControl * known = hvdc_regime(line);
    if (known != nullptr) return known;
    hvdc_regime_store_.push_back(HvdcRegimeControl(line));
    HvdcRegimeControl * res = &hvdc_regime_store_.back();
    hvdc_regime_of_[line] = res;
    return res;
}

PhaseShifterControl * OuterControls::reserve_phase_shifter(int trafo, bool solves_shift)
{
    if (!_supports_phase_shifters()) return nullptr;
    PhaseShifterControl * known = phase_shifter(trafo);
    if (known != nullptr) {
        known->solves_shift_ = known->solves_shift_ || solves_shift;
        return known;
    }
    phase_shifter_store_.push_back(PhaseShifterControl(this, trafo, solves_shift));
    PhaseShifterControl * res = &phase_shifter_store_.back();
    phase_shifter_of_[trafo] = res;
    phase_shifter_order_.push_back(res);
    return res;
}

PhaseShifterControl * OuterControls::phase_shifter(int trafo) { return find_control(phase_shifter_of_, trafo); }
const PhaseShifterControl * OuterControls::phase_shifter(int trafo) const { return find_control(phase_shifter_of_, trafo); }

bool PhaseShifterControl::handled() const { return owner_->_phase_handled(trafo_); }
real_type PhaseShifterControl::shift() const { return owner_->_phase_shift(trafo_); }
int PhaseShifterControl::position() const { return owner_->_phase_position(trafo_); }
void PhaseShifterControl::current(int side, real_type & i_pu, real_type & di_da) const
{
    owner_->_phase_current(trafo_, side, i_pu, di_da);
}

StandbySvcControl * OuterControls::standby_svc(int svc) { return find_control(standby_svc_of_, svc); }
const StandbySvcControl * OuterControls::standby_svc(int svc) const { return find_control(standby_svc_of_, svc); }
HvdcRegimeControl * OuterControls::hvdc_regime(int line) { return find_control(hvdc_regime_of_, line); }
const HvdcRegimeControl * OuterControls::hvdc_regime(int line) const { return find_control(hvdc_regime_of_, line); }

void OuterControls::clear_reservations()
{
    bus_voltage_of_.clear();
    bus_voltage_store_.clear();
    controller_hold_of_.clear();
    controller_hold_store_.clear();
    standby_svc_of_.clear();
    standby_svc_store_.clear();
    hvdc_regime_of_.clear();
    hvdc_regime_store_.clear();
    phase_shifter_of_.clear();
    phase_shifter_order_.clear();
    phase_shifter_store_.clear();
}

void OuterControls::reset_states()
{
    for (auto & control : bus_voltage_store_) control._reset();
    for (auto & hold : controller_hold_store_) hold.release();
    for (auto & svc : standby_svc_store_) svc._reset();
    for (auto & regime : hvdc_regime_store_) regime._reset();
    for (auto & shifter : phase_shifter_store_) shifter._reset();
    pending_vm_.clear();
    suspended_.clear();
}

void OuterControls::take_pending_vm(std::vector<int> & buses, std::vector<real_type> & vm)
{
    buses.clear();
    vm.clear();
    for (const auto & bv : pending_vm_) {
        buses.push_back(bv.first);
        vm.push_back(bv.second);
    }
    pending_vm_.clear();
}

std::vector<int> OuterControls::switchable_buses() const
{
    std::vector<int> res = caller_switchable_;
    for (const auto & bc : bus_voltage_of_) res.push_back(bc.first);
    return res;
}

std::vector<int> OuterControls::pinned_buses() const
{
    std::vector<int> res = caller_pinned_;
    for (const auto & bc : bus_voltage_of_) {
        if (!bc.second->is_pq()) res.push_back(bc.first);
    }
    return res;
}

}  // namespace ls2g

namespace ls2g {

std::vector<real_type> OuterControls::held_q() const
{
    std::vector<real_type> res;
    if (controller_hold_of_.empty() || controller_hold_of_.rbegin()->first < 0) return res;
    res.assign(static_cast<std::size_t>(controller_hold_of_.rbegin()->first) + 1,
               std::numeric_limits<real_type>::quiet_NaN());
    for (const auto & ch : controller_hold_of_) {
        if (ch.first >= 0) res[static_cast<std::size_t>(ch.first)] = ch.second->q();
    }
    return res;
}

std::vector<real_type> OuterControls::svc_target_vm() const
{
    std::vector<real_type> res;
    if (standby_svc_of_.empty() || standby_svc_of_.rbegin()->first < 0) return res;
    res.assign(static_cast<std::size_t>(standby_svc_of_.rbegin()->first) + 1,
               std::numeric_limits<real_type>::quiet_NaN());
    for (const auto & sc : standby_svc_of_) {
        if (sc.first >= 0) res[static_cast<std::size_t>(sc.first)] = sc.second->target_vm();
    }
    return res;
}

std::vector<int> OuterControls::hvdc_regimes() const
{
    std::vector<int> res;
    if (hvdc_regime_of_.empty() || hvdc_regime_of_.rbegin()->first < 0) return res;
    res.assign(static_cast<std::size_t>(hvdc_regime_of_.rbegin()->first) + 1, HvdcRegimeControl::KEEP);
    for (const auto & hr : hvdc_regime_of_) {
        if (hr.first >= 0) res[static_cast<std::size_t>(hr.first)] = hr.second->regime();
    }
    return res;
}

}  // namespace ls2g
