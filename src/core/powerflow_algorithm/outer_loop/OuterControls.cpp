// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "OuterControls.hpp"

namespace ls2g {

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

void OuterControls::clear_reservations()
{
    bus_voltage_of_.clear();
    bus_voltage_store_.clear();
}

void OuterControls::reset_states()
{
    for (auto & control : bus_voltage_store_) control._reset();
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
