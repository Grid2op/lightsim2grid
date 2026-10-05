// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "OuterLoopDriver.hpp"

#include <cmath>

#include "LSGrid.hpp"

namespace ls2g {

void OuterContext::record_bus(const char * action, bool taken, int solver_bus, LimitViolationType reason,
                              real_type value, real_type limit) const
{
    if (trace == nullptr) return;
    const GlobalBusIdVect & solver_to_me = grid->id_ac_solver_to_me();
    const int me = solver_bus >= 0 && solver_bus < static_cast<int>(solver_to_me.size())
                   ? solver_to_me[solver_bus].cast_int() : -1;
    record(action, taken, ViolationElementType::BUS, me, reason, value, limit);
}

bool is_state_unrealistic(const LSGrid & grid,
                          const CplxVect & V,
                          const std::vector<bool> & vm_unknown,
                          const OuterLoopDriverParams & params)
{
    const GlobalBusIdVect & solver_to_me = grid.id_ac_solver_to_me();
    Eigen::Ref<const RealVect> vn_kv = grid.get_bus_vn_kv();
    const Eigen::Index nb_bus = V.size();
    for(Eigen::Index bus = 0; bus < nb_bus; ++bus){
        if(static_cast<std::size_t>(bus) >= vm_unknown.size() || !vm_unknown[bus]) continue;
        const real_type vm = std::abs(V(bus));
        if(vm >= params.min_realistic_voltage && vm <= params.max_realistic_voltage) continue;
        const int me = solver_to_me[static_cast<int>(bus)].cast_int();
        if(vn_kv(me) >= params.min_nominal_voltage_realistic_check) return true;
    }
    return false;
}

}  // namespace ls2g
