// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BUS_Q_PLAN_H
#define BUS_Q_PLAN_H

#include <string>
#include <utility>
#include <vector>

namespace ls2g {

// The plan of the bus Q check, as plain data (see BusQCheck.hpp, which builds it): apart
// from the grid, so the outer loops' context can carry one without including LSGrid.
namespace bus_q_check {

/// One bus whose voltage is held by machines, and everything a row needs to check it.
/// Built once per compute() by `build_bus_q_plan`; nothing here varies from row to row.
struct BusQEntry
{
    int bus_solver = -1;          ///< the bus, solver numbering
    int bus_grid = -1;            ///< the same bus, grid numbering (what is reported)
    /// the generators holding this bus: a row that disconnects one takes its limits out of
    /// the sum (`min_q_mvar` / `max_q_mvar`, fixed)
    std::vector<int> gen_ids;
    /// the voltage-regulating storage units holding it (always through the PV path, so
    /// never part of `ctrl_pos`): a fixed [min_q, max_q] in MVAr, generator convention
    std::vector<int> storage_ids;
    /// the hvdc converter stations holding it, as (hvdc line id, side): also a fixed
    /// [min_q, max_q] in MVAr
    std::vector<std::pair<int, int> > station_ids;
    /// the voltage-mode SVCs holding it: a susceptance range, so their MVAr capability is
    /// re-evaluated at every row's own voltage
    std::vector<int> svc_ids;
    /// positions in the algorithm's controller list of the controllers standing on this
    /// bus (`BaseAlgo::get_controller_q`), empty for an ordinary PV bus
    std::vector<int> ctrl_pos;
    std::string sub_name;         ///< substation of the bus, empty if the grid has no names
};

struct BusQPlan
{
    std::vector<BusQEntry> buses;
    /// whether any bus of the plan needs the controller list at all: false on an ordinary
    /// grid of local PV machines, and then a row never asks for it
    bool needs_controller_q = false;

    bool empty() const { return buses.empty(); }
    void clear() { buses.clear(); needs_controller_q = false; }
};

}  // namespace bus_q_check

}  // namespace ls2g

#endif  // BUS_Q_PLAN_H
