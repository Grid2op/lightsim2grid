// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DefaultOuterLoops.hpp"

#include "LSGrid.hpp"
#include "DistributedSlackLoop.hpp"
#include "HvdcAcEmulationLimitsLoop.hpp"
#include "VoltageMonitoringLoop.hpp"

namespace ls2g {

std::vector<std::shared_ptr<BaseOuterLoop> > make_default_outer_loops(const LSGrid & /*grid*/)
{
    // OpenLoadFlow's order: DistributedSlack, (FreezingHvdcACEmulation), AcHvdcAcEmulationLimits,
    // (AreaInterchangeControl), (SecondaryVoltageControl), VoltageMonitoring, ReactiveLimits,
    // PhaseControl, TransformerVoltageControl, (TransformerReactivePowerControl),
    // ShuntVoltageControl, (AutomationSystem). Each loop joins this list as it is implemented.
    std::vector<std::shared_ptr<BaseOuterLoop> > res;
    res.push_back(std::make_shared<DistributedSlackLoop>());
    res.push_back(std::make_shared<HvdcAcEmulationLimitsLoop>());
    res.push_back(std::make_shared<VoltageMonitoringLoop>());
    return res;
}

}  // namespace ls2g
