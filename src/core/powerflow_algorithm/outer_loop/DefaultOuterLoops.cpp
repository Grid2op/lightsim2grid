// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DefaultOuterLoops.hpp"

#include <sstream>
#include <stdexcept>

#include "LSGrid.hpp"
#include "DistributedSlackLoop.hpp"
#include "HvdcAcEmulationLimitsLoop.hpp"
#include "PhaseControlLoop.hpp"
#include "ReactiveLimitsLoop.hpp"
#include "ShuntVoltageControlLoop.hpp"
#include "TransformerVoltageControlLoop.hpp"
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
    res.push_back(std::make_shared<ReactiveLimitsLoop>());
    return res;
}

namespace {
void check_param_count(const std::string & name, const AlgoConfig & params,
                       std::size_t nb_int, std::size_t nb_real)
{
    if (params.int_params.size() == nb_int && params.real_params.size() == nb_real) return;
    std::ostringstream exc_;
    exc_ << "make_outer_loop: a '" << name << "' loop has " << nb_int << " integer and " << nb_real
         << " real parameters, got " << params.int_params.size() << " and "
         << params.real_params.size() << ".";
    throw std::runtime_error(exc_.str());
}
}  // namespace

std::shared_ptr<BaseOuterLoop> make_outer_loop(const std::string & name, const AlgoConfig & params)
{
    // the inverse of each loop's _get_params
    if (name == "DistributedSlack") {
        check_param_count(name, params, 1, 3);
        DistributedSlackLoop::Params p;
        p.fail_on_residue = params.int_params[0] != 0;
        p.slack_bus_p_max_mismatch_mw = static_cast<real_type>(params.real_params[0]);
        p.p_residue_eps_mw = static_cast<real_type>(params.real_params[1]);
        p.moved_fraction = static_cast<real_type>(params.real_params[2]);
        return std::make_shared<DistributedSlackLoop>(p);
    }
    if (name == "ReactiveLimits") {
        check_param_count(name, params, 2, 6);
        ReactiveLimitsLoop::Params p;
        p.max_pq_pv_switch = params.int_params[0];
        p.robust_mode = params.int_params[1] != 0;
        p.max_reactive_power_mismatch = static_cast<real_type>(params.real_params[0]);
        p.mismatch_base_mva = static_cast<real_type>(params.real_params[1]);
        p.min_realistic_voltage = static_cast<real_type>(params.real_params[2]);
        p.max_realistic_voltage = static_cast<real_type>(params.real_params[3]);
        p.realistic_voltage_margin = static_cast<real_type>(params.real_params[4]);
        p.robust_restart_vm_pu = static_cast<real_type>(params.real_params[5]);
        return std::make_shared<ReactiveLimitsLoop>(p);
    }
    if (name == "TransformerVoltageControl") {
        check_param_count(name, params, 1, 2);
        TransformerVoltageControlLoop::Params p;
        p.use_initial_tap_position = params.int_params[0] != 0;
        p.max_controlled_nominal_voltage = static_cast<real_type>(params.real_params[0]);
        p.min_target_deadband_kv = static_cast<real_type>(params.real_params[1]);
        return std::make_shared<TransformerVoltageControlLoop>(p);
    }
    std::shared_ptr<BaseOuterLoop> res;
    if (name == "AcHvdcAcEmulationLimits") res = std::make_shared<HvdcAcEmulationLimitsLoop>();
    else if (name == "VoltageMonitoring") res = std::make_shared<VoltageMonitoringLoop>();
    else if (name == "PhaseControl") res = std::make_shared<PhaseControlLoop>();
    else if (name == "ShuntVoltageControl") res = std::make_shared<ShuntVoltageControlLoop>();
    else {
        std::ostringstream exc_;
        exc_ << "make_outer_loop: unknown outer loop '" << name << "'.";
        throw std::runtime_error(exc_.str());
    }
    check_param_count(name, params, 0, 0);
    return res;
}

void outer_detect(const std::vector<std::shared_ptr<BaseOuterLoop> > & loops,
                  OuterContext & ctx, std::vector<LimitViolation> & out)
{
    ctx.injections = nullptr;
    ctx.controls = nullptr;
    for (const auto & loop : loops) {
        // not is_needed: it reads what an outer-loop solve set up (the monitors a plan holds
        // for VoltageMonitoring, ...), which a plain solve does not; each detect decides
        // itself what the grid gives it to check
        if (loop == nullptr) continue;
        std::unique_ptr<BaseOuterLoop> fresh = loop->clone();
        fresh->initialize(ctx);
        fresh->detect(ctx, out);
    }
}

}  // namespace ls2g
