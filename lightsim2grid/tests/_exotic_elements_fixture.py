# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Shared "one of everything exotic" test grid: IEEE14 augmented with a static var
compensator, a storage unit, three HVDC lines (VSC with angle-droop, VSC without,
LCC without), and a phase-shifting transformer.

Built for PR #140 self-review: ``test_binary_serialization.py`` /
``test_pickleable.py`` used a full ``l2rpn_idf_2023`` grid2op environment (118
substations, slow to build) for every one of their per-element-type round-trip
tests, but that environment carries neither an SVC nor any HVDC line -- so
``test_binary_hvdc`` / ``test_binary_svcs`` (and their pickle counterparts) were
round-tripping empty containers, testing nothing about those element types. This
grid is built once (no grid2op environment at all -- ``init_from_pypowsybl``
straight into an ``LSGrid``) and reused wherever SVC / HVDC coverage is needed.

Deviation from a literal "1 HVDC LCC with droop, 1 HVDC LCC without droop, 1 HVDC
VSC" reading: lightsim2grid's ``HvdcLineContainer`` only supports the angle-droop
(AC emulation) extension for VSC-VSC lines (``HvdcLineContainer::init`` raises
otherwise) -- an LCC line with droop is not a configuration lightsim2grid can ever
load. Substituted for equivalent regime coverage: HVDC VSC-VSC *with* droop, HVDC
VSC-VSC *without* droop, and HVDC LCC-LCC (never has droop).
"""

import pandas as pd

import pypowsybl as pp

from lightsim2grid.network import init_from_pypowsybl


def build_exotic_elements_network() -> "pp.network.Network":
    """Return the raw pypowsybl ``Network`` (not yet solved / converted), for
    tests that want to inspect or tweak it before conversion."""
    n = pp.network.create_ieee14()

    # SVC, voltage-regulating with a slope, on the (otherwise load-only) leaf bus B14
    n.create_static_var_compensators(id="SVC1", voltage_level_id="VL14", bus_id="B14",
                                     regulation_mode="VOLTAGE", regulating=True,
                                     target_v=12.0, target_q=0.0, b_min=-0.02, b_max=0.02)
    n.create_extensions("voltagePerReactivePowerControl", id="SVC1", slope=0.01)

    # storage unit (battery), on the (otherwise load-only) leaf bus B13
    n.create_batteries(id="BATT1", voltage_level_id="VL13", bus_id="B13",
                       max_p=20.0, min_p=-20.0, target_p=5.0, target_q=1.0)

    # HVDC VSC-VSC, angle-droop (AC emulation) enabled
    n.create_vsc_converter_stations(id="VSC1", voltage_level_id="VL10", bus_id="B10", loss_factor=1.1,
                                    voltage_regulator_on=True, target_v=12.0, target_q=0.0)
    n.create_vsc_converter_stations(id="VSC2", voltage_level_id="VL11", bus_id="B11", loss_factor=1.1,
                                    voltage_regulator_on=True, target_v=12.0, target_q=0.0)
    n.create_hvdc_lines(id="HVDC_VSC_DROOP", converter_station1_id="VSC1", converter_station2_id="VSC2",
                        r=1.0, nominal_v=12.0, converters_mode="SIDE_1_RECTIFIER_SIDE_2_INVERTER",
                        target_p=2.0, max_p=20.0)
    n.create_extensions("hvdcAngleDroopActivePowerControl", id="HVDC_VSC_DROOP",
                        p0=1.0, droop=180.0, enabled=True)

    # HVDC VSC-VSC, fixed setpoint (no droop)
    n.create_vsc_converter_stations(id="VSC3", voltage_level_id="VL12", bus_id="B12", loss_factor=1.1,
                                    voltage_regulator_on=True, target_v=12.0, target_q=0.0)
    n.create_vsc_converter_stations(id="VSC4", voltage_level_id="VL9", bus_id="B9", loss_factor=1.1,
                                    voltage_regulator_on=True, target_v=12.0, target_q=0.0)
    n.create_hvdc_lines(id="HVDC_VSC_NODROOP", converter_station1_id="VSC3", converter_station2_id="VSC4",
                        r=1.0, nominal_v=12.0, converters_mode="SIDE_1_RECTIFIER_SIDE_2_INVERTER",
                        target_p=2.0, max_p=20.0)

    # HVDC LCC-LCC, fixed setpoint (LCC never supports droop, see module docstring)
    n.create_lcc_converter_stations(id="LCC1", voltage_level_id="VL14", bus_id="B14",
                                    loss_factor=1.1, power_factor=0.9)
    n.create_lcc_converter_stations(id="LCC2", voltage_level_id="VL13", bus_id="B13",
                                    loss_factor=1.1, power_factor=0.9)
    n.create_hvdc_lines(id="HVDC_LCC", converter_station1_id="LCC1", converter_station2_id="LCC2",
                        r=1.0, nominal_v=12.0, converters_mode="SIDE_1_RECTIFIER_SIDE_2_INVERTER",
                        target_p=1.0, max_p=20.0)

    # phase-shifting transformer: add a (non-regulating) phase tap changer to the
    # existing T4-7-1 2-winding transformer, 3 steps (-5deg / 0deg / +5deg)
    ptc_df = pd.DataFrame.from_records(
        index="id", columns=["id", "tap", "low_tap", "regulating", "regulation_mode"],
        data=[("T4-7-1", 1, 0, False, "CURRENT_LIMITER")])
    steps_df = pd.DataFrame.from_records(
        index="id", columns=["id", "b", "g", "r", "x", "rho", "alpha"],
        data=[("T4-7-1", 0., 0., 0., 0., 1.0, -5.0),
              ("T4-7-1", 0., 0., 0., 0., 1.0, 0.0),
              ("T4-7-1", 0., 0., 0., 0., 1.0, 5.0)])
    n.create_phase_tap_changers(ptc_df, steps_df)

    return n


def build_exotic_elements_grid(solve: bool = True):
    """Return the grid above, converted to an :class:`LSGrid` via
    ``init_from_pypowsybl`` (slack pinned at ``B1-G``). Solved (AC, flat start)
    by default, matching what every ``*.get_*()`` "results" accessor expects."""
    net = build_exotic_elements_network()
    grid = init_from_pypowsybl(net, gen_slack_id="B1-G")
    if solve:
        import numpy as np
        v0 = np.ones(grid.total_bus(), dtype=complex)
        conv = grid.ac_pf(v0, 20, 1e-7)
        if len(conv) == 0:
            raise RuntimeError("build_exotic_elements_grid: the ac_pf did not converge; "
                               "this fixture is expected to always converge, this is a bug.")
    return grid
