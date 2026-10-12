# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The outer loops (and the driver's parameters) OpenLoadFlow would run with some pypowsybl
loadflow parameters: which loops, in which order, with which parameters -- by default those
of the installed pypowsybl, which a grid loaded from it runs."""

import unittest
import warnings

from lightsim2grid.lightsim2grid_cpp import LSGrid
from lightsim2grid.algorithm import OuterLoopDriverParams, ReactiveLimits, TransformerVoltageControl

try:
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.network import init_from_pypowsybl, olf_outer_loops, set_olf_outer_loops
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False


def _names(loops):
    return [loop.name() for loop in loops]


def _quiet(parameters=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return olf_outer_loops(parameters)


class TestDefaultsWithoutPypowsybl(unittest.TestCase):
    """without pypowsybl, OpenLoadFlow's own defaults (OpenLoadFlowParameters,
    AcLoadFlowParameters, AbstractLoadFlowParameters, LfNetworkParameters at v2.3.0)"""

    def test_driver(self):
        self.assertEqual(LSGrid().get_outer_loop_driver_params(),
                         OuterLoopDriverParams(max_outer_iterations=20, voltage_remote_control_robust_mode=True,
                                               min_realistic_voltage=0.5, max_realistic_voltage=2.,
                                               min_nominal_voltage_realistic_check=0.))

    def test_loops(self):
        grid = LSGrid()
        # DefaultAcOuterLoopConfig's order, the loops its default parameters create
        self.assertEqual(_names(grid.get_outer_loops()),
                         ["DistributedSlack", "AcHvdcAcEmulationLimits", "VoltageMonitoring", "ReactiveLimits"])
        loop = ReactiveLimits()
        self.assertEqual((loop.min_realistic_voltage, loop.max_realistic_voltage), (0.5, 2.))
        self.assertEqual(loop.max_reactive_power_mismatch, 1e-4)
        self.assertEqual(loop.max_pq_pv_switch, 3)
        loop = TransformerVoltageControl()
        self.assertFalse(loop.use_initial_tap_position)
        self.assertEqual(loop.max_controlled_nominal_voltage, -1.)


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestOlfOuterLoops(unittest.TestCase):
    def test_installed_defaults(self):
        # whatever the build: the defaults it ships
        defaults = lf.get_provider_parameters("OpenLoadFlow")["default"]
        loops, driver = _quiet()
        self.assertEqual(driver.max_outer_iterations, int(defaults["maxOuterLoopIterations"]))
        self.assertEqual(driver.min_realistic_voltage, float(defaults["minRealisticVoltage"]))
        self.assertEqual(driver.max_realistic_voltage, float(defaults["maxRealisticVoltage"]))
        self.assertEqual(driver.min_nominal_voltage_realistic_check,
                         float(defaults["minNominalVoltageRealisticVoltageCheck"]))
        self.assertEqual(driver.voltage_remote_control_robust_mode,
                         defaults["voltageRemoteControlRobustMode"].lower() == "true")
        params = lf.Parameters()
        expected = (["DistributedSlack"] if params.distributed_slack else []) + \
                   (["AcHvdcAcEmulationLimits"] if params.hvdc_ac_emulation else []) + \
                   (["VoltageMonitoring"] if defaults["svcVoltageMonitoring"].lower() == "true" else []) + \
                   (["ReactiveLimits"] if params.use_reactive_limits else [])
        self.assertEqual(_names(loops)[:len(expected)], expected)
        reactive = next(loop for loop in loops if loop.name() == "ReactiveLimits")
        self.assertEqual(reactive.min_realistic_voltage, float(defaults["minRealisticVoltage"]))
        self.assertEqual(reactive.max_pq_pv_switch, int(defaults["reactiveLimitsMaxPqPvSwitch"]))

    def test_parameters_are_read(self):
        params = lf.Parameters(distributed_slack=True, use_reactive_limits=True, transformer_voltage_control_on=True,
                               shunt_compensator_voltage_control_on=True, phase_shifter_regulation_on=True)
        params.provider_parameters = {
            "maxOuterLoopIterations": "13", "voltageRemoteControlRobustMode": "false",
            "minRealisticVoltage": "0.7", "maxRealisticVoltage": "1.4", "minNominalVoltageRealisticVoltageCheck": "90",
            "reactiveLimitsMaxPqPvSwitch": "5", "slackBusPMaxMismatch": "0.25",
            "slackDistributionFailureBehavior": "LEAVE_ON_SLACK_BUS",
            "newtonRaphsonStoppingCriteriaType": "PER_EQUATION_TYPE_CRITERIA", "maxReactivePowerMismatch": "0.5",
            "transformerVoltageControlMode": "AFTER_GENERATOR_VOLTAGE_CONTROL",
            "transformerVoltageControlUseInitialTapPosition": "true", "generatorVoltageControlMinNominalVoltage": "63",
            "shuntVoltageControlMode": "WITH_GENERATOR_VOLTAGE_CONTROL",
            "phaseShifterControlMode": "CONTINUOUS_WITH_DISCRETISATION",
        }
        loops, driver = _quiet(params)
        self.assertEqual(driver, OuterLoopDriverParams(max_outer_iterations=13, voltage_remote_control_robust_mode=False,
                                                       min_realistic_voltage=0.7, max_realistic_voltage=1.4,
                                                       min_nominal_voltage_realistic_check=90.))
        # OpenLoadFlow's order, whatever order they were switched on in
        names = _names(loops)
        self.assertEqual(names[names.index("ReactiveLimits"):],
                         ["ReactiveLimits", "PhaseControl", "TransformerVoltageControl", "ShuntVoltageControl"])
        by_name = dict(zip(names, loops))
        self.assertEqual(by_name["DistributedSlack"].slack_bus_p_max_mismatch_mw, 0.25)
        self.assertFalse(by_name["DistributedSlack"].fail_on_residue)
        reactive = by_name["ReactiveLimits"]
        self.assertEqual(reactive.max_pq_pv_switch, 5)
        self.assertFalse(reactive.robust_mode)
        self.assertEqual((reactive.min_realistic_voltage, reactive.max_realistic_voltage), (0.7, 1.4))
        # per equation type: maxReactivePowerMismatch, MVar, in pu of 100 MVA
        self.assertAlmostEqual(reactive.max_reactive_power_mismatch, 0.5 / 100., places=15)
        self.assertTrue(by_name["TransformerVoltageControl"].use_initial_tap_position)
        self.assertEqual(by_name["TransformerVoltageControl"].max_controlled_nominal_voltage, 63.)

    def test_switched_off(self):
        params = lf.Parameters(distributed_slack=False, use_reactive_limits=False, hvdc_ac_emulation=False)
        params.provider_parameters = {"svcVoltageMonitoring": "false"}
        self.assertEqual(_quiet(params)[0], [])

    def test_explicit_names(self):
        params = lf.Parameters(distributed_slack=True, use_reactive_limits=True, transformer_voltage_control_on=True)
        # the list's own order; the name decides the mode, not transformerVoltageControlMode
        params.provider_parameters = {"outerLoopNames": "TransformerVoltageControl,ReactiveLimits,DistributedSlack",
                                      "transformerVoltageControlMode": "INCREMENTAL_VOLTAGE_CONTROL"}
        self.assertEqual(_names(_quiet(params)[0]), ["TransformerVoltageControl", "ReactiveLimits", "DistributedSlack"])
        # the name of later OpenLoadFlow versions
        params.provider_parameters = {"outerLoopNames": "HvdcAcEmulationLimits"}
        self.assertEqual(_names(_quiet(params)[0]), ["AcHvdcAcEmulationLimits"])
        # a loop its switch leaves off is not created, even when named
        params.use_reactive_limits = False
        params.provider_parameters = {"outerLoopNames": "ReactiveLimits,DistributedSlack"}
        self.assertEqual(_names(_quiet(params)[0]), ["DistributedSlack"])
        # given empty: no loop at all
        params.provider_parameters = {"outerLoopNames": ""}
        self.assertEqual(_quiet(params)[0], [])
        for bad in ("NoSuchLoop", "DistributedSlack,DistributedSlack"):
            params.provider_parameters = {"outerLoopNames": bad}
            with self.assertRaises(ValueError):
                olf_outer_loops(params)

    def test_missing_loops_warn(self):
        params = lf.Parameters(transformer_voltage_control_on=True)
        params.provider_parameters = {"transformerVoltageControlMode": "INCREMENTAL_VOLTAGE_CONTROL",
                                      "simulateAutomationSystems": "true"}
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loops, _ = olf_outer_loops(params)
        self.assertNotIn("TransformerVoltageControl", _names(loops))
        message = " ".join(str(w.message) for w in caught)
        self.assertIn("TransformerVoltageControl in mode INCREMENTAL_VOLTAGE_CONTROL", message)
        self.assertIn("AutomationSystem", message)

    def test_grid_loaded_from_pypowsybl(self):
        net = pp.network.create_ieee14()
        params = lf.Parameters(phase_shifter_regulation_on=True)
        params.provider_parameters = {"maxOuterLoopIterations": "7"}
        for given in (None, params):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                grid = init_from_pypowsybl(net, gen_slack_id="B1-G", olf_parameters=given)
            loops, driver = _quiet(given)
            self.assertEqual(_names(grid.get_outer_loops()), _names(loops))
            self.assertEqual(grid.get_outer_loop_driver_params(), driver)
        self.assertEqual(grid.get_outer_loop_driver_params().max_outer_iterations, 7)
        self.assertIn("PhaseControl", _names(grid.get_outer_loops()))
        # and on any grid
        other = LSGrid()
        set_olf_outer_loops(other, params)
        self.assertEqual(_names(other.get_outer_loops()), _names(grid.get_outer_loops()))
        self.assertEqual(other.get_outer_loop_driver_params().max_outer_iterations, 7)


if __name__ == "__main__":
    unittest.main()
