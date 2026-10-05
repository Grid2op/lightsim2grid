# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The NROuter_* algorithms (OpenLoadFlow's outer loops around a single-slack Newton), seen
from python. The driver itself is tested in C++ (src/tests/test_outer_loop_driver.cpp) with
scripted loops; each loop has its own test file."""

import copy
import os
import sys
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lightsim2grid.network import init_from_pandapower
from lightsim2grid.algorithm import ErrorType, OuterLoopStatus

try:
    from lightsim2grid.algorithm import NROuter_KLU  # noqa: F401
    OUTER, SING = "NROuter_KLU", "NRSing_KLU"
except ImportError:
    OUTER, SING = "NROuter_SparseLU", "NRSing_SparseLU"


def _grid(algo):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(pn.case118())
    grid.change_algorithm(algo)
    return grid


def _solve(grid):
    return grid.ac_pf(np.full(grid.total_bus(), 1.04 + 0j), 30, 1e-8)


class TestNROuter(unittest.TestCase):
    def test_without_loop_it_is_the_single_slack_newton(self):
        V_ref = _solve(_grid(SING))
        grid = _grid(OUTER)
        grid.clear_outer_loops()
        V = _solve(grid)
        self.assertEqual(V.shape, V_ref.shape)
        self.assertTrue(np.array_equal(V, V_ref))
        stats = grid.get_algo().get_outer_loop_stats()
        self.assertEqual(stats.status, OuterLoopStatus.STABLE)
        self.assertEqual(stats.nb_outer_iterations, 0)
        self.assertEqual(len(stats.nr_iterations), 1)
        self.assertTrue(grid.get_algo().supports_outer_loops())

    def test_one_analysis_over_several_solves(self):
        grid = _grid(OUTER)
        for _ in range(3):
            self.assertGreater(_solve(grid).shape[0], 0)
        self.assertEqual(grid.get_algo().get_linear_solver_stats().nb_analyze, 1)

    def test_loop_list(self):
        grid = _grid(OUTER)
        default = grid.get_outer_loops()
        grid.clear_outer_loops()
        self.assertEqual(grid.get_outer_loops(), [])
        grid.reset_outer_loops()
        self.assertEqual([l.name() for l in grid.get_outer_loops()], [l.name() for l in default])
        with self.assertRaises(RuntimeError):
            grid.add_outer_loop(None)

    def test_driver_parameters_round_trip(self):
        grid = _grid(OUTER)
        cfg = grid.get_ac_algo_config()
        # the Newton's 4 + 6 parameters, then the driver's 2 + 3
        self.assertEqual(len(cfg.int_params), 6)
        self.assertEqual(len(cfg.real_params), 9)
        self.assertEqual(cfg.int_params[4], 30)
        self.assertAlmostEqual(cfg.real_params[6], 0.8)
        self.assertAlmostEqual(cfg.real_params[7], 1.2)
        self.assertAlmostEqual(cfg.real_params[8], 180.)
        ip = list(cfg.int_params)
        ip[4] = 12
        cfg.int_params = ip
        grid.set_ac_algo_config(cfg)
        self.assertEqual(grid.get_ac_algo_config().int_params[4], 12)
        # a copy keeps the algorithm and its configuration
        grid2 = grid.copy()
        self.assertTrue(grid2.get_algo().supports_outer_loops())
        self.assertEqual(grid2.get_ac_algo_config().int_params[4], 12)
        self.assertEqual(copy.deepcopy(grid).get_ac_algo_config().int_params[4], 12)

    def test_unrealistic_state(self):
        grid = _grid(OUTER)
        cfg = grid.get_ac_algo_config()
        rp = list(cfg.real_params)
        rp[6] = 1.2  # nothing is realistic any more ...
        rp[7] = 1.3
        rp[8] = 0.   # ... at any nominal voltage
        cfg.real_params = rp
        grid.set_ac_algo_config(cfg)
        self.assertEqual(_solve(grid).shape[0], 0)
        self.assertEqual(grid.get_algo().get_error(), ErrorType.UnrealisticState)
        self.assertTrue(grid.get_algo().get_outer_loop_stats().unrealistic_state)

    def test_standalone_algorithm(self):
        from lightsim2grid import algorithm
        algo = getattr(algorithm, OUTER)()
        algo.max_outer_iterations = 7
        self.assertEqual(algo.max_outer_iterations, 7)
        with self.assertRaises(RuntimeError):
            algo.min_realistic_voltage = 2.  # above the max
        self.assertEqual(algo.get_outer_loop_stats().nb_outer_iterations, 0)


class TestDistributedSlack(unittest.TestCase):
    """OpenLoadFlow's DistributedSlack loop, from python, on pypowsybl's IEEE 14-bus grid built
    with OpenLoadFlow's loading rules (which flag the units taking part in the slack)."""

    @staticmethod
    def _net():
        import pypowsybl as pp
        net = pp.network.create_ieee14()
        gen = net.get_generators()
        net.update_generators(id=list(gen.index), min_p=[0.] * len(gen), max_p=[300.] * len(gen),
                              target_p=list(gen["target_p"].where(gen["target_p"] > 0., 10.)))
        return net

    def _grid(self, **kwargs):
        try:
            from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl
        except ImportError:
            self.skipTest("pypowsybl is not installed")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            return init_from_pypowsybl(self._net(), olf_rules=True, **kwargs)

    def test_in_the_default_list(self):
        from lightsim2grid.algorithm import DistributedSlack
        grid = self._grid(gen_slack_id="B1-G")
        self.assertEqual([loop.name() for loop in grid.get_outer_loops()],
                         ["DistributedSlack", "AcHvdcAcEmulationLimits", "VoltageMonitoring", "ReactiveLimits"])
        loop = DistributedSlack(slack_bus_p_max_mismatch_mw=2., fail_on_residue=False,
                                p_residue_eps_mw=1e-2, moved_fraction=0.5)
        self.assertEqual(loop.slack_bus_p_max_mismatch_mw, 2.)
        self.assertFalse(loop.fail_on_residue)
        self.assertEqual(loop.p_residue_eps_mw, 1e-2)
        self.assertEqual(loop.moved_fraction, 0.5)
        # set at construction only
        with self.assertRaises(AttributeError):
            loop.slack_bus_p_max_mismatch_mw = 3.
        with self.assertRaises(RuntimeError):
            DistributedSlack(p_residue_eps_mw=-1.)

    def test_same_as_the_newton_distributed_slack(self):
        from lightsim2grid.algorithm import DistributedSlack
        gen_ids = list(self._net().get_generators().index)
        # every unit has the same key here (max_p / default droop), as the in-Newton slack's
        ref = self._grid(gen_slack_id={gen_id: 1. for gen_id in gen_ids})
        ref.change_algorithm(SING.replace("NRSing", "NR"))
        V_ref = _solve(ref)
        grid = self._grid(gen_slack_id="B1-G")
        grid.change_algorithm(OUTER)
        grid.clear_outer_loops()
        grid.add_outer_loop(DistributedSlack(slack_bus_p_max_mismatch_mw=1e-6))
        V = _solve(grid)
        self.assertEqual(V.shape, V_ref.shape)
        stats = grid.get_algo().get_outer_loop_stats()
        self.assertEqual(stats.status, OuterLoopStatus.STABLE)
        self.assertGreater(stats.nb_outer_iterations, 0)
        self.assertEqual(grid.get_algo().get_linear_solver_stats().nb_analyze, 1)
        # the loop stops at OpenLoadFlow's residue (1e-3 MW)
        self.assertLess(np.abs(V - V_ref).max(), 1e-5)
        p = np.array([gen.res_p_mw for gen in grid.get_generators()])
        p_ref = np.array([gen.res_p_mw for gen in ref.get_generators()])
        self.assertLess(np.abs(p - p_ref).max(), 2e-3)

    def test_target_q_follows_the_moved_p(self):
        # a unit that does not regulate asks more than it can absorb: OpenLoadFlow clamps its
        # target Q into its limits at its CURRENT target P (forceTargetQInReactiveLimits),
        # which the slack moves
        import pypowsybl.loadflow as lf
        from lightsim2grid.algorithm import DistributedSlack
        from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl, OlfLoadingParameters
        from _olf_reference import reference_parameters

        def net_():
            net = self._net()
            net.update_generators(id="B3-G", voltage_regulator_on=False, target_q=-80.)
            net.create_curve_reactive_limits(id=["B3-G", "B3-G"], p=[0., 300.], min_q=[-10., -70.],
                                             max_q=[40., 40.])
            net.update_loads(id="B3-L", p0=net.get_loads().loc["B3-L", "p0"] + 150.)
            return net

        net = net_()
        params = reference_parameters(
            provider={"slackBusSelectionMode": "NAME", "slackBusesIds": net.get_generators().loc["B1-G", "bus_id"],
                      "outerLoopNames": "DistributedSlack", "newtonRaphsonConvEpsPerEq": "1e-10",
                      "slackBusPMaxMismatch": "1e-4"},
            distributed_slack=True, use_reactive_limits=True, read_slack_bus=False,
            balance_type=lf.BalanceType.PROPORTIONAL_TO_GENERATION_P_MAX)
        self.assertEqual(lf.run_ac(net, params)[0].status.name, "CONVERGED")
        olf = net.get_generators().loc["B3-G"]

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pypowsybl(net_(), gen_slack_id="B1-G", sort_index=False, buses_for_sub=False,
                                       olf_rules=OlfLoadingParameters(reactive_limits=True))
        grid.change_algorithm(OUTER)
        grid.clear_outer_loops()
        grid.add_outer_loop(DistributedSlack(slack_bus_p_max_mismatch_mw=1e-4))
        V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        gen = {g.name: g for g in grid.get_generators()}["B3-G"]
        self.assertAlmostEqual(gen.res_p_mw, -olf["p"], places=3)
        self.assertAlmostEqual(gen.res_q_mvar, -olf["q"], places=3)
        # at its min Q at the new P, not the initial one
        self.assertAlmostEqual(gen.res_q_mvar, grid.gen_limits_at(gen.id, gen.res_p_mw)[0], places=6)
        self.assertNotAlmostEqual(gen.res_q_mvar, gen.target_q_mvar, places=1)
        vm_olf = (net.get_buses()["v_mag"] / net.get_voltage_levels()["nominal_v"].reindex(
            net.get_buses()["voltage_level_id"]).to_numpy()).to_numpy()
        self.assertLess(np.abs(np.abs(V) - vm_olf).max(), 1e-6)

    def test_detection(self):
        from lightsim2grid.lightsim2grid_cpp import LimitViolationType
        grid = self._grid(gen_slack_id="B1-G")
        grid.change_algorithm(SING)
        self.assertGreater(_solve(grid).shape[0], 0)
        found = [v for v in grid.get_physical_violations()
                 if v.violation_type == LimitViolationType.SLACK_MISMATCH]
        self.assertEqual(len(found), 1)
        slack_gen = grid.get_generators()[0]
        self.assertAlmostEqual(found[0].value, slack_gen.res_p_mw - slack_gen.target_p_mw, places=6)
        # a single slack set up on purpose, with no unit flagged to share it, is not reported
        grid.set_gen_can_participate_slack([False] * len(list(grid.get_generators())),
                                           np.zeros(len(list(grid.get_generators()))))
        grid.set_storage_can_participate_slack([], np.zeros(0))
        self.assertGreater(_solve(grid).shape[0], 0)
        self.assertFalse([v for v in grid.get_physical_violations()
                          if v.violation_type == LimitViolationType.SLACK_MISMATCH])


class TestHvdcAcEmulationLimits(unittest.TestCase):
    """OpenLoadFlow's AcHvdcAcEmulationLimits loop, from python: an angle-droop hvdc line on
    case14 whose linear flow (from bus 3 to bus 9) is above what its converters transmit."""

    LF = 0.011

    def _model(self, **kwargs):
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from _aux_make_hvdc import make_case14_hvdc
        params = dict(loss_factor1=self.LF, loss_factor2=self.LF, droop_enabled=True, p0=10.,
                      droop_mw_per_deg=5., pmax12=15., pmax21=300.)
        params.update(kwargs)
        _, model = make_case14_hvdc(3, 9, **params)
        return model

    def _outer(self, **kwargs):
        from lightsim2grid.algorithm import HvdcAcEmulationLimits
        model = self._model(**kwargs)
        model.change_algorithm(OUTER)
        model.clear_outer_loops()
        model.add_outer_loop(HvdcAcEmulationLimits())
        return model

    def test_saturates_the_line(self):
        from _aux_make_hvdc import recv_mw
        model = self._outer()
        V = _solve(model)
        self.assertGreater(V.shape[0], 0)
        stats = model.get_algo().get_outer_loop_stats()
        self.assertEqual(stats.status, OuterLoopStatus.STABLE)
        self.assertEqual(stats.loop_iterations, [("AcHvdcAcEmulationLimits", 1)])
        self.assertEqual(model.get_algo().get_linear_solver_stats().nb_analyze, 1)
        line = model.get_dclines()[0]
        self.assertAlmostEqual(line.res_p1_mw, -15., places=8)
        self.assertAlmostEqual(line.res_p2_mw, recv_mw(15., self.LF, self.LF), places=8)
        # the regime is the algorithm's: the line's own is untouched
        self.assertEqual(model.get_status_droop_hvdc(0), 0)

        # the same as a single-slack Newton with the line saturated by hand
        ref = self._model()
        ref.change_algorithm(SING)
        ref.set_status_droop_hvdc(0, 1)
        V_ref = _solve(ref)
        self.assertLess(np.abs(V - V_ref).max(), 1e-8)

    def test_saturates_in_the_reverse_direction(self):
        model = self._outer(p0=-60., pmax12=300., pmax21=20.)
        self.assertGreater(_solve(model).shape[0], 0)
        self.assertEqual(model.get_algo().get_outer_loop_stats().status, OuterLoopStatus.STABLE)
        self.assertAlmostEqual(model.get_dclines()[0].res_p2_mw, -20., places=8)

    def test_inside_its_limits(self):
        model = self._outer(pmax12=300.)
        V = _solve(model)
        stats = model.get_algo().get_outer_loop_stats()
        self.assertEqual(stats.nb_outer_iterations, 0)
        ref = self._model(pmax12=300.)
        ref.change_algorithm(SING)
        self.assertTrue(np.array_equal(V, _solve(ref)))

    def test_not_needed_without_ac_emulation(self):
        model = self._outer(droop_enabled=False)
        self.assertGreater(_solve(model).shape[0], 0)
        self.assertEqual(model.get_algo().get_outer_loop_stats().loop_iterations, [])

    def test_detection(self):
        from lightsim2grid.lightsim2grid_cpp import LimitViolationType
        model = self._model()
        model.change_algorithm(SING)
        self.assertGreater(_solve(model).shape[0], 0)
        found = [v for v in model.get_physical_violations()
                 if v.violation_type == LimitViolationType.HIGH_P]
        self.assertEqual(len(found), 1)
        self.assertGreater(found[0].value, 15.)


class TestVoltageMonitoring(unittest.TestCase):
    """OpenLoadFlow's VoltageMonitoring loop, from python, against OpenLoadFlow itself: a
    standby SVC on load bus 9 of pypowsybl's IEEE 14-bus grid."""

    def _net(self, low_pu, high_pu, b0=0.):
        try:
            import pypowsybl as pp
        except ImportError:
            self.skipTest("pypowsybl is not installed")
        net = pp.network.create_ieee14()
        vn = float(net.get_voltage_levels().loc["VL9", "nominal_v"])
        net.create_static_var_compensators(id="SVC9", voltage_level_id="VL9", bus_id="B9",
                                           connectable_bus_id="B9", b_min=-0.5, b_max=0.5,
                                           regulation_mode="VOLTAGE", target_v=vn, target_q=0.,
                                           regulating=True)
        net.create_extensions("standbyAutomaton", id="SVC9", b0=b0, standby=True,
                              low_voltage_threshold=low_pu * vn, high_voltage_threshold=high_pu * vn,
                              low_voltage_setpoint=0.98 * vn, high_voltage_setpoint=1.02 * vn)
        return net

    @staticmethod
    def _olf(net):
        import pypowsybl.loadflow as lf
        from lightsim2grid.network.from_pypowsybl._olf_compare import iidm_bus_voltages
        from _olf_reference import reference_parameters
        params = reference_parameters(
            provider={"slackBusSelectionMode": "NAME",
                      "slackBusesIds": net.get_generators().loc["B1-G", "bus_id"],
                      "outerLoopNames": "VoltageMonitoring", "svcVoltageMonitoring": "true",
                      "newtonRaphsonConvEpsPerEq": "1e-12", "maxNewtonRaphsonIterations": "50"},
            distributed_slack=False, use_reactive_limits=False, read_slack_bus=False)
        res = lf.run_ac(net, params)[0]
        return res, iidm_bus_voltages(net)["vm_pu"], net.get_static_var_compensators().loc["SVC9", "q"]

    def _ls(self, net):
        from lightsim2grid.algorithm import VoltageMonitoring
        from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl, OlfLoadingParameters
        from lightsim2grid.network.from_pypowsybl._olf_compare import lightsim_bus_to_iidm
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pypowsybl(net, gen_slack_id="B1-G", sort_index=False, buses_for_sub=False,
                                       olf_rules=OlfLoadingParameters(reactive_limits=False))
        grid.change_algorithm(OUTER)
        grid.clear_outer_loops()
        grid.add_outer_loop(VoltageMonitoring())
        V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-12)
        self.assertGreater(V.shape[0], 0)
        self.assertEqual(grid.get_algo().get_linear_solver_stats().nb_analyze, 1)
        to_iidm = lightsim_bus_to_iidm(grid, net)
        vm = {to_iidm[i]: abs(V[i]) for i in range(V.shape[0]) if i in to_iidm}
        svc = next(svc for svc in grid.get_svcs() if svc.name == "SVC9")
        return grid, vm, svc

    def _same_as_olf(self, low_pu, high_pu, b0=0.):
        net = self._net(low_pu, high_pu, b0)
        grid, vm, svc = self._ls(net)
        res, olf_vm, olf_q = self._olf(self._net(low_pu, high_pu, b0))
        self.assertEqual(res.status.name, "CONVERGED")
        self.assertLess(max(abs(vm[b] - olf_vm[b]) for b in vm if b in olf_vm.index), 1e-9)
        # OpenLoadFlow reports the SVC's Q, b0 included, in the load convention
        self.assertAlmostEqual(svc.res_q_mvar, -olf_q, places=5)
        return grid, vm, svc

    def test_switched_on_above_the_high_threshold(self):
        grid, vm, svc = self._same_as_olf(0.90, 1.00)
        self.assertAlmostEqual(vm["VL9_0"], 1.02, places=9)
        stats = grid.get_algo().get_outer_loop_stats()
        self.assertEqual(stats.loop_iterations, [("VoltageMonitoring", 1)])
        # the SVC itself is untouched
        self.assertEqual(svc.regulation_mode, 0)
        self.assertTrue(svc.standby)

    def test_switched_on_below_the_low_threshold(self):
        _, vm, _ = self._same_as_olf(1.10, 1.20)
        self.assertAlmostEqual(vm["VL9_0"], 0.98, places=9)

    def test_idle_inside_the_thresholds(self):
        grid, _, svc = self._same_as_olf(0.90, 1.20)
        self.assertEqual(grid.get_algo().get_outer_loop_stats().nb_outer_iterations, 0)
        self.assertAlmostEqual(svc.res_q_mvar, 0., places=9)

    def test_b0(self):
        # a fixed susceptance carried by the SVC, idle and switched on
        _, _, svc = self._same_as_olf(0.90, 1.20, b0=0.05)
        self.assertGreater(svc.b0_pu, 0.)
        self.assertGreater(svc.res_q_mvar, 0.)
        self._same_as_olf(0.90, 1.00, b0=0.05)

    def test_in_the_default_list(self):
        grid, _, _ = self._ls(self._net(0.90, 1.20))
        grid.reset_outer_loops()
        self.assertEqual([loop.name() for loop in grid.get_outer_loops()],
                         ["DistributedSlack", "AcHvdcAcEmulationLimits", "VoltageMonitoring", "ReactiveLimits"])


class TestReactiveLimits(unittest.TestCase):
    """OpenLoadFlow's ReactiveLimits loop, from python, against OpenLoadFlow itself on
    pypowsybl's IEEE 300-bus grid (every unit regulates its own bus)."""

    def _net(self):
        try:
            import pypowsybl as pp
        except ImportError:
            self.skipTest("pypowsybl is not installed")
        return pp.network.create_ieee300()

    def test_same_as_olf(self):
        import pypowsybl.loadflow as lf
        from lightsim2grid.algorithm import ReactiveLimits
        from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl, OlfLoadingParameters
        from lightsim2grid.network.from_pypowsybl._olf_compare import iidm_bus_voltages
        from lightsim2grid.network.from_pypowsybl._aux_add_slack import _default_distributed_slack
        net = self._net()
        gen = net.get_generators()
        # the unit the comparison harness (utils/olf_outer_compare.py) takes as the slack
        slack_gen = next(iter(_default_distributed_slack(net, gen)))
        from _olf_reference import reference_parameters
        params = reference_parameters(
            provider={"slackBusSelectionMode": "NAME", "slackBusesIds": gen.loc[slack_gen, "bus_id"],
                      "outerLoopNames": "ReactiveLimits",
                      "newtonRaphsonConvEpsPerEq": "1e-9", "maxNewtonRaphsonIterations": "50"},
            distributed_slack=False, use_reactive_limits=True, read_slack_bus=False)
        res = lf.run_ac(net, params)[0]
        self.assertEqual(res.status.name, "CONVERGED")
        olf_vm = iidm_bus_voltages(net)["vm_pu"]
        olf_q = net.get_generators()["q"]

        net = self._net()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pypowsybl(net, gen_slack_id=slack_gen, sort_index=False, buses_for_sub=False,
                                       keep_half_open_lines=True, fuse_zero_impedance_branches=True,
                                       olf_rules=OlfLoadingParameters(reactive_limits=True))
        grid.change_algorithm(OUTER)
        grid.clear_outer_loops()
        grid.add_outer_loop(ReactiveLimits(max_reactive_power_mismatch=1e-9))
        # OpenLoadFlow's DC_VALUES start: the angles of a DC load flow, 1 pu
        flat = np.full(grid.total_bus(), 1.0 + 0j)
        V0 = np.exp(1j * np.angle(grid.dc_pf(flat, 30, 1e-12)))
        V = grid.ac_pf(V0, 30, 1e-9)
        self.assertGreater(V.shape[0], 0)
        stats = grid.get_algo().get_outer_loop_stats()
        self.assertGreater(stats.nb_outer_iterations, 0)
        self.assertEqual(grid.get_algo().get_linear_solver_stats().nb_analyze, 1)
        # through the result view: a bus merged by fuse_zero_impedance_branches reads the
        # voltage of the bus it was merged into
        from lightsim2grid.network.from_pypowsybl import LightsimResultNetwork
        ls_bus = LightsimResultNetwork(grid, net).get_buses()
        nominal = net.get_voltage_levels()["nominal_v"]
        ls_vm = (ls_bus["v_mag"] / ls_bus["voltage_level_id"].map(nominal)).dropna()
        common = ls_vm.index.intersection(olf_vm.index)
        self.assertGreater(len(common), 290)
        self.assertLess((ls_vm[common] - olf_vm[common]).abs().max(), 1e-8)
        # the units frozen at a limit publish it
        compared = 0
        for g in grid.get_generators():
            if g.name in olf_q.index and g.connected and np.isfinite(olf_q[g.name]):
                self.assertAlmostEqual(g.res_q_mvar, -olf_q[g.name], places=4)
                compared += 1
        self.assertGreater(compared, 50)

    def test_detection_is_olf_first_round(self):
        # a plain solve reports exactly the buses OpenLoadFlow switches PV -> PQ first
        import re
        import pypowsybl.loadflow as lf
        import pypowsybl.report as rp
        from lightsim2grid.lightsim2grid_cpp import LimitViolationType
        from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl, OlfLoadingParameters
        from lightsim2grid.network.from_pypowsybl._aux_add_slack import _default_distributed_slack
        from lightsim2grid.network.from_pypowsybl._olf_compare import lightsim_bus_to_iidm
        from _olf_reference import reference_parameters
        net = self._net()
        slack_gen = next(iter(_default_distributed_slack(net, net.get_generators())))
        params = reference_parameters(
            provider={"slackBusSelectionMode": "NAME", "slackBusesIds": net.get_generators().loc[slack_gen, "bus_id"],
                      "outerLoopNames": "ReactiveLimits", "newtonRaphsonConvEpsPerEq": "1e-9"},
            distributed_slack=False, use_reactive_limits=True, read_slack_bus=False)
        report = rp.ReportNode()
        lf.run_ac(net, params, report_node=report)
        first = re.split(r"\n\s*\+ (?=\d+ buses switched PV -> PQ)", str(report))
        self.assertGreater(len(first), 1)
        olf = set(re.findall(r"Switch bus '([^']+)' PV -> PQ", first[1]))

        net = self._net()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pypowsybl(net, gen_slack_id=slack_gen, sort_index=False, buses_for_sub=False,
                                       keep_half_open_lines=True, fuse_zero_impedance_branches=True,
                                       olf_rules=OlfLoadingParameters(reactive_limits=True))
        grid.change_algorithm(SING)
        flat = np.full(grid.total_bus(), 1.0 + 0j)
        V0 = np.exp(1j * np.angle(grid.dc_pf(flat, 30, 1e-12)))
        self.assertGreater(grid.ac_pf(V0, 30, 1e-9).shape[0], 0)
        to_iidm = lightsim_bus_to_iidm(grid, net)
        # OpenLoadFlow's maxReactivePowerMismatch (1e-9 pu of 100 MVA)
        ours = {to_iidm[v.element_id] for v in grid.get_physical_violations(True, 1e-7, 0.)
                if v.violation_type in (LimitViolationType.LOW_Q, LimitViolationType.HIGH_Q)}
        self.assertGreater(len(olf), 0)
        self.assertEqual(ours, olf)

    def test_curve_limits(self):
        # a generator's limits read on its capability curve, as OpenLoadFlow reads them
        grid = _grid(SING)
        grid.set_gen_capability_curves(np.array([0, 0, 0], dtype=np.int32), np.array([0., 50., 100.]),
                                       np.array([-10., -20., -10.]), np.array([10., 30., 10.]))
        self.assertEqual(grid.gen_limits_at(0, 25.), (-15., 20.))     # interpolated
        self.assertEqual(grid.gen_limits_at(0, 120.), (-6., 2.))      # extrapolated
        self.assertEqual(grid.gen_limits_at(0, 150.), (-5., -5.))     # crossed: their mean
        self.assertEqual(grid.gen_limits_at(1, 25.), (None, None))    # no curve
        self.assertEqual(grid.copy().gen_limits_at(0, 25.), (-15., 20.))

    def test_parameters(self):
        from lightsim2grid.algorithm import ReactiveLimits
        loop = ReactiveLimits(max_pq_pv_switch=5, max_reactive_power_mismatch=1e-6, mismatch_base_mva=10.,
                              realistic_voltage_margin=1.05, robust_restart_vm_pu=0.95)
        self.assertEqual(loop.max_pq_pv_switch, 5)
        self.assertAlmostEqual(loop.max_reactive_power_mismatch, 1e-6)
        self.assertEqual(loop.mismatch_base_mva, 10.)
        self.assertEqual(loop.realistic_voltage_margin, 1.05)
        self.assertEqual(loop.robust_restart_vm_pu, 0.95)
        # set at construction only
        with self.assertRaises(AttributeError):
            loop.max_pq_pv_switch = 2
        with self.assertRaises(RuntimeError):
            ReactiveLimits(max_pq_pv_switch=-1)
        with self.assertRaises(RuntimeError):
            ReactiveLimits(mismatch_base_mva=0.)


class TestReactiveDispatch(unittest.TestCase):
    """LSGrid.set_reactive_dispatch_olf: the reactive power of a bus held by several units is
    split between them as OpenLoadFlow does, against OpenLoadFlow itself (no outer loop)."""

    def _net(self):
        try:
            import pypowsybl as pp
        except ImportError:
            self.skipTest("pypowsybl is not installed")
        net = pp.network.create_ieee14()
        gen = net.get_generators()
        bus = net.get_bus_breaker_view_buses().query("voltage_level_id == 'VL2'").index[0]
        # a high set-point: the bus needs a lot of reactive power
        target_v = 1.06 * net.get_voltage_levels().loc["VL2", "nominal_v"]
        net.update_generators(id="B2-G", target_v=target_v)
        # three units on bus 2, one with a narrow range it has to stop at
        for gen_id in ("B2-G2", "B2-G3"):
            net.create_generators(id=gen_id, voltage_level_id="VL2", bus_id=bus, target_p=10., min_p=0.,
                                  max_p=100., target_v=target_v, voltage_regulator_on=True)
        # the narrow unit's range is lopsided: its share of the bus exceeds what it can produce
        net.create_minmax_reactive_limits(id=["B2-G", "B2-G2", "B2-G3"], min_q=[-40., -10., -30.],
                                          max_q=[50., 2., 60.])
        return net

    def test_same_split_as_olf(self):
        import pypowsybl.loadflow as lf
        from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl, OlfLoadingParameters
        net = self._net()
        from _olf_reference import reference_parameters
        params = reference_parameters(
            provider={"slackBusSelectionMode": "NAME", "slackBusesIds": net.get_generators().loc["B1-G", "bus_id"],
                      "outerLoopNames": "", "newtonRaphsonConvEpsPerEq": "1e-10"},
            distributed_slack=False, use_reactive_limits=True, read_slack_bus=False)
        res = lf.run_ac(net, params)[0]
        self.assertEqual(res.status.name, "CONVERGED")
        olf_q = net.get_generators()["q"]

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pypowsybl(self._net(), gen_slack_id="B1-G", sort_index=False, buses_for_sub=False,
                                       olf_rules=OlfLoadingParameters(reactive_limits=True))
        self.assertTrue(grid.get_reactive_dispatch_olf())
        grid.change_algorithm(SING)
        self.assertGreater(grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10).shape[0], 0)
        ls_q = {g.name: g.res_q_mvar for g in grid.get_generators()}
        for gen_id in ("B2-G", "B2-G2", "B2-G3"):
            self.assertAlmostEqual(ls_q[gen_id], -olf_q[gen_id], places=5)
        # the narrow unit stopped at its limit
        self.assertAlmostEqual(abs(ls_q["B2-G2"]), 2., places=6)

        # off: the former split, by reactive range at the target P, no limit
        grid.set_reactive_dispatch_olf(False)
        self.assertGreater(grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10).shape[0], 0)
        self.assertGreater(abs({g.name: g.res_q_mvar for g in grid.get_generators()}["B2-G2"]), 2.)


if __name__ == "__main__":
    unittest.main()
