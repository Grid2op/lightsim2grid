# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""OpenLoadFlow's ShuntVoltageControl outer loop (WITH_GENERATOR_VOLTAGE_CONTROL) against
OpenLoadFlow itself: two switched shunts added to a bus of pypowsybl's IEEE 14-bus grid,
regulating its voltage."""

import unittest
import warnings

import numpy as np

try:
    import pandas as pd
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl, OlfLoadingParameters
    from lightsim2grid.network.from_pypowsybl._olf_compare import iidm_bus_voltages
    from lightsim2grid.network.from_pypowsybl._result_network import LightsimResultNetwork
    from lightsim2grid.algorithm import ShuntVoltageControl
    from lightsim2grid.lightsim2grid_cpp import LimitViolationType, ViolationCategory, ViolationElementType
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

SLACK_GEN = "B1-G"
SHUNTS = ["SH-A", "SH-B"]


def _net(bus, target_kv, sections=(3, 1), with_trafo_control=False):
    """IEEE 14 with two linear shunts on `bus` (10 sections of 0.02 S, 4 of 0.05 S), regulating
    it at `target_kv`; with `with_trafo_control`, T4-9-1 also regulates bus 9 at 12.6 kV."""
    net = pp.network.create_ieee14()
    vl = "VL" + bus[1:]
    sh = pd.DataFrame.from_records(
        index="id", columns=["id", "model_type", "section_count", "target_v", "target_deadband", "voltage_level_id",
                             "bus_id", "connectable_bus_id"],
        data=[("SH-A", "LINEAR", sections[0], target_kv, 0., vl, bus, bus),
              ("SH-B", "LINEAR", sections[1], target_kv, 0., vl, bus, bus)])
    model = pd.DataFrame.from_records(index="id", columns=["id", "g_per_section", "b_per_section", "max_section_count"],
                                      data=[("SH-A", 0., 0.02, 10), ("SH-B", 0., 0.05, 4)])
    net.create_shunt_compensators(sh, model)
    net.update_shunt_compensators(id=SHUNTS, voltage_regulation_on=[True, True])
    if with_trafo_control:
        rtc = pd.DataFrame.from_records(
            index="id", columns=["id", "target_deadband", "target_v", "oltc", "low_tap", "tap", "regulating", "regulated_side"],
            data=[("T4-9-1", 0.1, 12.6, True, 0, 10, False, "TWO")])
        steps = pd.DataFrame.from_records(index="id", columns=["id", "b", "g", "r", "x", "rho"],
                                          data=[("T4-9-1", 0., 0., 0., 0., 0.9 + 0.01 * k) for k in range(21)])
        net.create_ratio_tap_changers(rtc, steps)
        net.update_ratio_tap_changers(id="T4-9-1", regulating=True)
    return net


def _olf(net, with_trafo_control):
    from _olf_reference import reference_parameters
    loops = ("TransformerVoltageControl," if with_trafo_control else "") + "ShuntVoltageControl"
    params = reference_parameters(
        provider={"slackBusSelectionMode": "NAME", "slackBusesIds": net.get_generators().loc[SLACK_GEN, "bus_id"],
                  "newtonRaphsonConvEpsPerEq": "1e-10", "outerLoopNames": loops},
        distributed_slack=False, use_reactive_limits=False, read_slack_bus=False,
        shunt_compensator_voltage_control_on=True, transformer_voltage_control_on=with_trafo_control)
    res = lf.run_ac(net, params)[0]
    assert res.status.name == "CONVERGED", res.status_text
    counts = [int(c) for c in net.get_shunt_compensators(all_attributes=True).loc[SHUNTS, "solved_section_count"]]
    return counts, params


def _grid(net, olf_parameters=None):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pypowsybl(net, gen_slack_id=SLACK_GEN, sort_index=False, buses_for_sub=False,
                                   olf_rules=OlfLoadingParameters(reactive_limits=False),
                                   olf_parameters=olf_parameters)


class TestShuntVoltageControl(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is not installed")

    def _compare(self, bus, target_kv, sections=(3, 1), with_trafo_control=False):
        net = _net(bus, target_kv, sections, with_trafo_control)
        olf_counts, params = _olf(net, with_trafo_control)

        # the loops OpenLoadFlow ran, read off the same parameters
        grid = _grid(_net(bus, target_kv, sections, with_trafo_control), olf_parameters=params)
        self.assertEqual([l.name() for l in grid.get_outer_loops()],
                         (["TransformerVoltageControl"] if with_trafo_control else []) + ["ShuntVoltageControl"])
        grid.change_algorithm("NROuter_KLU")
        V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        self.assertEqual(grid.get_algo().get_outer_loop_stats().status.name, "STABLE")
        self.assertEqual(grid.get_algo().get_linear_solver_stats().nb_analyze, 1)

        res = {s.name: s for s in grid.get_shunts()}
        self.assertEqual([res[s].res_section_count for s in SHUNTS], olf_counts)
        self.assertEqual([res[s].section_count for s in SHUNTS], list(sections))  # inputs untouched
        olf_v = iidm_bus_voltages(net)
        ls = LightsimResultNetwork(grid, net).get_buses()
        vm = ls["v_mag"] / ls["voltage_level_id"].map(net.get_voltage_levels()["nominal_v"])
        common = olf_v.index.intersection(vm.dropna().index)
        self.assertLess((olf_v.loc[common, "vm_pu"] - vm.loc[common]).abs().max(), 1e-9)
        return olf_counts

    def test_sections(self):
        # the susceptance shared out, the largest shunt first
        self.assertEqual(self._compare("B9", 12.8), [2, 0])
        self.assertEqual(self._compare("B9", 13.2), [6, 3])
        self.assertEqual(self._compare("B9", 12.4), [0, 0])
        self.assertEqual(self._compare("B10", 12.8, sections=(5, 0)), [2, 1])

    def test_hidden_by_a_generator(self):
        # bus 6 is a generator's: the shunts keep their sections, not even reshuffled
        self.assertEqual(self._compare("B6", 12.8, sections=(0, 2)), [0, 2])

    def test_hidden_by_a_transformer(self):
        # a transformer regulating the same bus takes precedence
        self.assertEqual(self._compare("B9", 13.2, sections=(0, 2), with_trafo_control=True), [0, 2])

    def test_detection(self):
        grid = _grid(_net("B9", 13.2))
        grid.change_algorithm("NRSing_KLU")
        grid.clear_outer_loops()
        grid.add_outer_loop(ShuntVoltageControl())
        self.assertGreater(grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10).shape[0], 0)
        viols = [v for v in grid.get_physical_violations(True, 1e-3, 1e-4)
                 if v.violation_type == LimitViolationType.SHUNT_VOLTAGE_CONTROL]
        self.assertEqual(len(viols), 1)
        self.assertEqual(viols[0].element_type, ViolationElementType.BUS)
        self.assertEqual(viols[0].category, ViolationCategory.CONTROL)
        self.assertAlmostEqual(viols[0].limit, 13.2, places=9)

    def test_not_in_the_default_list(self):
        grid = _grid(_net("B9", 13.2))
        self.assertNotIn("ShuntVoltageControl", [l.name() for l in grid.get_outer_loops()])

    def test_a_contingency_trips_a_regulating_shunt(self):
        from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP
        import copy
        grid = _grid(_net("B9", 13.2))
        grid.change_algorithm("NROuter_KLU")
        shunt_id = [s.name for s in grid.get_shunts()].index("SH-A")
        v0 = np.full(grid.total_bus(), 1.0 + 0j)

        # no ShuntVoltageControl loop (OpenLoadFlow's default): nothing owns the shunt's
        # sections, a contingency may take it out
        self.assertNotIn("ShuntVoltageControl", [l.name() for l in grid.get_outer_loops()])
        ca = ContingencyAnalysisCPP(grid)
        self.assertIn("NROuter", ca.get_algo_name())
        ca.add_contingency_shunt(shunt_id)
        ca.compute(v0, 30, 1e-10)
        ref = copy.deepcopy(grid)
        ref.deactivate_shunt(shunt_id)
        V_ref = ref.ac_pf(v0.copy(), 30, 1e-10)
        self.assertGreater(V_ref.shape[0], 0)
        np.testing.assert_allclose(ca.get_voltages()[0], V_ref, atol=1e-8, rtol=0.)

        # the loop in the list: not supported by the batch, it owns the shunt's sections
        grid.add_outer_loop(ShuntVoltageControl())
        ca = ContingencyAnalysisCPP(grid)
        ca.add_contingency_shunt(shunt_id)
        ca.compute(v0, 30, 1e-10)
        self.assertFalse(ca.converged_mask()[0])
        self.assertIn("ShuntVoltageControl", ca.get_batch_unsupported_reasons()[0])

        # ... and solved after the batch by the fallback, as a one-off solve with the loop does
        ca = ContingencyAnalysisCPP(grid)
        ca.batch_unsupported_fallback = True
        ca.add_contingency_shunt(shunt_id)
        ca.compute(v0, 30, 1e-10)
        self.assertTrue(ca.converged_mask()[0])
        self.assertIn("ShuntVoltageControl", ca.get_batch_unsupported_reasons()[0])
        ref = copy.deepcopy(grid)
        ref.deactivate_shunt(shunt_id)
        V_ref = ref.ac_pf(v0.copy(), 30, 1e-10)
        self.assertGreater(V_ref.shape[0], 0)
        np.testing.assert_allclose(ca.get_voltages()[0], V_ref, atol=1e-8, rtol=0.)


class TestShuntReconnected(unittest.TestCase):
    """The outer loops reserve a shunt's control only when the Newton rebuilds, and a shunt
    going on or off only changes the diagonal of Ybus, which asks for no rebuild. A regulating
    shunt that was off when the controls were reserved came back with no control: the solve
    said STABLE, the sections did not move and the bus missed its target.

    Every test compares a grid that went through the sequence with a fresh one solved once in
    the same final state: the history must not show in the result."""

    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is not installed")

    def _grid(self, bus="B9", target_kv=13.2, sections=(3, 1), off=(), moved=()):
        grid = _grid(_net(bus, target_kv, sections))
        grid.change_algorithm("NROuter_KLU")
        grid.clear_outer_loops()
        grid.add_outer_loop(ShuntVoltageControl())
        for shunt_id in self._ids(grid, off):
            grid.deactivate_shunt(shunt_id)
        self._move(grid, moved)
        return grid

    def _move(self, grid, names):
        """Move shunts to a bus next to the one they started on, energized by a line."""
        for shunt_id in self._ids(grid, names):
            bus = grid.get_shunts()[shunt_id].regulated_bus  # they all start on the bus they regulate
            grid.change_bus_shunt(shunt_id, next(l.bus2_id for l in grid.get_lines() if l.bus1_id == bus))

    @staticmethod
    def _ids(grid, names):
        all_names = [s.name for s in grid.get_shunts()]
        return [all_names.index(name) for name in names]

    def _solve(self, grid, v0=None):
        v0 = np.full(grid.total_bus(), 1.0 + 0j) if v0 is None else v0
        V = grid.ac_pf(v0.copy(), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        self.assertEqual(grid.get_algo().get_outer_loop_stats().status.name, "STABLE")
        return V

    @staticmethod
    def _analyzes(grid):
        return grid.get_algo().get_linear_solver_stats().nb_analyze

    @staticmethod
    def _counts(grid):
        res = {s.name: s.res_section_count for s in grid.get_shunts()}
        return [res[name] for name in SHUNTS]

    def _assert_same_as_fresh(self, grid, V, off=(), **kwargs):
        fresh = self._grid(off=off, **kwargs)
        V_fresh = self._solve(fresh)
        self.assertEqual(self._counts(grid), self._counts(fresh))
        np.testing.assert_allclose(V, V_fresh, atol=1e-8, rtol=0.)

    def test_all_shunts_off_at_the_first_solve(self):
        # no group at all was reserved: the bus had no control to come back to
        grid = self._grid(off=SHUNTS)
        self._solve(grid)
        for shunt_id in self._ids(grid, SHUNTS):
            grid.reactivate_shunt(shunt_id)
        self._assert_same_as_fresh(grid, self._solve(grid))

    def test_one_shunt_off_at_the_first_solve(self):
        # a group was reserved for the bus, without the shunt that comes back
        for off in SHUNTS:
            with self.subTest(off=off):
                grid = self._grid(off=[off])
                self._solve(grid)
                grid.reactivate_shunt(self._ids(grid, [off])[0])
                self._assert_same_as_fresh(grid, self._solve(grid))

    def test_reconnected_from_the_previous_voltages(self):
        # a warm start, as an agent solving step after step does
        grid = self._grid(off=SHUNTS)
        V = self._solve(grid)
        for shunt_id in self._ids(grid, SHUNTS):
            grid.reactivate_shunt(shunt_id)
        self._assert_same_as_fresh(grid, self._solve(grid, V))

    def test_several_cycles(self):
        grid = self._grid(off=SHUNTS)
        ids = self._ids(grid, SHUNTS)
        for _ in range(3):
            self._solve(grid)
            for shunt_id in ids:
                grid.reactivate_shunt(shunt_id)
            self._assert_same_as_fresh(grid, self._solve(grid))
            for shunt_id in ids:
                grid.deactivate_shunt(shunt_id)
            self._assert_same_as_fresh(grid, self._solve(grid), off=SHUNTS)

    def test_one_shunt_tripped_after_a_solve(self):
        # the other direction: the group was reserved with a shunt that is now gone
        for off in SHUNTS:
            with self.subTest(off=off):
                grid = self._grid()
                self._solve(grid)
                grid.deactivate_shunt(self._ids(grid, [off])[0])
                self._assert_same_as_fresh(grid, self._solve(grid), off=[off])

    def test_hidden_by_a_generator_stays_hidden(self):
        # bus 6 is a generator's: reconnected, the shunts still have nothing to do
        kwargs = dict(bus="B6", target_kv=12.8, sections=(0, 2))
        grid = self._grid(off=SHUNTS, **kwargs)
        self._solve(grid)
        for shunt_id in self._ids(grid, SHUNTS):
            grid.reactivate_shunt(shunt_id)
        V = self._solve(grid)
        self.assertEqual(self._counts(grid), [0, 2])
        self._assert_same_as_fresh(grid, V, **kwargs)

    def test_rebuilds_once(self):
        # one rebuild for the reconnection, none for the solves that follow it
        grid = self._grid(off=SHUNTS)
        self._solve(grid)
        before = self._analyzes(grid)
        for shunt_id in self._ids(grid, SHUNTS):
            grid.reactivate_shunt(shunt_id)
        self._solve(grid)
        after = self._analyzes(grid)
        self.assertEqual(after, before + 1)
        self._solve(grid)
        self.assertEqual(self._analyzes(grid), after)

    def test_all_shunts_moved_after_a_solve(self):
        # moving a shunt between two energized buses renumbers no bus and adds no Ybus entry: the
        # controls reserved at the previous bus stayed, and the bus missed its target
        grid = self._grid()
        self._solve(grid)
        self._move(grid, SHUNTS)
        self._assert_same_as_fresh(grid, self._solve(grid), moved=SHUNTS)

    def test_one_shunt_moved_after_a_solve(self):
        # the shunts of the bus become two controllers of the same regulated bus
        for name in SHUNTS:
            with self.subTest(moved=name):
                grid = self._grid()
                self._solve(grid)
                self._move(grid, [name])
                self._assert_same_as_fresh(grid, self._solve(grid), moved=[name])

    def test_moved_from_the_previous_voltages(self):
        grid = self._grid()
        V = self._solve(grid)
        self._move(grid, SHUNTS)
        self._assert_same_as_fresh(grid, self._solve(grid, V), moved=SHUNTS)

    def test_moved_there_and_back(self):
        grid = self._grid()
        V_start = self._solve(grid)
        ids = self._ids(grid, SHUNTS)
        home = [grid.get_shunts()[k].bus_id for k in ids]
        self._move(grid, SHUNTS)
        self._assert_same_as_fresh(grid, self._solve(grid), moved=SHUNTS)
        for shunt_id, bus in zip(ids, home):
            grid.change_bus_shunt(shunt_id, bus)
        V = self._solve(grid)
        self._assert_same_as_fresh(grid, V)
        np.testing.assert_allclose(V, V_start, atol=1e-8, rtol=0.)

    def test_move_rebuilds_once(self):
        grid = self._grid()
        self._solve(grid)
        before = self._analyzes(grid)
        self._move(grid, SHUNTS)
        self._solve(grid)
        after = self._analyzes(grid)
        self.assertEqual(after, before + 1)
        self._solve(grid)
        self.assertEqual(self._analyzes(grid), after)

    def test_move_to_the_same_bus_is_no_rebuild(self):
        grid = self._grid()
        self._solve(grid)
        before = self._analyzes(grid)
        for shunt_id in self._ids(grid, SHUNTS):
            grid.change_bus_shunt(shunt_id, grid.get_shunts()[shunt_id].bus_id)
        self._solve(grid)
        self.assertEqual(self._analyzes(grid), before)

    def test_no_rebuild_when_moving_a_shunt_that_does_not_regulate(self):
        grid = self._grid()
        shunt_id = self._ids(grid, ["SH-A"])[0]
        shunt = grid.get_shunts()[shunt_id]
        bus = shunt.regulated_bus
        grid.set_shunt_section_regulation(shunt_id, False, shunt.target_vm_pu, shunt.target_deadband_pu, bus)
        self._solve(grid)
        before = self._analyzes(grid)
        grid.change_bus_shunt(shunt_id, next(l.bus2_id for l in grid.get_lines() if l.bus1_id == bus))
        self._solve(grid)
        self.assertEqual(self._analyzes(grid), before)

    def test_no_rebuild_for_a_shunt_that_does_not_regulate(self):
        # switching such a shunt changes no control: the Newton keeps its factorization
        grid = self._grid()
        for shunt_id, shunt in zip(self._ids(grid, SHUNTS), [grid.get_shunts()[k] for k in self._ids(grid, SHUNTS)]):
            grid.set_shunt_section_regulation(shunt_id, False, shunt.target_vm_pu, shunt.target_deadband_pu,
                                              shunt.regulated_bus)
        self._solve(grid)
        before = self._analyzes(grid)
        shunt_id = self._ids(grid, ["SH-A"])[0]
        grid.deactivate_shunt(shunt_id)
        self._solve(grid)
        grid.reactivate_shunt(shunt_id)
        self._solve(grid)
        self.assertEqual(self._analyzes(grid), before)


if __name__ == "__main__":
    unittest.main()
