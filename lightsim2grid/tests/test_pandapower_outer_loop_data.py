# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The data the NROuter_* outer loops read, as the pandapower converter fills it (issue 273):
the ratio tap changers and what they regulate, the shunt sections and what they regulate,
and the units that can participate in the slack."""

import unittest
import warnings

import numpy as np
import pandapower as pp
from pandapower.control import DiscreteTapControl, ContinuousTapControl

from lightsim2grid.network import init_from_pandapower
from lightsim2grid.algorithm import DistributedSlack, TransformerVoltageControl, ShuntVoltageControl
from lightsim2grid.algorithm._outer_loop_checks import distributed_slack_is_idle, warn_if_distributed_slack_idle

try:
    from pandapower.control import DiscreteShuntController
    HAS_SHUNT_CONTROLLER = True
except ImportError:
    HAS_SHUNT_CONTROLLER = False

PP_ORIG_FILES = ("pandapower_v2", "pandapower_v3")


def _trafo_net(tap_pos=0, tap_neutral=0, tap_side="hv", shunt_step=1, shunt_vn_kv=20.):
    """ext_grid -- 110/20 kV transformer -- line -- load + shunt"""
    net = pp.create_empty_network()
    b0 = pp.create_bus(net, 110.)
    b1 = pp.create_bus(net, 20.)
    b2 = pp.create_bus(net, 20.)
    pp.create_ext_grid(net, b0, vm_pu=1.02)
    t = pp.create_transformer(net, b0, b1, std_type="25 MVA 110/20 kV")
    net.trafo.loc[t, ["tap_pos", "tap_neutral", "tap_side", "tap_min", "tap_max"]] = \
        [tap_pos, tap_neutral, tap_side, tap_neutral - 9, tap_neutral + 9]
    net.trafo["shift_degree"] = 0.  # keeps a flat start
    pp.create_line(net, b1, b2, 3., "NA2XS2Y 1x95 RM/25 12/20 kV")
    pp.create_load(net, b2, 15., 5.)
    pp.create_shunt(net, b2, q_mvar=-1.5, p_mw=0.01, step=shunt_step, max_step=4, vn_kv=shunt_vn_kv)
    return net


def _convert(net, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return init_from_pandapower(net, **kwargs)


def _vm(grid, max_it=30, tol=1e-10):
    V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), max_it, tol)
    assert V.shape[0] > 0, "the powerflow diverged"
    return np.abs(V)


def _runpp(net, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pp.runpp(net, lightsim2grid=False, numba=False, tolerance_mva=1e-10, **kwargs)


def _with_loops(grid, *loops):
    grid.change_algorithm("NROuter_KLU" if "NROuter_KLU" in grid.available_algorithm_names() else "NROuter_SparseLU")
    for loop in loops:
        grid.add_outer_loop(loop)
    return grid


class TestTapChangers(unittest.TestCase):
    def test_table_spans_tap_min_tap_max(self):
        grid = _convert(_trafo_net(tap_pos=2, tap_neutral=3))
        taps = grid.get_trafos()[0]
        self.assertEqual(taps.ratio_low_tap, -6)
        self.assertEqual(taps.ratio_high_tap, 12)
        self.assertEqual(taps.ratio_tap_position, 2)

    def test_moved_tap_is_a_fresh_conversion(self):
        """moving the tap gives the transformer the converter builds at that tap"""
        for orig in PP_ORIG_FILES:
            for side in ("hv", "lv"):
                for neutral in (0, 3):
                    grid = _convert(_trafo_net(neutral + 2, neutral, side), pp_orig_file=orig)
                    for pos in (-9, -1, 0, 4, 9):
                        pos += neutral
                        with self.subTest(orig=orig, side=side, neutral=neutral, pos=pos):
                            grid.change_trafo_ratio_tap(0, pos)
                            fresh = _convert(_trafo_net(pos, neutral, side), pp_orig_file=orig)
                            np.testing.assert_allclose(_vm(grid), _vm(fresh), atol=1e-12, rtol=0.)

    def test_tap_neutral_matches_pandapower(self):
        """tap_neutral != 0 is no longer refused, and is counted from as pandapower counts"""
        for side in ("hv", "lv"):
            for pos in (-6, 3, 10):
                with self.subTest(side=side, pos=pos):
                    net = _trafo_net(pos, 3, side)
                    vm = _vm(_convert(net, pp_orig_file="pandapower_v2"))
                    _runpp(net)
                    np.testing.assert_allclose(vm, net.res_bus["vm_pu"].to_numpy(), atol=1e-8, rtol=0.)

    def test_no_table_without_a_range(self):
        net = _trafo_net()
        net.trafo["tap_min"] = np.nan
        grid = _convert(net)
        self.assertFalse(grid.get_trafos()[0].has_ratio_tap_changer)

    def test_no_table_with_a_tap_step_degree(self):
        net = _trafo_net()
        net.trafo["tap_step_degree"] = 2.
        grid = _convert(net)
        self.assertFalse(grid.get_trafos()[0].has_ratio_tap_changer)

    def test_no_table_with_a_non_integer_tap(self):
        net = _trafo_net(tap_pos=1.5)
        grid = _convert(net)
        self.assertFalse(grid.get_trafos()[0].has_ratio_tap_changer)

    def test_no_controller_no_regulation(self):
        """the tap only moves with a controller"""
        grid = _with_loops(_convert(_trafo_net()), TransformerVoltageControl())
        _vm(grid)
        self.assertEqual(grid.get_trafos()[0].res_ratio_tap_position, 0)

    def test_continuous_tap_control(self):
        """the loop solves the ratio continuously then rounds it: the closest tap to pandapower's continuous one"""
        for side in ("hv", "lv"):
            for target in (1.04, 0.97):
                with self.subTest(side=side, target=target):
                    net = _trafo_net(0, 0, side)
                    ContinuousTapControl(net, 0, vm_set_pu=target, side="lv", tol=1e-6)
                    grid = _with_loops(_convert(net), TransformerVoltageControl())
                    _vm(grid, tol=1e-8)
                    _runpp(net, run_control=True)
                    expected = int(np.round(net.trafo["tap_pos"].iloc[0]))
                    self.assertNotEqual(expected, 0)
                    self.assertEqual(grid.get_trafos()[0].res_ratio_tap_position, expected)

    def test_discrete_tap_control(self):
        """the voltage ends within the band of the controller (the tap may differ from
        pandapower's, which stops at the first tap within the band)"""
        for side in ("hv", "lv"):
            for low, high in ((1.02, 1.04), (0.96, 0.98), (1.04, 1.07)):
                with self.subTest(side=side, band=(low, high)):
                    net = _trafo_net(0, 0, side)
                    DiscreteTapControl(net, 0, vm_lower_pu=low, vm_upper_pu=high, side="lv")
                    grid = _with_loops(_convert(net), TransformerVoltageControl())
                    vm = _vm(grid, tol=1e-8)
                    self.assertNotEqual(grid.get_trafos()[0].res_ratio_tap_position, 0)
                    self.assertGreaterEqual(vm[1], low)
                    self.assertLessEqual(vm[1], high)

    def test_controller_out_of_service(self):
        net = _trafo_net()
        ContinuousTapControl(net, 0, vm_set_pu=1.04, side="lv", in_service=False)
        grid = _with_loops(_convert(net), TransformerVoltageControl())
        _vm(grid)
        self.assertEqual(grid.get_trafos()[0].res_ratio_tap_position, 0)


class TestShunts(unittest.TestCase):
    def test_step_and_vn_kv_match_pandapower(self):
        """the power of a shunt is q_mvar * step * (vn_bus / vn_kv) ** 2, as pandapower's"""
        for step, vn_kv in ((1, 20.), (2, 21.), (0, 20.), (3, 19.)):
            with self.subTest(step=step, vn_kv=vn_kv):
                net = _trafo_net(shunt_step=step, shunt_vn_kv=vn_kv)
                vm = _vm(_convert(net, pp_orig_file="pandapower_v2"))
                _runpp(net)
                np.testing.assert_allclose(vm, net.res_bus["vm_pu"].to_numpy(), atol=1e-8, rtol=0.)

    def test_sections(self):
        net = _trafo_net(shunt_step=2, shunt_vn_kv=21.)
        grid = _convert(net)
        shunt = grid.get_shunts()[0]
        self.assertEqual(shunt.section_count, 2)
        self.assertEqual(shunt.max_section_count, 4)
        ratio = (20. / 21.) ** 2
        self.assertAlmostEqual(shunt.target_q_mvar, -3. * ratio)
        grid.change_shunt_section_count(0, 4)
        self.assertAlmostEqual(grid.get_shunts()[0].target_q_mvar, -6. * ratio)

    def test_no_sections_without_max_step(self):
        net = _trafo_net()
        net.shunt["max_step"] = np.nan
        grid = _convert(net)
        self.assertEqual(grid.get_shunts()[0].max_section_count, 0)

    @unittest.skipUnless(HAS_SHUNT_CONTROLLER, "this pandapower has no shunt controller")
    def test_shunt_controller_matches_pandapower(self):
        net = _trafo_net(shunt_step=0)
        DiscreteShuntController(net, 0, vm_set_pu=0.96, tol=0.005)
        grid = _with_loops(_convert(net), ShuntVoltageControl())
        shunt = grid.get_shunts()[0]
        self.assertTrue(shunt.regulating)
        self.assertAlmostEqual(shunt.target_vm_pu, 0.96)
        self.assertAlmostEqual(shunt.target_deadband_pu, 0.01)
        self.assertEqual(shunt.regulated_bus, 2)
        _vm(grid, tol=1e-8)
        _runpp(net, run_control=True)
        self.assertEqual(grid.get_shunts()[0].res_section_count, int(net.shunt["step"].iloc[0]))


def _slack_net():
    """a slack=True generator and a second one with a slack weight, 1:3"""
    net = pp.create_empty_network()
    b = [pp.create_bus(net, 110.) for _ in range(3)]
    pp.create_gen(net, b[0], p_mw=5., vm_pu=1.02, slack=True, slack_weight=1.)
    pp.create_gen(net, b[2], p_mw=10., vm_pu=1.01, slack_weight=3.)
    pp.create_line(net, b[0], b[1], 30., "149-AL1/24-ST1A 110.0")
    pp.create_line(net, b[1], b[2], 30., "149-AL1/24-ST1A 110.0")
    pp.create_load(net, b[1], 40., 10.)
    return net


class TestSlackParticipation(unittest.TestCase):
    def _solve(self, grid):
        grid.change_algorithm("NROuter_SparseLU")
        grid.clear_outer_loops()
        grid.add_outer_loop(DistributedSlack(slack_bus_p_max_mismatch_mw=1e-4))
        _vm(grid)
        return np.array([el.res_p_mw for el in grid.get_generators()])

    def test_default_single_slack_is_idle(self):
        """pandapower 2: only the slack=True generator participates, the loop has nothing to do"""
        grid = _convert(_slack_net(), pp_orig_file="pandapower_v2")
        self.assertFalse(grid.get_generators()[1].can_participate_slack)
        grid.change_algorithm("NROuter_SparseLU")
        self.assertTrue(distributed_slack_is_idle(grid))
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            self.assertTrue(warn_if_distributed_slack_idle(grid))
        self.assertTrue(any("DistributedSlack" in str(el.message) for el in w))
        p = self._solve(grid)
        self.assertAlmostEqual(p[1], 10.)

    def test_no_warning_with_another_algorithm(self):
        grid = _convert(_slack_net(), pp_orig_file="pandapower_v2")
        grid.change_algorithm("NR_SparseLU")
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.assertFalse(warn_if_distributed_slack_idle(grid))

    def test_slack_weight_shares_as_pandapower(self):
        net = _slack_net()
        _runpp(net, distributed_slack=True)
        ref = net.res_gen["p_mw"].to_numpy()
        for orig in PP_ORIG_FILES:
            with self.subTest(orig=orig):
                grid = _convert(_slack_net(), pp_orig_file=orig, slack_weight_can_participate=True)
                gen = grid.get_generators()[1]
                self.assertTrue(gen.can_participate_slack)
                self.assertAlmostEqual(gen.can_participate_slack_weight, 3.)
                grid.change_algorithm("NROuter_SparseLU")
                self.assertFalse(distributed_slack_is_idle(grid))
                np.testing.assert_allclose(self._solve(grid), ref, atol=1e-3, rtol=0.)

    def test_option_off_keeps_newton_slack(self):
        """the option does not change the slack the Newton itself distributes"""
        for orig in PP_ORIG_FILES:
            with self.subTest(orig=orig):
                ref = _convert(_slack_net(), pp_orig_file=orig)
                grid = _convert(_slack_net(), pp_orig_file=orig, slack_weight_can_participate=True)
                np.testing.assert_array_equal([el.is_slack for el in grid.get_generators()],
                                              [el.is_slack for el in ref.get_generators()])
                np.testing.assert_allclose(_vm(grid), _vm(ref), atol=1e-12, rtol=0.)


class TestBackend(unittest.TestCase):
    def setUp(self):
        try:
            import grid2op  # noqa: F401
        except ImportError:
            self.skipTest("grid2op is not installed")

    def _make(self, **backend_kwargs):
        import grid2op
        from lightsim2grid import LightSimBackend
        return grid2op.make("l2rpn_case14_sandbox", test=True,
                            backend=LightSimBackend(**backend_kwargs),
                            _add_to_name=f"_{type(self).__name__}")

    def test_warns_with_nrouter(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            env = self._make(algo_type="NROuter_SparseLU")
        try:
            self.assertTrue(any("DistributedSlack" in str(el.message) for el in w))
            env.step(env.action_space())
        finally:
            env.close()

    def test_no_warning_with_default_algorithm(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            env = self._make()
        env.close()
        self.assertFalse(any("DistributedSlack" in str(el.message) for el in w))

    def test_no_warning_without_distributed_slack(self):
        """a list of loops without DistributedSlack has nothing to warn about"""
        from lightsim2grid.algorithm import ReactiveLimits
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            env = self._make(algo_type="NROuter_SparseLU", outer_loops=[ReactiveLimits()])
        env.close()
        self.assertFalse(any("DistributedSlack" in str(el.message) for el in w))

    def test_set_outer_loops_warns(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            env = self._make(algo_type="NROuter_SparseLU", outer_loops=[])
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                env.backend.set_outer_loops([DistributedSlack()])
            self.assertTrue(any("DistributedSlack" in str(el.message) for el in w))
        finally:
            env.close()

    def test_loader_kwargs(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            env = self._make(loader_kwargs={"slack_weight_can_participate": True})
        env.close()


if __name__ == "__main__":
    unittest.main()
