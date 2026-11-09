# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Follow-up of the pandapower converter work of #289: the shunts whose power per step comes
from ``net.shunt_characteristic_table`` (``step_dependency_table``), against ``pp.runpp``.

The table is deliberately not proportional to the step, and has steps changing only one of
p and q: moving a shunt must update what the solver cached for either of them."""

import copy
import unittest
import warnings

import numpy as np
import pandas as pd
import pandapower as pp

from lightsim2grid.network import init_from_pandapower
from global_var_tests import MAX_PP2_DATAREADER, CURRENT_PP_VERSION
from test_pp_converter_trafo_vs_runpp import make_net, pp_voltages, ls_voltages, TOL_V

try:
    from pandapower.control import DiscreteShuntController
    HAS_SHUNT_CONTROLLER = True
except ImportError:
    HAS_SHUNT_CONTROLLER = False

# step: (p_mw, q_mvar). 1 -> 2 moves only p, 2 -> 3 moves only q
ROWS = {0: (0., 0.), 1: (0.02, 2.), 2: (0.05, 2.), 3: (0.05, 5.), 4: (0.1, 9.)}
LOAD_BUS = 2


def set_table(net, rows=None, id_characteristic=0):
    """add a characteristic to the table of the net"""
    rows = ROWS if rows is None else rows
    new = pd.DataFrame([{"id_characteristic": id_characteristic, "step": step, "p_mw": p, "q_mvar": q}
                        for step, (p, q) in rows.items()])
    if "shunt_characteristic_table" in net and net.shunt_characteristic_table.shape[0]:
        new = pd.concat([net.shunt_characteristic_table, new], ignore_index=True)
    net["shunt_characteristic_table"] = new


def add_tabular_shunt(net, step=2, max_step=4, id_characteristic=0, bus=LOAD_BUS, **kwargs):
    """a shunt reading the table; its own p_mw / q_mvar are not what pandapower uses"""
    return pp.create_shunt(net, bus, q_mvar=1., p_mw=0.003, step=step, max_step=max_step,
                           step_dependency_table=True, id_characteristic_table=id_characteristic, **kwargs)


def tabular_net(step=2, max_step=4, rows=None, **trafo_kw):
    net = make_net(**trafo_kw)
    set_table(net, rows)
    add_tabular_shunt(net, step=step, max_step=max_step)
    return net


def convert(net):
    net = copy.deepcopy(net)
    net._options = dict(trafo_model="t")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return init_from_pandapower(net, pp_orig_file="pandapower_v3")


def at_step(net, step, shunt=0):
    net = copy.deepcopy(net)
    net.shunt.loc[net.shunt.index[shunt], "step"] = step
    return net


class _Base(unittest.TestCase):
    def setUp(self):
        if CURRENT_PP_VERSION <= MAX_PP2_DATAREADER:
            self.skipTest("pandapower 2 has no step_dependency_table")

    def assert_same_as_pp(self, net, model=None):
        """same complex voltages as pp.runpp; `model` when it was moved rather than converted"""
        V_pp = pp_voltages(net)
        V_ls = ls_voltages(convert(net) if model is None else model, V_pp.shape[0])
        np.testing.assert_allclose(V_ls, V_pp, rtol=0., atol=TOL_V)


class TestTabularShuntVsRunpp(_Base):
    def test_every_step(self):
        """the power of a tabular shunt is the row of its step, not q_mvar * step"""
        for step in ROWS:
            with self.subTest(step=step):
                self.assert_same_as_pp(tabular_net(step=step))

    def test_rated_voltage_scales_the_table(self):
        """as p_mw and q_mvar, the rows are given at the rated voltage of the shunt"""
        net = make_net()
        set_table(net)
        add_tabular_shunt(net, step=3, vn_kv=21.)
        self.assert_same_as_pp(net)

    def test_mixed_with_linear_shunts(self):
        """each shunt reads its own characteristic, or none, by position in the table"""
        net = make_net()
        set_table(net, id_characteristic=7)
        set_table(net, {1: (0.01, -1.), 2: (0.04, -3.)}, id_characteristic=3)
        pp.create_shunt(net, LOAD_BUS, q_mvar=-0.5, p_mw=0.01, step=2, max_step=4, index=12)
        add_tabular_shunt(net, step=3, id_characteristic=7, index=5)
        add_tabular_shunt(net, step=2, max_step=2, id_characteristic=3, bus=1, index=8)
        self.assert_same_as_pp(net)
        model = convert(net)
        self.assertEqual([sh.max_section_count for sh in model.get_shunts()], [4, 4, 2])

    def test_out_of_service(self):
        net = tabular_net(step=3)
        net.shunt["in_service"] = False
        self.assert_same_as_pp(net)

    def test_missing_power_is_zero(self):
        """pandapower reads a NaN in the table as 0"""
        rows = dict(ROWS)
        rows[2] = (float("nan"), 2.)
        self.assert_same_as_pp(tabular_net(step=2, rows=rows))

    def test_dc(self):
        """the active power of the row also enters the DC powerflow"""
        for step in (1, 4):
            with self.subTest(step=step):
                net = tabular_net(step=step, tap_pos=0)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    pp.rundcpp(net, calculate_voltage_angles=True)
                model = convert(net)
                V = model.dc_pf(np.full(model.total_bus(), 1. + 0j), 1, 1e-8)
                self.assertGreater(V.shape[0], 0)
                np.testing.assert_allclose(np.angle(V[:3], deg=True), net.res_bus["va_degree"].to_numpy(),
                                           rtol=0., atol=1e-8)


class TestMovedStepVsRunpp(_Base):
    """converted at one step, moved to another with `change_shunt_section_count`: the grid is
    what pp.runpp computes at the new step"""

    def test_sections_are_the_table(self):
        shunt = convert(tabular_net(step=3)).get_shunts()[0]
        self.assertEqual(shunt.section_count, 3)
        self.assertEqual(shunt.max_section_count, 4)
        self.assertAlmostEqual(shunt.target_p_mw, 0.05)
        self.assertAlmostEqual(shunt.target_q_mvar, 5.)

    def test_moved_step(self):
        """one model moved through every step, in an order going through p-only and q-only
        changes and through step 0, so that a stale Ybus would show"""
        net = tabular_net(step=2)
        model = convert(net)
        for step in (1, 2, 3, 0, 4, 3, 2, 1, 0):
            with self.subTest(step=step):
                model.change_shunt_section_count(0, step)
                self.assert_same_as_pp(at_step(net, step), model)

    def test_moved_step_scaled_by_rated_voltage(self):
        net = make_net()
        set_table(net)
        add_tabular_shunt(net, step=1, vn_kv=19.)
        model = convert(net)
        for step in (4, 0, 2):
            with self.subTest(step=step):
                model.change_shunt_section_count(0, step)
                self.assert_same_as_pp(at_step(net, step), model)

    def test_table_starting_at_step_one(self):
        """no row for step 0 (CIM non linear compensators): steps 1 .. n are still the sections"""
        rows = {k: v for k, v in ROWS.items() if k != 0}
        net = tabular_net(step=1, rows=rows)
        model = convert(net)
        for step in (4, 2):
            with self.subTest(step=step):
                model.change_shunt_section_count(0, step)
                self.assert_same_as_pp(at_step(net, step), model)

    def test_moved_with_the_transformer_tap(self):
        """the tap of the transformer and the step of the shunt, moved together on one model"""
        net = tabular_net(step=1)
        model = convert(net)
        for tap, step in ((-4, 3), (6, 0), (-9, 4), (3, 2)):
            with self.subTest(tap=tap, step=step):
                model.change_trafo_ratio_tap(0, tap)
                model.change_shunt_section_count(0, step)
                moved = at_step(net, step)
                moved.trafo["tap_pos"] = tap
                self.assert_same_as_pp(moved, model)

    def test_max_step_bounds_the_sections(self):
        """the sections are what max_step allows, even when the table goes further"""
        net = tabular_net(step=1, max_step=2)
        model = convert(net)
        self.assertEqual(model.get_shunts()[0].max_section_count, 2)
        model.change_shunt_section_count(0, 2)
        self.assert_same_as_pp(at_step(net, 2), model)
        with self.assertRaises(RuntimeError):
            model.change_shunt_section_count(0, 3)

    def test_table_shorter_than_max_step(self):
        """a step without a row does not exist: the sections stop with the table"""
        rows = {k: v for k, v in ROWS.items() if k <= 2}
        shunt = convert(tabular_net(step=1, max_step=4, rows=rows)).get_shunts()[0]
        self.assertEqual(shunt.max_section_count, 2)


class TestNoSections(_Base):
    """what a section table cannot say: the shunt is right where it is, cannot be moved, and
    the user is told"""

    def _no_sections(self, net):
        net = copy.deepcopy(net)
        net._options = dict(trafo_model="t")
        with self.assertWarnsRegex(UserWarning, "no sections"):
            model = init_from_pandapower(net, pp_orig_file="pandapower_v3")
        self.assertEqual(model.get_shunts()[0].max_section_count, 0)
        self.assert_same_as_pp(net)

    def test_power_at_step_zero(self):
        """no section on means no power, in a section table"""
        rows = dict(ROWS)
        rows[0] = (0.01, 1.)
        self._no_sections(tabular_net(step=2, rows=rows))
        self._no_sections(tabular_net(step=0, rows=rows))

    def test_hole_in_the_table(self):
        rows = {k: v for k, v in ROWS.items() if k != 3}
        self._no_sections(tabular_net(step=2, rows=rows))

    def test_step_above_max_step(self):
        """pp.runpp only reads the table, max_step does not limit it"""
        self._no_sections(tabular_net(step=4, max_step=2))

    def test_no_max_step(self):
        net = tabular_net(step=2)
        net.shunt["max_step"] = np.nan
        self._no_sections(net)


class TestRefused(_Base):
    """the tables pp.runpp cannot read either: no silent guess at the power of the shunt"""

    def _refused(self, net, msg):
        with self.assertRaisesRegex(RuntimeError, msg):
            convert(net)

    def test_step_not_in_the_table(self):
        self._refused(tabular_net(step=7), "no row")

    def test_no_table(self):
        net = tabular_net()
        del net["shunt_characteristic_table"]
        self._refused(net, "shunt_characteristic_table")

    def test_no_characteristic(self):
        net = tabular_net()
        net.shunt["id_characteristic_table"] = pd.NA
        self._refused(net, "id_characteristic_table")

    def test_duplicated_step(self):
        net = tabular_net(step=2)
        set_table(net, {2: (0.5, 0.5)})
        self._refused(net, "several rows")


@unittest.skipUnless(HAS_SHUNT_CONTROLLER, "this pandapower has no shunt controller")
class TestController(_Base):
    def test_controller_regulates_through_the_table(self):
        """a DiscreteShuntController on a tabular shunt: the outer loop walks the table, and
        lands on the step pandapower's controller finds"""
        from lightsim2grid.algorithm import ShuntVoltageControl
        rows = {0: (0., 0.), 1: (0., -1.), 2: (0., -2.5), 3: (0., -6.), 4: (0., -8.)}
        net = tabular_net(step=0, rows=rows)
        net.shunt["q_mvar"] = -1.  # the controller takes its direction from the sign of this column
        ref = copy.deepcopy(net)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pp.runpp(ref, tolerance_mva=1e-10, lightsim2grid=False, numba=False)
        target = float(ref.res_bus.at[LOAD_BUS, "vm_pu"]) + 0.01
        DiscreteShuntController(net, 0, vm_set_pu=target, tol=0.005)
        model = convert(net)
        shunt = model.get_shunts()[0]
        self.assertTrue(shunt.regulating)
        self.assertEqual(shunt.max_section_count, 4)
        model.change_algorithm("NROuter_KLU" if "NROuter_KLU" in model.available_algorithm_names()
                               else "NROuter_SparseLU")
        model.add_outer_loop(ShuntVoltageControl())
        V = model.ac_pf(np.full(model.total_bus(), 1. + 0j), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pp.runpp(net, run_control=True, tolerance_mva=1e-10, lightsim2grid=False, numba=False)
        self.assertEqual(model.get_shunts()[0].res_section_count, int(net.shunt["step"].iloc[0]))
        self.assertGreater(int(net.shunt["step"].iloc[0]), 0)
        self.assertAlmostEqual(abs(V[LOAD_BUS]), net.res_bus.at[LOAD_BUS, "vm_pu"], delta=1e-8)


if __name__ == "__main__":
    unittest.main()
