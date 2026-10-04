# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""``set_keep_vinit_at_group_controlled_buses``: the buses a voltage-control group
regulates (a remote-regulating generator, a voltage-mode SVC) keep the starting
magnitude instead of being set to their set-point before the solve.

Their magnitude is a Newton-Raphson unknown held by the group's voltage row, so the
option may change the path to the solution, never the solution itself. What it
changes is observable directly on a DC powerflow, which echoes the seeded magnitude
back as its result.
"""

import unittest
import numpy as np

from lightsim2grid.lightsim2grid_cpp import TimeSeriesCPP, ContingencyAnalysisCPP
from lightsim2grid.algorithm import AlgorithmType
from test_voltage_control_batch import (_make, _base_injections, KINDS, GEN_ID, REG_BUS,
                                        SVC_V_SET, MAX_IT, TOL)

# a starting magnitude distinct from every set-point of case14 and of the fixtures
VINIT_VM = 0.97
CONT_LINES = [0, 5, 10]


def _vset(model, kind):
    return model.get_generators()[GEN_ID].target_vm_pu if kind == "remote_gen" else SVC_V_SET


class TestKeepVinitLSGrid(unittest.TestCase):
    def test_default_and_accessors(self):
        _, model = _make("remote_gen")
        self.assertFalse(model.get_keep_vinit_at_group_controlled_buses())
        model.set_keep_vinit_at_group_controlled_buses(True)
        self.assertTrue(model.get_keep_vinit_at_group_controlled_buses())
        self.assertTrue(model.copy().get_keep_vinit_at_group_controlled_buses())
        model.set_keep_vinit_at_group_controlled_buses(False)
        self.assertFalse(model.get_keep_vinit_at_group_controlled_buses())

    def test_seed_at_the_regulated_bus(self):
        # the sloped SVC is left out: a sloped regulator never sets its bus's magnitude
        for kind in ("remote_gen", "svc_flat"):
            with self.subTest(kind=kind):
                net, model = _make(kind)
                vinit = np.full(net.bus.shape[0], VINIT_VM, dtype=np.complex128)
                Vdc = model.dc_pf(vinit, MAX_IT, TOL)
                self.assertAlmostEqual(abs(Vdc[REG_BUS]), _vset(model, kind), places=12)

                model.set_keep_vinit_at_group_controlled_buses(True)
                Vdc = model.dc_pf(vinit, MAX_IT, TOL)
                self.assertAlmostEqual(abs(Vdc[REG_BUS]), VINIT_VM, places=12)

    def test_pv_bus_is_still_set_to_its_target(self):
        # an ordinary PV bus has a fixed magnitude: the option must not touch it
        net, model = _make("remote_gen")
        model.set_keep_vinit_at_group_controlled_buses(True)
        vinit = np.full(net.bus.shape[0], VINIT_VM, dtype=np.complex128)
        Vdc = model.dc_pf(vinit, MAX_IT, TOL)
        gen = model.get_generators()[0]
        self.assertNotEqual(gen.bus_id, REG_BUS)
        self.assertAlmostEqual(abs(Vdc[gen.bus_id]), gen.target_vm_pu, places=12)

    def test_same_ac_solution(self):
        for kind in KINDS:
            with self.subTest(kind=kind):
                net, model = _make(kind)
                vinit = np.full(net.bus.shape[0], VINIT_VM, dtype=np.complex128)
                V_ref = model.ac_pf(vinit, MAX_IT, TOL)
                self.assertGreater(V_ref.shape[0], 0)

                _, model_keep = _make(kind)
                model_keep.set_keep_vinit_at_group_controlled_buses(True)
                V_keep = model_keep.ac_pf(vinit, MAX_IT, TOL)
                self.assertGreater(V_keep.shape[0], 0)
                self.assertLess(np.abs(V_keep - V_ref).max(), 1e-8)


class TestKeepVinitBatch(unittest.TestCase):
    def test_inherited_from_the_grid(self):
        _, model = _make("remote_gen")
        self.assertFalse(ContingencyAnalysisCPP(model).keep_vinit_at_group_controlled_buses)
        self.assertFalse(TimeSeriesCPP(model).keep_vinit_at_group_controlled_buses)
        model.set_keep_vinit_at_group_controlled_buses(True)
        self.assertTrue(ContingencyAnalysisCPP(model).keep_vinit_at_group_controlled_buses)
        self.assertTrue(TimeSeriesCPP(model).keep_vinit_at_group_controlled_buses)

    def test_setter_leaves_the_source_grid_alone(self):
        _, model = _make("remote_gen")
        ca = ContingencyAnalysisCPP(model)
        ca.keep_vinit_at_group_controlled_buses = True
        self.assertTrue(ca.keep_vinit_at_group_controlled_buses)
        self.assertFalse(model.get_keep_vinit_at_group_controlled_buses())

    def test_contingency_analysis_same_results(self):
        for kind in KINDS:
            with self.subTest(kind=kind):
                net, model = _make(kind)
                vinit = np.full(net.bus.shape[0], VINIT_VM, dtype=np.complex128)
                res = []
                for keep in (False, True):
                    ca = ContingencyAnalysisCPP(model)
                    ca.keep_vinit_at_group_controlled_buses = keep
                    for line_id in CONT_LINES:
                        ca.add_n1(line_id)
                    ca.compute(vinit, MAX_IT, TOL)
                    res.append(np.array(ca.get_voltages()))
                self.assertGreater(np.abs(res[0]).max(), 0.)
                self.assertLess(np.abs(res[1] - res[0]).max(), 1e-8)

    def test_time_series_same_results(self):
        for kind in KINDS:
            with self.subTest(kind=kind):
                net, model = _make(kind)
                gen_p, load_p, load_q = _base_injections(model)
                scales = np.array([1.0, 1.08])
                gen_mat = np.outer(scales, gen_p)
                lp_mat = np.outer(scales, load_p)
                lq_mat = np.outer(scales, load_q)
                sgen_mat = np.zeros((len(scales), 0))
                vinit = np.full(net.bus.shape[0], VINIT_VM, dtype=np.complex128)
                res = []
                for keep in (False, True):
                    ts = TimeSeriesCPP(model)
                    ts.keep_vinit_at_group_controlled_buses = keep
                    status = ts.compute_Vs(gen_mat, sgen_mat, lp_mat, lq_mat, vinit, MAX_IT, TOL)
                    self.assertEqual(status, 1)
                    res.append(np.array(ts.get_voltages()))
                self.assertLess(np.abs(res[1] - res[0]).max(), 1e-8)

    def test_seed_at_the_regulated_bus_dc(self):
        # a DC batch echoes the seeded magnitude back, as dc_pf does
        for keep in (False, True):
            with self.subTest(keep=keep):
                net, model = _make("remote_gen")
                vinit = np.full(net.bus.shape[0], VINIT_VM, dtype=np.complex128)
                ca = ContingencyAnalysisCPP(model)
                ca.change_algorithm(AlgorithmType.DC_SparseLU)
                ca.keep_vinit_at_group_controlled_buses = keep
                for line_id in CONT_LINES:
                    ca.add_n1(line_id)
                ca.compute(vinit, MAX_IT, TOL)
                Vs = np.array(ca.get_voltages())
                expected = VINIT_VM if keep else _vset(model, "remote_gen")
                solved = [k for k in range(Vs.shape[0]) if np.abs(Vs[k]).max() > 0.]
                self.assertGreater(len(solved), 0)
                for k in solved:
                    self.assertAlmostEqual(abs(Vs[k][REG_BUS]), expected, places=12)

    def test_changing_it_drops_the_base_case(self):
        net, model = _make("remote_gen")
        vinit = np.full(net.bus.shape[0], VINIT_VM, dtype=np.complex128)
        ca = ContingencyAnalysisCPP(model)
        for line_id in CONT_LINES:
            ca.add_n1(line_id)
        ca.compute(vinit, MAX_IT, TOL)
        ca.compute(vinit, MAX_IT, TOL)
        self.assertTrue(ca.base_case_was_reused())
        ca.keep_vinit_at_group_controlled_buses = True
        ca.compute(vinit, MAX_IT, TOL)
        self.assertFalse(ca.base_case_was_reused())
        # setting the value it already has is not a change
        ca.keep_vinit_at_group_controlled_buses = True
        ca.compute(vinit, MAX_IT, TOL)
        self.assertTrue(ca.base_case_was_reused())


if __name__ == "__main__":
    unittest.main()
