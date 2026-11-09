# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Contingencies of ``ContingencyAnalysisCPP`` tripping what is neither a branch, a shunt, a
load nor a generator: a batch row against the one-off powerflow of the grid with the element
really deactivated."""

import copy
import unittest
import warnings

import numpy as np
import pandapower as pp
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP
from lightsim2grid.algorithm import AlgorithmType

MAX_IT = 30
TOL = 1e-10
LEAF_BUS = 7       # case14: a leaf, behind one transformer
GEN_REMOTE = 3     # the generator on the leaf bus, made a remote regulator where a row must be refused
REG_BUS = 9


def _leaf_branch(grid):
    n_line = len(grid.get_lines())
    branch = [n_line + t.id for t in grid.get_trafos() if LEAF_BUS in (t.bus1_id, t.bus2_id)]
    assert len(branch) == 1
    return branch[0]


class _ContingencyMoreElements(unittest.TestCase):
    """What every kind of element has to do. A subclass says how to build the grid (`_net`
    gives the pandapower network, `_setup_grid` what can only be set on the grid) and how
    to trip an element of it: `ARG` names the list of `add_contingency`, `_deactivate` does it
    on a grid."""
    ARG = None

    @staticmethod
    def _net():
        return pn.case14()

    def _setup_grid(self, grid):
        pass

    def _deactivate(self, grid, el_id):
        raise NotImplementedError

    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(self._net())
        self._setup_grid(self.grid)
        self.n_line = len(self.grid.get_lines())
        self.v0 = np.ones(self.grid.total_bus(), dtype=complex)

    def _reference(self, el_ids=(), branches=(), gens=(), dc=False, grid=None):
        ref = copy.deepcopy(self.grid if grid is None else grid)
        for br in branches:
            if br < self.n_line:
                ref.deactivate_powerline(int(br))
            else:
                ref.deactivate_trafo(int(br) - self.n_line)
        for gen_id in gens:
            ref.deactivate_gen(int(gen_id))
        for el_id in el_ids:
            self._deactivate(ref, int(el_id))
        ref.tell_solver_need_reset()
        V = ref.dc_pf(1.0 * self.v0, MAX_IT, TOL) if dc else ref.ac_pf(1.0 * self.v0, MAX_IT, TOL)
        self.assertGreater(V.shape[0], 0)
        solved = np.zeros(V.shape[0], dtype=bool)
        solved[np.asarray(ref.id_dc_solver_to_me() if dc else ref.id_ac_solver_to_me())] = True
        V[~solved] = 0.
        return V

    def _analysis(self, conts, grid=None, dc=False, **attributes):
        ca = ContingencyAnalysisCPP(self.grid if grid is None else grid)
        if dc:
            ca.change_algorithm(AlgorithmType.DC_SparseLU)
        for name, value in attributes.items():
            setattr(ca, name, value)
        for el_ids, branches in conts:
            ca.add_contingency(branch_ids=list(branches), **{self.ARG: list(el_ids)})
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        return ca

    def _row(self, ca, el_ids, branches=()):
        pos = self._arg_position()
        for row, cont in enumerate(ca.my_contingencies()):
            if cont[pos] == sorted(el_ids) and cont[0] == sorted(branches):
                return row
        raise AssertionError(f"no row for {el_ids}, {branches}")

    def _arg_position(self):
        return ["branch_ids", "shunt_ids", "load_ids", "gen_ids", "branch_ids_side1", "branch_ids_side2",
                "sgen_ids", "storage_ids", "svc_ids", "hvdc_ids"].index(self.ARG)


class TestStaticGenerators(_ContingencyMoreElements):
    ARG = "sgen_ids"
    TRIPS = [([0], []), ([1], []), ([0, 1], []), ([0], [5])]

    @staticmethod
    def _net():
        net = pn.case14()
        for bus, p_mw, q_mvar in ((3, 8., 2.), (LEAF_BUS, 5., 0.)):
            pp.create_sgen(net, bus, p_mw=p_mw, q_mvar=q_mvar)
        return net

    def _deactivate(self, grid, el_id):
        grid.deactivate_sgen(el_id)

    def test_matches_a_one_off_powerflow(self):
        ca = self._analysis(self.TRIPS)
        V = np.asarray(ca.get_voltages())
        for el_ids, branches in self.TRIPS:
            row = self._row(ca, el_ids, branches)
            self.assertTrue(ca.converged_mask()[row])
            np.testing.assert_allclose(V[row], self._reference(el_ids, branches), atol=1e-8, rtol=0.)

    def test_dc(self):
        ca = self._analysis(self.TRIPS, dc=True)
        V = np.asarray(ca.get_voltages())
        for el_ids, branches in self.TRIPS:
            np.testing.assert_allclose(V[self._row(ca, el_ids, branches)],
                                       self._reference(el_ids, branches, dc=True), atol=1e-8, rtol=0.)

    def test_registration(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_contingency(sgen_ids=[1, 0])
        ca.add_contingency(sgen_ids=[0, 1])
        self.assertEqual([c[6] for c in ca.my_contingencies()], [[0, 1]])
        self.assertTrue(ca.remove_contingency(sgen_ids=[0, 1]))
        self.assertEqual(len(ca.my_contingencies()), 0)
        with self.assertRaises(RuntimeError):
            ca.add_contingency(sgen_ids=[len(self.grid.get_static_generators())])

    def test_redistribute_slack_shares_what_it_loses(self):
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 1.)
        ca = self._analysis([([0], [])], redistribute_slack=True)
        ref = copy.deepcopy(self.grid)
        ref.deactivate_sgen(0)
        ref.redistribute_active_power(8.)
        ref.tell_solver_need_reset()
        V_ref = ref.ac_pf(1.0 * self.v0, MAX_IT, TOL)
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], V_ref, atol=1e-6, rtol=0.)

    def test_a_tripped_element_is_not_lost_twice_when_stranded(self):
        # sgen 1 stands on the leaf bus: stranded by the branch, tripped or not, it is lost once
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 1.)
        branch = _leaf_branch(self.grid)
        ca = self._analysis([([1], [branch])], handle_disconnected_grid=True, redistribute_slack=True)
        ca_ref = ContingencyAnalysisCPP(self.grid)
        ca_ref.handle_disconnected_grid = True
        ca_ref.redistribute_slack = True
        ca_ref.add_n1(branch)
        ca_ref.compute(1.0 * self.v0, MAX_IT, TOL)
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], np.asarray(ca_ref.get_voltages())[0],
                                   atol=1e-8, rtol=0.)

    def test_fallback(self):
        # a remote regulator takes the row out of the batch: the fallback trips the rest too
        self.grid.set_gen_regulated_bus(GEN_REMOTE, REG_BUS)
        self.grid.tell_solver_need_reset()
        ca = ContingencyAnalysisCPP(self.grid)
        ca.batch_unsupported_fallback = True
        ca.add_contingency(gen_ids=[GEN_REMOTE], sgen_ids=[0])
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertTrue(ca.converged_mask()[0])
        self.assertIn("remote bus", ca.get_batch_unsupported_reasons()[0])
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], self._reference([0], gens=[GEN_REMOTE]),
                                   atol=1e-8, rtol=0.)


if __name__ == "__main__":
    unittest.main()
