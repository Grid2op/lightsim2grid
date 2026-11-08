# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A contingency the batch does not support is left UNSOLVED, the others still run (#280).

Tripping a generator that regulates a remote bus is not supported by the batch yet (its
voltage control lives in Jacobian rows and columns the batch does not reserve per row). It
used to make ``compute()`` raise for the whole batch; that row alone is now reported
``LimitViolationType.UNSOLVED``, with its reason in ``get_unsolved_reasons()``.

Fixture: pandapower case14, generator 3 (on bus 7) regulating bus 9 remotely, as in
``test_voltage_control_batch.py``.
"""

import copy
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.network import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import (ContingencyAnalysisCPP, LimitViolationType,  # noqa: E402
                                             ViolationCategory, ViolationElementType)

MAX_IT = 30
TOL = 1e-11
GEN_REMOTE = 3    # on bus 7 ...
REG_BUS = 9       # ... regulating bus 9
TRAFO_BEHIND = 3  # buses 6-7: bus 7 is behind it, alone
LOCAL_GEN = 1     # a generator regulating its own bus


def _grid():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(pn.case14())
    grid.set_gen_regulated_bus(GEN_REMOTE, REG_BUS)
    grid.tell_solver_need_reset()
    return grid


class TestUnsolvedContingencies(unittest.TestCase):
    def setUp(self):
        self.grid = _grid()
        self.n_line = len(self.grid.get_lines())
        self.v0 = np.ones(self.grid.total_bus(), dtype=complex)
        self.buses = np.asarray(self.grid.id_ac_solver_to_me(), dtype=int)

    def _reference(self, branches=(), gens=()):
        ref = copy.deepcopy(self.grid)
        for br in branches:
            if br < self.n_line:
                ref.deactivate_powerline(int(br))
            else:
                ref.deactivate_trafo(int(br - self.n_line))
        for g in gens:
            ref.deactivate_gen(int(g))
        return ref.ac_pf(1.0 * self.v0, MAX_IT, TOL)

    def _analysis(self, nb_thread=1, handle_disconnected_grid=False):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.compute_limit_violations = True
        ca.handle_disconnected_grid = handle_disconnected_grid
        ca.nb_thread = nb_thread
        for line_id in (0, 3, 7):
            ca.add_n1(line_id)
        ca.add_contingency_gen(GEN_REMOTE)
        ca.add_contingency_gen(LOCAL_GEN)
        ca.add_contingency(branch_ids=[5], gen_ids=[GEN_REMOTE])
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        return ca

    @staticmethod
    def _row(ca, branches=(), gens=()):
        for row, cont in enumerate(ca.my_contingencies()):
            if sorted(cont[0]) == sorted(branches) and sorted(cont[3]) == sorted(gens):
                return row
        raise AssertionError(f"no row for branches {branches}, generators {gens}")

    def _check(self, ca):
        reasons = ca.get_unsolved_reasons()
        conts = ca.my_contingencies()
        self.assertEqual(len(reasons), len(conts))
        unsolved = {self._row(ca, gens=[GEN_REMOTE]), self._row(ca, branches=[5], gens=[GEN_REMOTE])}
        for row, cont in enumerate(conts):
            with self.subTest(row=row, contingency=cont):
                if row in unsolved:
                    self.assertIn("remote bus", reasons[row])
                    self.assertIn(f"generator {GEN_REMOTE}", reasons[row])
                    self.assertFalse(ca.converged()[row])
                    self.assertFalse(ca.converged_mask()[row])
                    viol = ca.get_violations()[row]
                    self.assertEqual(len(viol), 1)
                    self.assertEqual(viol[0].violation_type, LimitViolationType.UNSOLVED)
                    self.assertEqual(viol[0].element_type, ViolationElementType.GRID)
                    self.assertEqual(viol[0].category, ViolationCategory.SOLVER)
                    self.assertTrue(np.all(ca.get_voltages()[row] == 0.))
                else:
                    self.assertEqual(reasons[row], "")
                    self.assertTrue(ca.converged()[row])
                    ref = self._reference(branches=cont[0], gens=cont[3])
                    self.assertGreater(ref.shape[0], 0)
                    np.testing.assert_allclose(ca.get_voltages()[row][self.buses], ref[self.buses],
                                               rtol=1e-8, atol=1e-8)

    def test_compute_does_not_raise(self):
        self._check(self._analysis())

    def test_multithreaded(self):
        self._check(self._analysis(nb_thread=3))

    def test_handle_disconnected_grid(self):
        self._check(self._analysis(handle_disconnected_grid=True))

    def test_flows_of_an_unsolved_row_are_nan(self):
        """as those of a row the connectivity check skips"""
        ca = self._analysis()
        amps = ca.compute_flows()
        row = self._row(ca, gens=[GEN_REMOTE])
        self.assertTrue(np.all(np.isnan(amps[row])))
        other = self._row(ca, gens=[LOCAL_GEN])
        self.assertTrue(np.any(amps[other] != 0.))

    def test_a_split_grid_stays_not_simulated(self):
        """The connectivity check wins: a contingency that strands the remote controller's own
        bus is NOT_SIMULATED, whatever else it trips -- no one-off solve can do better."""
        ca = ContingencyAnalysisCPP(self.grid)
        ca.compute_limit_violations = True
        ca.add_contingency(branch_ids=[self.n_line + TRAFO_BEHIND], gen_ids=[GEN_REMOTE])
        ca.add_n1(0)
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        row = self._row(ca, branches=[self.n_line + TRAFO_BEHIND], gens=[GEN_REMOTE])
        self.assertEqual(ca.get_unsolved_reasons()[row], "")
        self.assertEqual([v.violation_type for v in ca.get_violations()[row]],
                         [LimitViolationType.NOT_SIMULATED])

    def test_reasons_follow_the_contingency_list(self):
        """Registering again (L2 cleared) and recomputing settles the reasons anew."""
        ca = self._analysis()
        self.assertEqual(sum(bool(r) for r in ca.get_unsolved_reasons()), 2)
        ca.remove_contingency(gen_ids=[GEN_REMOTE])
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertEqual(sum(bool(r) for r in ca.get_unsolved_reasons()), 1)
        self.assertEqual(len(ca.get_unsolved_reasons()), len(ca.my_contingencies()))

    def test_without_unsupported_contingency_nothing_changes(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_n1(0)
        ca.add_contingency_gen(LOCAL_GEN)
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertEqual(list(ca.get_unsolved_reasons()), ["", ""])
        self.assertTrue(all(ca.converged_mask()))


if __name__ == "__main__":
    unittest.main()
