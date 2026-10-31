# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The result arrays of the batch classes co-own their buffer.

An array returned by get_voltages / get_flows / get_power_flows / compute_flows /
compute_power_flows / get_sbuses keeps the values of the run that produced it, whatever
the batch object does next: another compute(), clear(), or being deleted. Before, they
were views on buffers the next compute() freed (and compute_flows' did not even keep
its object alive).
"""

import gc
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

from lightsim2grid.network import init_from_pandapower
from lightsim2grid.algorithm import AlgorithmType
from lightsim2grid.lightsim2grid_cpp import (TimeSeriesCPP, InjectionSweepCPP, ScenarioSweepCPP,
                                             ContingencyAnalysisCPP, ContinuationSweepCPP)

MAX_IT = 10
TOL = 1e-8


def _grid():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pandapower(pn.case14())


def _churn():
    """Allocate and free, so a freed buffer is likely reused and overwritten."""
    gc.collect()
    junk = [np.random.rand(500) for _ in range(2000)]
    del junk


class _Mixin:
    """For one batch class: build it, and run it with a given load scaling."""
    nb_rows = 3
    second_run_differs = True

    def make(self):
        raise NotImplementedError

    def solve(self, obj, scale):
        raise NotImplementedError

    def setUp(self):
        self.grid = _grid()
        self.v0 = np.ones(self.grid.total_bus(), dtype=complex)
        self.load_p = np.array([el.target_p_mw for el in self.grid.get_loads()])

    def _getters(self, obj):
        obj.compute_flows()
        obj.compute_power_flows()
        return {"get_voltages": obj.get_voltages(), "get_flows": obj.get_flows(),
                "get_power_flows": obj.get_power_flows()}

    def test_arrays_survive_the_next_compute(self):
        obj = self.make()
        self.solve(obj, 1.0)
        held = self._getters(obj)
        copies = {k: v.copy() for k, v in held.items()}
        self.solve(obj, 1.2)
        new = self._getters(obj)
        _churn()
        for k in held:
            np.testing.assert_array_equal(held[k], copies[k], err_msg=k)
        # the next run did change the results: the arrays are snapshots, not views
        if self.second_run_differs:
            a, b = new["get_flows"], copies["get_flows"]
            self.assertFalse(a.shape == b.shape and np.allclose(a, b))

    def test_arrays_survive_clear_and_deletion(self):
        obj = self.make()
        self.solve(obj, 1.0)
        held = self._getters(obj)
        held["compute_flows"] = obj.compute_flows()
        held["compute_power_flows"] = obj.compute_power_flows()
        copies = {k: v.copy() for k, v in held.items()}
        obj.clear()
        del obj
        _churn()
        for k in held:
            np.testing.assert_array_equal(held[k], copies[k], err_msg=k)

    def test_arrays_do_not_hold_the_object(self):
        obj = self.make()
        self.solve(obj, 1.0)
        for arr in self._getters(obj).values():
            self.assertNotIsInstance(arr.base, type(obj))
            self.assertFalse(arr.flags.writeable)
            self.assertEqual(arr.ndim, 2)

    def test_flows_buffer_reused_when_nobody_holds_it(self):
        obj = self.make()
        self.solve(obj, 1.0)
        first = obj.compute_flows()
        addr, vals = first.ctypes.data, first.copy()
        # held: a second call writes elsewhere and leaves the first alone
        second = obj.compute_flows()
        self.assertNotEqual(second.ctypes.data, addr)
        np.testing.assert_array_equal(first, vals)
        np.testing.assert_array_equal(second, vals)
        # dropped: the next call writes in place again, without allocating
        addr = second.ctypes.data
        del first, second
        third = obj.compute_flows()
        self.assertEqual(third.ctypes.data, addr)
        np.testing.assert_array_equal(third, vals)


class _InjectionMixin(_Mixin):
    def solve(self, obj, scale):
        obj.modify_load_p(np.tile(self.load_p * scale, (self.nb_rows, 1)))
        obj.compute(self.v0, MAX_IT, TOL)
        self.assertEqual(obj.get_status(), 1)


class TestTimeSeriesViews(_InjectionMixin, unittest.TestCase):
    def make(self):
        return TimeSeriesCPP(self.grid)

    def test_get_sbuses_survives_the_next_compute(self):
        obj = self.make()
        self.solve(obj, 1.0)
        sb = obj.get_sbuses()
        copy = sb.copy()
        self.assertFalse(sb.flags.writeable)
        self.solve(obj, 1.2)
        new = obj.get_sbuses()
        del obj
        _churn()
        np.testing.assert_array_equal(sb, copy)
        self.assertFalse(np.allclose(new, copy))


class TestInjectionSweepViews(_InjectionMixin, unittest.TestCase):
    def make(self):
        return InjectionSweepCPP(self.grid)


class TestScenarioSweepViews(_InjectionMixin, unittest.TestCase):
    def make(self):
        obj = ScenarioSweepCPP(self.grid)
        mask = np.zeros((self.nb_rows, len(self.grid.get_lines())), dtype=bool)
        mask[1, 0] = True
        obj.set_contingency_lines(mask)
        return obj


class TestScenarioSweepDCViews(TestScenarioSweepViews):
    """The DC fast path rebuilds the voltages lazily, inside get_voltages()."""
    def make(self):
        obj = super().make()
        obj.change_algorithm(AlgorithmType.DC_SparseLU)
        return obj


class TestContingencyAnalysisViews(_Mixin, unittest.TestCase):
    def make(self):
        obj = ContingencyAnalysisCPP(self.grid)
        for l_id in range(3):
            obj.add_n1(l_id)
        return obj

    # it copies the grid, so a second run solves the same cases again
    second_run_differs = False

    def solve(self, obj, scale):
        obj.compute(self.v0, MAX_IT, TOL)


class TestContinuationSweepViews(_Mixin, unittest.TestCase):
    def make(self):
        return ContinuationSweepCPP(self.grid)

    def solve(self, obj, scale):
        obj.set_target_load_p(self.load_p * 2.0 * scale)
        obj.compute(self.v0, MAX_IT, TOL)


if __name__ == "__main__":
    unittest.main()
