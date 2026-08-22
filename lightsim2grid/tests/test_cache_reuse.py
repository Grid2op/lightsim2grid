# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Python-visible half of the cache-reuse tests.

The exhaustive half -- every mutating method of ``LSGrid`` checked against a grid that
rebuilds from scratch -- lives in ``src/tests/test_cache_reuse.cpp`` (tag ``[cache_reuse]``),
where it also runs under the valgrind and ASan/UBSan passes of the C++ CI. What is checked
here is the python-facing API: the getters and setters, the defaults, and that the sequences
that used to segfault the interpreter no longer do.
"""

import os
import pickle
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# dependency-free grid builder (no grid2op / pandapower / pypowsybl needed)
from _exotic_elements_case_ls import build_exotic_elements_case_grid  # noqa: E402


class _Base(unittest.TestCase):
    def _grid(self, reuse=True):
        grid = build_exotic_elements_case_grid()
        grid.allow_cache_reuse(reuse)
        return grid, np.ones(grid.total_bus(), dtype=complex) * grid.get_init_vm_pu()

    def _ac(self, grid, v_init):
        return np.asarray(grid.ac_pf(v_init, 30, 1e-10))

    def _dc(self, grid, v_init):
        return np.asarray(grid.dc_pf(v_init, 1, 1e-10))


class TestCacheReuseApi(_Base):
    def test_on_by_default(self):
        grid, _ = self._grid()
        self.assertTrue(grid.get_allow_cache_reuse())
        self.assertTrue(grid.get_allow_ac_cache_reuse())
        self.assertTrue(grid.get_allow_dc_cache_reuse())

    def test_per_family_setters(self):
        grid, _ = self._grid()
        grid.allow_ac_cache_reuse(False)
        self.assertFalse(grid.get_allow_ac_cache_reuse())
        self.assertTrue(grid.get_allow_dc_cache_reuse())
        self.assertFalse(grid.get_allow_cache_reuse())  # "both" is false as soon as one is

        grid.allow_dc_cache_reuse(False)
        self.assertFalse(grid.get_allow_dc_cache_reuse())

        grid.allow_cache_reuse(True)
        self.assertTrue(grid.get_allow_ac_cache_reuse())
        self.assertTrue(grid.get_allow_dc_cache_reuse())
        self.assertTrue(grid.get_allow_cache_reuse())

    def test_marking_is_automatic_and_per_family(self):
        grid, v_init = self._grid()
        self.assertTrue(grid.get_ac_algo_controler().need_reset_solver())
        self.assertTrue(grid.get_dc_algo_controler().need_reset_solver())

        self.assertGreater(self._ac(grid, v_init).shape[0], 0)
        # no unset_changes() anywhere: the AC family marked itself
        self.assertFalse(grid.get_ac_algo_controler().need_reset_solver())
        self.assertTrue(grid.get_dc_algo_controler().need_reset_solver())

        self.assertGreater(self._dc(grid, v_init).shape[0], 0)
        self.assertFalse(grid.get_dc_algo_controler().need_reset_solver())
        self.assertFalse(grid.get_ac_algo_controler().need_reset_solver())

    def test_prevent_cache_reuse_is_per_family_and_one_shot(self):
        grid, v_init = self._grid()
        self._ac(grid, v_init)
        self._dc(grid, v_init)

        grid.prevent_ac_cache_reuse()
        self.assertTrue(grid.get_ac_algo_controler().need_reset_solver())
        self.assertFalse(grid.get_dc_algo_controler().need_reset_solver())

        grid.prevent_dc_cache_reuse()
        self.assertTrue(grid.get_dc_algo_controler().need_reset_solver())

        # one-shot: the family caches again on the next powerflow
        self._ac(grid, v_init)
        self.assertFalse(grid.get_ac_algo_controler().need_reset_solver())
        self.assertTrue(grid.get_allow_cache_reuse())

        grid.prevent_cache_reuse()
        self.assertTrue(grid.get_ac_algo_controler().need_reset_solver())
        self.assertTrue(grid.get_dc_algo_controler().need_reset_solver())

    def test_tell_solver_need_reset_still_works(self):
        grid, v_init = self._grid()
        self._ac(grid, v_init)
        grid.tell_solver_need_reset()  # historical name of prevent_cache_reuse
        self.assertTrue(grid.get_ac_algo_controler().need_reset_solver())
        self.assertTrue(grid.get_dc_algo_controler().need_reset_solver())

    def test_unset_changes_is_a_noop_when_reuse_is_on(self):
        grid, v_init = self._grid()
        self._ac(grid, v_init)
        grid.prevent_cache_reuse()
        grid.unset_changes()  # both families are automatic: must do nothing
        self.assertTrue(grid.get_ac_algo_controler().need_reset_solver())
        self.assertTrue(grid.get_dc_algo_controler().need_reset_solver())

    def test_per_family_solver_side_getters(self):
        grid, v_init = self._grid()
        self.assertEqual(grid.get_ac_pv_solver().shape[0], 0)
        self.assertEqual(grid.get_dc_pv_solver().shape[0], 0)

        self._ac(grid, v_init)
        self.assertGreater(grid.get_ac_pv_solver().shape[0], 0)
        self.assertGreater(grid.get_ac_slack_weights_solver().shape[0], 0)
        self.assertEqual(grid.get_dc_pv_solver().shape[0], 0)  # untouched by the AC solve

        ac_pv = grid.get_ac_pv_solver().copy()
        ac_pq = grid.get_ac_pq_solver().copy()
        ac_sw = np.array(grid.get_ac_slack_weights_solver(), copy=True)

        self._dc(grid, v_init)
        self.assertGreater(grid.get_dc_pv_solver().shape[0], 0)
        # the DC solve did not write over the AC family's copies
        np.testing.assert_array_equal(grid.get_ac_pv_solver(), ac_pv)
        np.testing.assert_array_equal(grid.get_ac_pq_solver(), ac_pq)
        np.testing.assert_allclose(grid.get_ac_slack_weights_solver(), ac_sw, rtol=0., atol=0.)


class TestCacheReuseIsInvisible(_Base):
    def _sequence(self, grid, v_init, dc=False):
        out = []
        run = self._dc if dc else self._ac
        nb_load = len(grid.get_loads())
        for step in range(8):
            grid.change_p_load(step % nb_load, 10. + step)
            if step == 2:
                grid.deactivate_powerline(4)
            if step == 4:
                grid.change_bus_gen(2, 9)
            if step == 6:
                grid.reactivate_powerline(4)
            out.append(run(grid, v_init))
        return out

    def test_same_answer_with_and_without_reuse(self):
        for dc in (False, True):
            with self.subTest(dc=dc):
                cached, v_init = self._grid(reuse=True)
                rebuilt, _ = self._grid(reuse=False)
                for step, (a, b) in enumerate(zip(self._sequence(cached, v_init, dc),
                                                  self._sequence(rebuilt, v_init, dc))):
                    self.assertEqual(a.shape, b.shape, f"step {step}")
                    np.testing.assert_allclose(a, b, rtol=0., atol=1e-9,
                                               err_msg=f"step {step}")

    def test_dc_magnitudes_follow_the_v_init_given(self):
        """DC solves angles only: the magnitudes it reports are the ones handed to it.

        They used to be cached across calls (guarded by a flag about the GRID's setpoints,
        not about the vector passed in), so a second dc_pf with a different v_init returned
        the first call's magnitudes for every bus not pinned by a controller.
        """
        grid, _ = self._grid()
        v_high = np.ones(grid.total_bus(), dtype=complex) * 1.06
        v_low = np.ones(grid.total_bus(), dtype=complex) * 1.01
        self._dc(grid, v_high)
        got = self._dc(grid, v_low)

        ref_grid, _ = self._grid()
        expected = self._dc(ref_grid, v_low)
        np.testing.assert_allclose(np.abs(got), np.abs(expected), rtol=0., atol=1e-12)


class TestDeserializedGridStartsCold(_Base):
    """Nothing the solvers cache ever crosses a serialization boundary.

    The bus labelling, Ybus / Sbus, the pv-pq split and the slack weights are all derived
    from the elements, so a grid restored from a pickle (or a binary file) rebuilds them
    rather than reading them back. That is a security property, not a performance one: a
    cache read from a file is a second copy of state that cannot be checked against the
    elements it claims to describe.
    """

    def _assert_cold(self, grid):
        self.assertTrue(grid.get_ac_algo_controler().need_reset_solver())
        self.assertTrue(grid.get_dc_algo_controler().need_reset_solver())
        self.assertEqual(grid.get_ac_pv_solver().shape[0], 0)
        self.assertEqual(grid.get_dc_pv_solver().shape[0], 0)
        self.assertEqual(grid.get_ac_pq_solver().shape[0], 0)
        self.assertEqual(grid.get_dc_pq_solver().shape[0], 0)
        self.assertEqual(grid.get_ac_slack_weights_solver().shape[0], 0)
        self.assertEqual(grid.get_dc_slack_weights_solver().shape[0], 0)

    def test_pickle_round_trip_is_cold(self):
        source, v_init = self._grid()
        self.assertGreater(self._ac(source, v_init).shape[0], 0)  # source carries a live cache
        self.assertGreater(self._dc(source, v_init).shape[0], 0)

        restored = pickle.loads(pickle.dumps(source))
        self._assert_cold(restored)

        ref, ref_v = self._grid(reuse=False)
        np.testing.assert_allclose(self._ac(restored, v_init), self._ac(ref, ref_v),
                                   rtol=0., atol=1e-9)

    def test_binary_round_trip_is_cold(self):
        source, v_init = self._grid()
        self.assertGreater(self._ac(source, v_init).shape[0], 0)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "cold.lsb")
            source.save_binary(path)
            restored = type(source).load_binary(path)
        self._assert_cold(restored)

        ref, ref_v = self._grid(reuse=False)
        np.testing.assert_allclose(self._ac(restored, v_init), self._ac(ref, ref_v),
                                   rtol=0., atol=1e-9)


class TestCopyStartsCold(_Base):
    """`LSGrid.copy()` does not carry the cache either: the copy constructor resets the
    solver state, so the copy rebuilds on its first powerflow. Unlike serialization this
    is a choice, not a safety requirement -- a copy is the same grid in the same process,
    so its cache would be valid. It is pinned here so the behaviour is deliberate.
    """

    def test_copy_is_cold_but_correct(self):
        source, v_init = self._grid()
        self.assertGreater(self._ac(source, v_init).shape[0], 0)
        source.change_p_load(1, 61.)

        copied = source.copy()
        self.assertTrue(copied.get_ac_algo_controler().need_reset_solver())
        self.assertEqual(copied.get_ac_pv_solver().shape[0], 0)
        # ... and the settings do travel with it
        self.assertTrue(copied.get_allow_cache_reuse())

        ref, ref_v = self._grid(reuse=False)
        ref.change_p_load(1, 61.)
        np.testing.assert_allclose(self._ac(copied, v_init), self._ac(ref, ref_v),
                                   rtol=0., atol=1e-9)


class TestNoLongerSegfaults(_Base):
    """`unset_changes()` used to be an unchecked assertion: these sequences of documented
    API calls each killed the interpreter with SIGSEGV. It is kept for backward
    compatibility and still marks when automatic reuse is off, so they stay reachable."""

    def test_before_any_powerflow(self):
        grid, v_init = self._grid(reuse=False)
        grid.unset_changes()
        self.assertGreater(self._ac(grid, v_init).shape[0], 0)

    def test_ac_after_dc(self):
        grid, v_init = self._grid(reuse=False)
        self._dc(grid, v_init)
        grid.unset_changes()
        got = self._ac(grid, v_init)
        ref, ref_v = self._grid()
        np.testing.assert_allclose(got, self._ac(ref, ref_v), rtol=0., atol=1e-9)

    def test_dc_after_ac(self):
        grid, v_init = self._grid(reuse=False)
        self._ac(grid, v_init)
        grid.unset_changes()
        got = self._dc(grid, v_init)
        ref, ref_v = self._grid()
        np.testing.assert_allclose(got, self._dc(ref, ref_v), rtol=0., atol=1e-9)

    def test_alternating_families(self):
        grid, v_init = self._grid(reuse=False)
        ref_ac, ref_v = self._grid()
        v_ac_ref = self._ac(ref_ac, ref_v)
        ref_dc, ref_v = self._grid()
        v_dc_ref = self._dc(ref_dc, ref_v)
        for _ in range(3):
            np.testing.assert_allclose(self._ac(grid, v_init), v_ac_ref, rtol=0., atol=1e-9)
            grid.unset_changes()
            np.testing.assert_allclose(self._dc(grid, v_init), v_dc_ref, rtol=0., atol=1e-9)
            grid.unset_changes()


if __name__ == "__main__":
    unittest.main()
