# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Python-visible half of the ``unset_changes()`` cache-consistency guard tests
(the C++ half, which also runs under valgrind / ASan, is in
``src/tests/test_lsgrid.cpp``, tag ``[unset_changes]``).

``unset_changes()`` tells the grid "the cached solver-side data matches me", for
BOTH solver families at once. Nothing used to check that claim, so it could be
made about a family whose data had never been built: the next powerflow of that
family reused an empty id map / Ybus / Sbus and indexed them with bus ids in the
hundreds. Every sequence below segfaulted the interpreter -- an out-of-bounds
read, not an exception -- and all of them are documented usage of the API.
"""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# dependency-free grid builder (no grid2op / pandapower / pypowsybl needed)
from _exotic_elements_case_ls import build_exotic_elements_case_grid  # noqa: E402


class TestUnsetChangesCacheGuard(unittest.TestCase):
    def _grid(self):
        grid = build_exotic_elements_case_grid()
        v_init = np.ones(grid.total_bus(), dtype=complex)
        return grid, v_init

    def _ac(self, grid, v_init):
        return np.asarray(grid.ac_pf(v_init, 20, 1e-8))

    def _dc(self, grid, v_init):
        return np.asarray(grid.dc_pf(v_init, 1, 1e-8))

    def test_unset_changes_before_any_powerflow(self):
        """`unset_changes()` on a grid that was never solved."""
        grid, v_init = self._grid()
        grid.unset_changes()  # the claim is simply false: nothing was built yet
        v = self._ac(grid, v_init)
        self.assertGreater(v.shape[0], 0, "ac_pf diverged")

        ref_grid, ref_v_init = self._grid()
        np.testing.assert_allclose(v, self._ac(ref_grid, ref_v_init), rtol=0., atol=1e-10)

    def test_ac_after_dc(self):
        """`unset_changes()` clears the AC flags too, though no AC solve ever ran."""
        grid, v_init = self._grid()
        self.assertGreater(self._dc(grid, v_init).shape[0], 0)
        grid.unset_changes()
        v = self._ac(grid, v_init)
        self.assertGreater(v.shape[0], 0, "ac_pf diverged")

        ref_grid, ref_v_init = self._grid()
        np.testing.assert_allclose(v, self._ac(ref_grid, ref_v_init), rtol=0., atol=1e-10)

    def test_dc_after_ac(self):
        """... and symmetrically for the DC family."""
        grid, v_init = self._grid()
        self.assertGreater(self._ac(grid, v_init).shape[0], 0)
        grid.unset_changes()
        v = self._dc(grid, v_init)
        self.assertGreater(v.shape[0], 0, "dc_pf diverged")

        ref_grid, ref_v_init = self._grid()
        np.testing.assert_allclose(v, self._dc(ref_grid, ref_v_init), rtol=0., atol=1e-10)

    def test_alternating_ac_dc(self):
        """`slack_weights_`, `bus_pv_` and `bus_pq_` are shared by the two families."""
        grid, v_init = self._grid()
        ref_ac_grid, ref_v = self._grid()
        v_ac_ref = self._ac(ref_ac_grid, ref_v)
        ref_dc_grid, ref_v = self._grid()
        v_dc_ref = self._dc(ref_dc_grid, ref_v)

        for _ in range(3):
            np.testing.assert_allclose(self._ac(grid, v_init), v_ac_ref, rtol=0., atol=1e-10)
            grid.unset_changes()
            np.testing.assert_allclose(self._dc(grid, v_init), v_dc_ref, rtol=0., atol=1e-10)
            grid.unset_changes()

    def test_cached_path_is_lossless(self):
        """Reusing the cache must not change the answer (the point of the feature)."""
        cached, v_init = self._grid()
        rebuilt, _ = self._grid()
        nb_load = len(cached.get_loads())
        for step in range(5):
            cached.change_p_load(step % nb_load, 10. + step)
            rebuilt.change_p_load(step % nb_load, 10. + step)
            v_cached = self._ac(cached, v_init)
            v_rebuilt = self._ac(rebuilt, v_init)
            self.assertGreater(v_cached.shape[0], 0, f"ac_pf diverged at step {step}")
            np.testing.assert_allclose(v_cached, v_rebuilt, rtol=0., atol=1e-10)
            cached.unset_changes()  # only `cached` claims its cache is in sync


if __name__ == "__main__":
    unittest.main()
