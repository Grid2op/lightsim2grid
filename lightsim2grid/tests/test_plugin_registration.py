# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""End-to-end tests for load_algorithm_plugin().

These load a *real* compiled plugin (the examples/external_algorithm/ dummy
solver) and assert that both the happy path and every error path behave: the
error paths must raise a normal, catchable Python exception, and -- the whole
point of the entry-point design -- must NOT abort the interpreter. If a bad
plugin load crashed the process, the test runner would die here rather than
report a failure, so simply reaching the assertions after each error load is
itself the crash-safety check.

The plugin must be built first (see examples/external_algorithm/CMakeLists.txt);
the plugin_registration CI workflow does this. When it is not available the
whole case is skipped so the suite still runs in a plain source checkout.
"""

import os
import sys
import shutil
import tempfile
import unittest


def _find_plugin():
    env = os.environ.get("LS2G_TEST_PLUGIN")
    if env and os.path.exists(env):
        return os.path.abspath(env)
    base = os.path.join(os.path.dirname(__file__), "..", "..",
                        "examples", "external_algorithm", "build")
    if sys.platform == "win32":
        candidates = [
            os.path.join(base, "Release", "dummy_solver.dll"),
            os.path.join(base, "Debug", "dummy_solver.dll"),
            os.path.join(base, "dummy_solver.dll"),
        ]
    else:
        candidates = [os.path.join(base, "libdummy_solver.so")]
    for cand in candidates:
        cand = os.path.abspath(cand)
        if os.path.exists(cand):
            return cand
    return None


class TestPluginRegistration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plugin = _find_plugin()
        if cls.plugin is None:
            raise unittest.SkipTest(
                "example solver plugin not built; build examples/external_algorithm/ "
                "or set LS2G_TEST_PLUGIN to its .so/.dll path")

    def _names(self):
        from lightsim2grid.lightsim2grid_cpp import LSGrid
        return LSGrid().available_algorithm_names()

    def test_1_load_registers_solver(self):
        import lightsim2grid
        from lightsim2grid.lightsim2grid_cpp import LSGrid, AlgorithmType

        if "DummyExternal" not in self._names():
            lightsim2grid.load_algorithm_plugin(self.plugin)

        self.assertIn("DummyExternal", self._names())
        gm = LSGrid()
        gm.change_algorithm("DummyExternal")
        self.assertEqual(gm.get_algo_type(), AlgorithmType.Custom)

    def test_2_duplicate_load_raises_without_crashing(self):
        import lightsim2grid

        if "DummyExternal" not in self._names():
            lightsim2grid.load_algorithm_plugin(self.plugin)

        # Loading the same plugin again re-runs its registration, which now
        # collides with the already-registered name: a catchable error, not a
        # crash (previously this could std::terminate() the interpreter).
        with self.assertRaises(Exception):
            lightsim2grid.load_algorithm_plugin(self.plugin)
        # Still alive, and the registry is intact (atomic rejection).
        self.assertIn("DummyExternal", self._names())

    def test_3_duplicate_from_a_copy_raises_without_crashing(self):
        import lightsim2grid

        if "DummyExternal" not in self._names():
            lightsim2grid.load_algorithm_plugin(self.plugin)

        # A different file on disk that registers the same name must also be
        # rejected cleanly (exercises the collision path via a fresh dlopen).
        suffix = os.path.splitext(self.plugin)[1] or ".so"
        tmp = tempfile.NamedTemporaryFile(prefix="ls2g_dupe_", suffix=suffix, delete=False)
        tmp.close()
        try:
            shutil.copyfile(self.plugin, tmp.name)
            with self.assertRaises(Exception):
                lightsim2grid.load_algorithm_plugin(tmp.name)
            self.assertIn("DummyExternal", self._names())
        finally:
            os.unlink(tmp.name)

    def test_4_library_without_entry_point_raises(self):
        import lightsim2grid
        from lightsim2grid import lightsim2grid_cpp

        # The compiled extension itself is a perfectly loadable shared library,
        # but it exports no ls2g_register_plugin entry point.
        not_a_plugin = lightsim2grid_cpp.__file__
        with self.assertRaises(Exception):
            lightsim2grid.load_algorithm_plugin(not_a_plugin)

    def test_5_missing_file_raises(self):
        import lightsim2grid
        with self.assertRaises(Exception):
            lightsim2grid.load_algorithm_plugin(
                os.path.join(tempfile.gettempdir(), "ls2g_no_such_plugin_xyz.so"))


if __name__ == "__main__":
    unittest.main()
