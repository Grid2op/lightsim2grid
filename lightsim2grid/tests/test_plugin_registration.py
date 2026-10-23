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
from unittest import mock


def _candidates_in(build_dir):
    if sys.platform == "win32":
        names = ["Release/dummy_solver.dll", "Debug/dummy_solver.dll", "dummy_solver.dll"]
    else:
        names = ["libdummy_solver.so"]  # .so on macOS too (see the plugin CMakeLists)
    return [os.path.join(build_dir, n) for n in names]


def _find_plugin():
    # Full path to the built plugin, when known (set by CI).
    env = os.environ.get("LS2G_TEST_PLUGIN")
    if env and os.path.exists(env):
        return os.path.abspath(env)

    build_dirs = []
    env_dir = os.environ.get("LS2G_TEST_PLUGIN_DIR")
    if env_dir:
        build_dirs.append(env_dir)
    # source-tree fallback (in-repo dev run, not the installed package, whose
    # __file__ is in site-packages and cannot reach the source-tree examples/).
    build_dirs.append(os.path.join(os.path.dirname(__file__), "..", "..",
                                   "examples", "external_algorithm", "build"))

    for build_dir in build_dirs:
        for cand in _candidates_in(build_dir):
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

    def test_6_found_by_name_through_env(self):
        import lightsim2grid
        with mock.patch.dict(os.environ, {"LS2G_PLUGIN_PATH": os.path.dirname(self.plugin)}):
            self.assertEqual(os.path.normcase(lightsim2grid.find_algorithm_plugin("dummy_solver")),
                             os.path.normcase(self.plugin))


class TestFindAlgorithmPlugin(unittest.TestCase):
    """find_algorithm_plugin's lookup, on empty files: nothing is loaded."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="ls2g_find_plugin_")
        # no examples/ of a checkout unless a test sets one
        patcher = mock.patch("lightsim2grid._checkout_examples_dir", return_value=None)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def _lib(self, folder, name):
        os.makedirs(folder, exist_ok=True)
        file_name = f"{name}.dll" if sys.platform == "win32" else f"lib{name}.so"
        path = os.path.join(folder, file_name)
        open(path, "wb").close()
        return path

    def _find(self, name, env_dirs=()):
        import lightsim2grid
        with mock.patch.dict(os.environ, {"LS2G_PLUGIN_PATH": os.pathsep.join(env_dirs)}):
            return lightsim2grid.find_algorithm_plugin(name)

    def test_env_dirs_in_order(self):
        first = self._lib(os.path.join(self.tmp, "a"), "my_algorithm")
        self._lib(os.path.join(self.tmp, "b"), "my_algorithm")
        dirs = [os.path.join(self.tmp, "missing"), os.path.join(self.tmp, "a"), os.path.join(self.tmp, "b")]
        self.assertEqual(self._find("my_algorithm", dirs), first)

    def test_algorithm_suffix_optional(self):
        lib = self._lib(self.tmp, "my_algorithm")
        self.assertEqual(self._find("my", [self.tmp]), lib)

    def test_release_subfolder(self):
        lib = self._lib(os.path.join(self.tmp, "Release"), "my_algorithm")
        self.assertEqual(self._find("my_algorithm", [self.tmp]), lib)

    def test_checkout_examples(self):
        examples = os.path.join(self.tmp, "examples")
        for d in ("foo_algorithm", "two_libs_algorithm"):
            os.makedirs(os.path.join(examples, d))
            open(os.path.join(examples, d, "CMakeLists.txt"), "w").close()
        foo = self._lib(os.path.join(examples, "foo_algorithm", "build"), "bar")
        two = self._lib(os.path.join(examples, "two_libs_algorithm", "build"), "two_libs_algorithm")
        self._lib(os.path.join(examples, "two_libs_algorithm", "build"), "other")
        with mock.patch("lightsim2grid._checkout_examples_dir", return_value=examples):
            self.assertEqual(self._find("foo"), foo)  # by folder name, whatever the library is called
            self.assertEqual(self._find("bar"), foo)  # by library name
            self.assertEqual(self._find("two_libs"), two)  # two libraries: the one named after the folder
            # the environment comes first
            env_lib = self._lib(os.path.join(self.tmp, "env"), "bar")
            self.assertEqual(self._find("bar", [os.path.join(self.tmp, "env")]), env_lib)

    def test_not_found_raises(self):
        import lightsim2grid
        with self.assertRaises(RuntimeError) as ctx:
            self._find("ls2g_no_such_plugin_xyz", [self.tmp])
        self.assertIn("ls2g_no_such_plugin_xyz", str(ctx.exception))
        # load_algorithm_plugin raises the same error for a name it cannot find
        with mock.patch.dict(os.environ, {"LS2G_PLUGIN_PATH": self.tmp}):
            with self.assertRaises(RuntimeError):
                lightsim2grid.load_algorithm_plugin("ls2g_no_such_plugin_xyz")


if __name__ == "__main__":
    unittest.main()
