# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""LightSimBackend's `uniform_voltage_start` (LSGrid.set_uniform_voltage_start, #284): set
through the backend or its constructor, it survives env.reset(), backend.copy(), a change
of algorithm, pickling and a backend rebuilt from grid2op's `_my_kwargs`."""

import os
import pickle
import tempfile
import unittest
import warnings

import grid2op

from lightsim2grid import LightSimBackend


def _make(**kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend(**kwargs))


class TestBackendUniformVoltageStart(unittest.TestCase):
    def setUp(self) -> None:
        self.env = _make()

    def tearDown(self) -> None:
        self.env.close()

    def _check(self, backend, value):
        self.assertEqual(backend.get_uniform_voltage_start(), value)
        self.assertEqual(backend._grid.get_uniform_voltage_start(), value)

    def test_off_by_default(self):
        self._check(self.env.backend, False)

    def test_set_applies_now(self):
        self.env.backend.set_uniform_voltage_start(True)
        self._check(self.env.backend, True)
        obs, reward, done, info = self.env.step(self.env.action_space())
        self.assertFalse(done)

    def test_survives_reset(self):
        # the reset grid is rebuilt from the one taken at load time, before the setter ran
        self.env.backend.set_uniform_voltage_start(True)
        self.env.reset()
        self._check(self.env.backend, True)
        obs, reward, done, info = self.env.step(self.env.action_space())
        self.assertFalse(done)
        self.env.backend.set_uniform_voltage_start(False)
        self.env.reset()
        self._check(self.env.backend, False)

    def test_survives_copy(self):
        self.env.backend.set_uniform_voltage_start(True)
        backend_copy = self.env.backend.copy()
        self._check(backend_copy, True)
        backend_copy.reset(self.env.get_path_env())
        self._check(backend_copy, True)

    def test_survives_set_algo_type(self):
        self.env.backend.set_uniform_voltage_start(True)
        self.env.backend.set_algo_type("NR_SparseLU")
        self.env.reset()
        self._check(self.env.backend, True)

    def test_backend_picklable(self):
        self.env.backend.set_uniform_voltage_start(True)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "backend.pickle")
            with open(path, "wb") as f:
                pickle.dump(self.env.backend, f)
            with open(path, "rb") as f:
                backend_1 = pickle.load(f)
        self._check(backend_1, True)
        backend_1.reset(self.env.get_path_env())
        self._check(backend_1, True)

    def test_kwargs_for_rebuilt_backend(self):
        # grid2op rebuilds a backend from `_my_kwargs` when it does not copy it (runner,
        # observation_backend_kwargs): a value set later has to be there too
        self.env.backend.set_uniform_voltage_start(True)
        kwargs = self.env.backend._my_kwargs
        self.assertTrue(kwargs["uniform_voltage_start"])
        self.assertTrue(LightSimBackend(**kwargs).get_uniform_voltage_start())

    def test_same_result_as_without(self):
        # only the path to the solution changes: same observation after a few steps
        ref = _make()
        try:
            self.env.backend.set_uniform_voltage_start(True)
            for env in (self.env, ref):
                env.seed(0)
                env.set_id(0)
            obs, obs_ref = self.env.reset(), ref.reset()
            for _ in range(3):
                obs, *_ = self.env.step(self.env.action_space())
                obs_ref, *_ = ref.step(ref.action_space())
            self.assertLess(abs(obs.v_or - obs_ref.v_or).max(), 1e-4)
            self.assertLess(abs(obs.p_or - obs_ref.p_or).max(), 1e-3)
        finally:
            ref.close()


class TestBackendUniformVoltageStartKwargs(unittest.TestCase):
    def test_constructor(self):
        env = _make(uniform_voltage_start=True)
        try:
            self.assertTrue(env.backend._grid.get_uniform_voltage_start())
            env.reset()
            self.assertTrue(env.backend._grid.get_uniform_voltage_start())
            # the backend of obs.simulate is a copy
            obs = env.reset()
            _, _, done, _ = obs.simulate(env.action_space())
            self.assertFalse(done)
            self.assertTrue(env.backend.copy().get_uniform_voltage_start())
        finally:
            env.close()


if __name__ == "__main__":
    unittest.main()
