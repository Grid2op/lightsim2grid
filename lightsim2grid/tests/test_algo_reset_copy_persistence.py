# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Regression tests for a bug found while auditing the documentation of
`LightSimBackend.set_algo_type` / the `AlgoConfig` mechanism:

  - `set_algo_type` used to reject plain strings, so a string-only built-in
    algorithm (eg `NRRefactorRetry_KLU`, which has no `AlgorithmType` enum value)
    or a plugin could only be selected through the private
    `env.backend._grid.change_algorithm(...)`, which does not survive
    `env.reset()` / `backend.copy()`. `set_algo_type` now accepts a plain string
    too, and persists it across both.
  - there was no public, persistent way to customize the AC/DC `AlgoConfig`
    (scaling / refactor policy parameters): only the private
    `env.backend._grid.set_ac_algo_config(...)` existed, which is silently
    reverted on the next `env.reset()`. `LightSimBackend.set_ac_algo_config` /
    `set_dc_algo_config` are new public, persistent equivalents.
"""

import unittest
import warnings

import grid2op

from lightsim2grid import LightSimBackend
from lightsim2grid.algorithm import AlgorithmType, ScalingPolicyType


class TestAlgoAndConfigPersistence(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox",
                                     test=True,
                                     backend=LightSimBackend())

    def tearDown(self) -> None:
        self.env.close()

    def test_algo_type_survives_reset(self):
        """A string-only built-in algorithm (no AlgorithmType enum value, reported
        as AlgorithmType.Custom) selected through `set_algo_type` must still be
        active after env.reset()."""
        self.env.backend.set_algo_type("NRRefactorRetry_KLU")
        self.assertEqual(self.env.backend._grid.get_algo_type(), AlgorithmType.Custom)

        self.env.reset()

        self.assertEqual(self.env.backend._grid.get_algo_type(), AlgorithmType.Custom)
        # and the grid should still be usable
        obs, reward, done, info = self.env.step(self.env.action_space())
        self.assertFalse(done)

    def test_algo_type_survives_copy(self):
        """Same as above, but through LightSimBackend.copy() (used eg by grid2op's
        runner / multi-process environments)."""
        self.env.backend.set_algo_type("NRRefactorRetry_KLU")

        backend_copy = self.env.backend.copy()

        self.assertEqual(backend_copy._grid.get_algo_type(), AlgorithmType.Custom)

    def test_algo_config_survives_reset(self):
        """A custom AlgoConfig (eg a non-default ScalingPolicyType) set through
        `set_ac_algo_config` must still be active after env.reset()."""
        self.env.backend.set_algo_type(AlgorithmType.NR_SparseLU)
        config = self.env.backend.get_ac_algo_config()
        int_params = list(config.int_params)
        self.assertNotEqual(int_params[0], int(ScalingPolicyType.LineSearch))
        int_params[0] = int(ScalingPolicyType.LineSearch)
        config.int_params = int_params
        self.env.backend.set_ac_algo_config(config)

        self.env.reset()

        after = self.env.backend.get_ac_algo_config()
        self.assertEqual(int(after.int_params[0]), int(ScalingPolicyType.LineSearch))
        # and the grid should still be usable
        obs, reward, done, info = self.env.step(self.env.action_space())
        self.assertFalse(done)

    def test_algo_config_survives_copy(self):
        """Same as above, but through LightSimBackend.copy()."""
        self.env.backend.set_algo_type(AlgorithmType.NR_SparseLU)
        config = self.env.backend.get_ac_algo_config()
        int_params = list(config.int_params)
        int_params[0] = int(ScalingPolicyType.LineSearch)
        config.int_params = int_params
        self.env.backend.set_ac_algo_config(config)

        backend_copy = self.env.backend.copy()

        after = backend_copy.get_ac_algo_config()
        self.assertEqual(int(after.int_params[0]), int(ScalingPolicyType.LineSearch))


if __name__ == "__main__":
    unittest.main()
