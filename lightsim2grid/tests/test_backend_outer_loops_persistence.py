# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The outer-loop configuration of a LightSimBackend (the list of loops, the driver's
parameters, resume_outer_state) set through the backend survives env.reset(),
backend.copy() (so obs.simulate) and pickling, and can be given to the constructor."""

import copy
import math
import os
import pickle
import tempfile
import unittest
import warnings

import grid2op
from grid2op.Exceptions import BackendError

from lightsim2grid import LightSimBackend
from lightsim2grid.algorithm import (DistributedSlack, ReactiveLimits, VoltageMonitoring,
                                     OuterLoopDriverParams, make_outer_loop)

try:
    from lightsim2grid.algorithm import NROuter_KLU  # noqa: F401
    OUTER, OTHER_OUTER = "NROuter_KLU", "NROuter_SparseLU"
except ImportError:
    OUTER, OTHER_OUTER = "NROuter_SparseLU", "NROuter_SparseLU"


def _loops_state(loops):
    # NaN (an unset parameter) compares unequal to itself: compared as None
    return [(loop.name(), list(loop.get_params().int_params),
             [None if math.isnan(el) else el for el in loop.get_params().real_params])
            for loop in loops]


def _my_loops():
    return [DistributedSlack(slack_bus_p_max_mismatch_mw=2.0),
            ReactiveLimits(max_pq_pv_switch=5)]


def _my_driver():
    return OuterLoopDriverParams(max_outer_iterations=30)


class TestOuterLoopObjects(unittest.TestCase):
    """What makes the backend's state possible: loops rebuilt from name + parameters."""
    def test_make_outer_loop_round_trip(self):
        for loop in _my_loops() + [VoltageMonitoring()]:
            params = loop.get_params()
            rebuilt = make_outer_loop(loop.name(), params.int_params, params.real_params)
            self.assertIs(type(rebuilt), type(loop))
            self.assertEqual(_loops_state([rebuilt]), _loops_state([loop]))

    def test_make_outer_loop_unknown_name(self):
        with self.assertRaises(Exception):
            make_outer_loop("NotALoop")

    def test_loops_pickle_and_deepcopy(self):
        loops = _my_loops()
        for other in (pickle.loads(pickle.dumps(loops)), copy.deepcopy(loops)):
            self.assertEqual([type(el) for el in other], [type(el) for el in loops])
            self.assertEqual(_loops_state(other), _loops_state(loops))
            self.assertIsNot(other[0], loops[0])

    def test_driver_params_pickle_and_deepcopy(self):
        driver = OuterLoopDriverParams(max_outer_iterations=30, voltage_remote_control_robust_mode=False,
                                       min_realistic_voltage=0.6, max_realistic_voltage=1.4,
                                       min_nominal_voltage_realistic_check=10.)
        self.assertEqual(pickle.loads(pickle.dumps(driver)), driver)
        self.assertEqual(copy.deepcopy(driver), driver)


class TestBackendOuterLoopsPersistence(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True,
                                    backend=LightSimBackend(algo_type=OUTER))
        self.env.backend.set_outer_loops(_my_loops())
        self.env.backend.set_outer_loop_driver_params(_my_driver())

    def tearDown(self) -> None:
        self.env.close()

    def _check(self, backend):
        self.assertEqual(_loops_state(backend._grid.get_outer_loops()), _loops_state(_my_loops()))
        self.assertEqual(backend._grid.get_outer_loop_driver_params(), _my_driver())
        self.assertEqual(_loops_state(backend.get_outer_loops()), _loops_state(_my_loops()))
        self.assertEqual(backend.get_outer_loop_driver_params(), _my_driver())

    def test_set_applies_now(self):
        self._check(self.env.backend)
        obs, reward, done, info = self.env.step(self.env.action_space())
        self.assertFalse(done)

    def test_survives_reset(self):
        self.env.reset()
        self._check(self.env.backend)
        obs, reward, done, info = self.env.step(self.env.action_space())
        self.assertFalse(done)

    def test_survives_copy(self):
        backend_copy = self.env.backend.copy()
        self._check(backend_copy)
        # the copy's own reset keeps it too (its grid is rebuilt from the initial one)
        backend_copy.reset(self.env.get_path_env())
        self._check(backend_copy)

    def test_survives_set_algo_type(self):
        self.env.backend.set_algo_type(OTHER_OUTER)
        self.env.reset()
        self._check(self.env.backend)

    def test_empty_list(self):
        self.env.backend.set_outer_loops([])
        self.env.reset()
        self.assertEqual(self.env.backend._grid.get_outer_loops(), [])

    def test_reset_outer_loops(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ref = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend(algo_type=OUTER))
        default_loops = _loops_state(ref.backend._grid.get_outer_loops())
        ref.close()
        self.env.backend.reset_outer_loops()
        self.assertEqual(_loops_state(self.env.backend._grid.get_outer_loops()), default_loops)
        self.env.reset()
        self.assertEqual(_loops_state(self.env.backend._grid.get_outer_loops()), default_loops)
        # the driver's parameters are not touched
        self.assertEqual(self.env.backend._grid.get_outer_loop_driver_params(), _my_driver())

    def test_kwargs_for_rebuilt_backend(self):
        """grid2op rebuilds a backend from `_my_kwargs` when it does not copy it (runner,
        observation_backend_kwargs): the configuration set later must be there too."""
        kwargs = self.env.backend._my_kwargs
        self.assertEqual(_loops_state(kwargs["outer_loops"]), _loops_state(_my_loops()))
        self.assertEqual(kwargs["outer_loop_driver_params"], _my_driver())
        rebuilt = LightSimBackend(**kwargs)
        self.assertEqual(_loops_state(rebuilt.get_outer_loops()), _loops_state(_my_loops()))
        self.assertEqual(rebuilt.get_outer_loop_driver_params(), _my_driver())

    def test_backend_picklable(self):
        self.env.backend.set_resume_outer_state(True)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "backend.pickle")
            with open(path, "wb") as f:
                pickle.dump(self.env.backend, f)
            with open(path, "rb") as f:
                backend_1 = pickle.load(f)
        self._check(backend_1)
        backend_1.reset(self.env.get_path_env())
        self._check(backend_1)
        self.assertTrue(backend_1.get_resume_outer_state())

    def test_resume_outer_state(self):
        self.assertFalse(self.env.backend.get_resume_outer_state())
        self.env.backend.set_resume_outer_state(True)
        self.assertTrue(self.env.backend.get_resume_outer_state())
        self.env.reset()
        self.assertTrue(self.env.backend.get_resume_outer_state())
        self.assertTrue(self.env.backend.copy().get_resume_outer_state())
        self.env.backend.set_algo_type(OTHER_OUTER)
        self.assertTrue(self.env.backend.get_resume_outer_state())
        obs, reward, done, info = self.env.step(self.env.action_space())
        self.assertFalse(done)
        self.env.backend.set_resume_outer_state(False)
        self.env.reset()
        self.assertFalse(self.env.backend.get_resume_outer_state())

    def test_invalid_inputs(self):
        with self.assertRaises(BackendError):
            self.env.backend.set_outer_loops(ReactiveLimits())  # not a list
        with self.assertRaises(BackendError):
            self.env.backend.set_outer_loops([ReactiveLimits(), ReactiveLimits()])
        with self.assertRaises(BackendError):
            self.env.backend.set_outer_loops(["ReactiveLimits"])
        with self.assertRaises(BackendError):
            self.env.backend.set_outer_loop_driver_params(30)
        # nothing was changed by the refused calls
        self.env.reset()
        self._check(self.env.backend)


class TestBackendOuterLoopsKwargs(unittest.TestCase):
    """Given to the constructor: in place before the first powerflow, and on every copy
    grid2op makes (the backend of obs.simulate)."""
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True,
                                    backend=LightSimBackend(algo_type=OUTER,
                                                            outer_loops=_my_loops(),
                                                            outer_loop_driver_params=_my_driver(),
                                                            resume_outer_state=True))

    def tearDown(self) -> None:
        self.env.close()

    def _check(self, backend):
        self.assertEqual(_loops_state(backend._grid.get_outer_loops()), _loops_state(_my_loops()))
        self.assertEqual(backend._grid.get_outer_loop_driver_params(), _my_driver())
        self.assertTrue(backend.get_resume_outer_state())
        self.assertEqual(backend._grid.get_ac_algo_config().int_params[5], 1)

    def test_after_make(self):
        self._check(self.env.backend)

    def test_after_reset(self):
        self.env.reset()
        self._check(self.env.backend)

    def test_simulate(self):
        obs = self.env.reset()
        sim_obs, sim_r, sim_done, sim_info = obs.simulate(self.env.action_space())
        self.assertFalse(sim_done)
        self._check(self.env.observation_space.obs_env.backend)

    def test_set_before_load_grid(self):
        backend = LightSimBackend(algo_type=OUTER)
        self.assertIsNone(backend.get_outer_loops())
        backend.set_outer_loops(_my_loops())
        backend.set_outer_loop_driver_params(_my_driver())
        backend.set_resume_outer_state(True)
        self.assertEqual(_loops_state(backend.get_outer_loops()), _loops_state(_my_loops()))
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=backend)
        try:
            self._check(env.backend)
            env.reset()
            self._check(env.backend)
        finally:
            env.close()


if __name__ == "__main__":
    unittest.main()
