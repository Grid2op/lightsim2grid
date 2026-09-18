# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings

import numpy as np

import grid2op
from grid2op.Parameters import Parameters
from lightsim2grid import LightSimBackend

class TestEnvDN(unittest.TestCase):
    """
    This class tests the new functionality introduced in version 0.5.5, where the Ybus matrix
    is not recomputed when the topology does not change. This is an optimization that can
    drastically improve the performances in such cases. 
    """
    def _aux_test_do_nothing(self, is_dc):
        class LightSimBackend_Test(LightSimBackend):
            def __init__(self, detailed_infos_for_cascading_failures=False):
                super().__init__(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
                self.is_dc = False
            def runpf(self, is_dc=False):
                self.is_dc = is_dc
                return super().runpf(is_dc=is_dc)

        param = Parameters()
        param.ENV_DC = is_dc
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make("l2rpn_case14_sandbox",
                               test=True, 
                               param=param,
                               backend=LightSimBackend_Test())
        obs = env.reset()
        # print("\n\n\n")
        # print("##################")
        # print("env.reset() called")
        obs, *_ = env.step(env.action_space())
        # print("\n\n\n")
        # print("###################")
        # print("1st env.step called")
        assert env.backend.is_dc == is_dc
        obs, *_ = env.step(env.action_space())
        # print("\n\n\n")
        # print("###################")
        # print("2nd env.step called")
        assert env.backend.is_dc == is_dc
        obs, *_ = env.step(env.action_space())
        # print("\n\n\n")
        # print("###################")
        # print("3rd env.step called")
        assert env.backend.is_dc == is_dc
        obs, *_ = env.step(env.action_space())
        # print("\n\n\n")
        # print("###################")
        # print("4th env.step called")
        assert env.backend.is_dc == is_dc
    
    def test_dc(self):
        self._aux_test_do_nothing(is_dc=True)
    
    def test_ac(self):
        self._aux_test_do_nothing(is_dc=False)


class TestDcThenAcSameBackend(unittest.TestCase):
    """A DC powerflow followed by an AC one, on the same backend instance.

    This used to be unsafe: the two solver families shared their cached bus
    labelling / matrices / pv-pq split, so an AC powerflow run after a DC one
    could pick up the DC family's data -- an out-of-bounds read, ie a segfault.
    `LightSimBackend.runpf` carried a `self._last_dc` /
    `tell_solver_need_reset()` workaround for it ("otherwise might segfault").
    Since 1.0.0 the families cache independently and the workaround is gone;
    this pins the behaviour it used to protect.
    """
    def _make_backend(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        env.reset()
        return env

    def _ac_flows(self, backend):
        conv, exc_ = backend.runpf(is_dc=False)
        assert conv, f"AC powerflow diverged: {exc_}"
        p_or, q_or, v_or, a_or = backend.lines_or_info()
        return np.concatenate((p_or, q_or, v_or, a_or))

    def _dc_flows(self, backend):
        conv, exc_ = backend.runpf(is_dc=True)
        assert conv, f"DC powerflow diverged: {exc_}"
        p_or, q_or, v_or, a_or = backend.lines_or_info()
        return np.concatenate((p_or, q_or, v_or, a_or))

    def test_ac_after_dc(self):
        """the AC result must not depend on a DC powerflow having run before it"""
        ref_env = self._make_backend()
        try:
            ref = self._ac_flows(ref_env.backend)
        finally:
            ref_env.close()

        env = self._make_backend()
        try:
            self._dc_flows(env.backend)
            np.testing.assert_allclose(self._ac_flows(env.backend), ref, rtol=0., atol=1e-9)
        finally:
            env.close()

    def test_dc_after_ac(self):
        """... and symmetrically for the DC result"""
        ref_env = self._make_backend()
        try:
            ref = self._dc_flows(ref_env.backend)
        finally:
            ref_env.close()

        env = self._make_backend()
        try:
            self._ac_flows(env.backend)
            np.testing.assert_allclose(self._dc_flows(env.backend), ref, rtol=0., atol=1e-9)
        finally:
            env.close()

    def test_alternating(self):
        """alternating the two families on one backend, several times over"""
        ref_env = self._make_backend()
        try:
            ac_ref = self._ac_flows(ref_env.backend)
        finally:
            ref_env.close()
        ref_env = self._make_backend()
        try:
            dc_ref = self._dc_flows(ref_env.backend)
        finally:
            ref_env.close()

        env = self._make_backend()
        try:
            for _ in range(3):
                np.testing.assert_allclose(self._ac_flows(env.backend), ac_ref, rtol=0., atol=1e-9)
                np.testing.assert_allclose(self._dc_flows(env.backend), dc_ref, rtol=0., atol=1e-9)
        finally:
            env.close()


if __name__ == "__main__":
    unittest.main()
