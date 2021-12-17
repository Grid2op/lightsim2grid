# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings
import grid2op
import numpy as np

from grid2op.tests.helper_path_test import PATH_DATA_TEST_PP, PATH_DATA_TEST

from grid2op.Space import GridObjects  # lazy import
__has_storage = hasattr(GridObjects, "n_storage")

from grid2op.tests.helper_path_test import HelperTests
from grid2op.tests.BaseBackendTest import BaseTestNames, BaseTestLoadingCase, BaseTestLoadingBackendFunc
from grid2op.tests.BaseBackendTest import BaseTestTopoAction, BaseTestEnvPerformsCorrectCascadingFailures
from grid2op.tests.BaseBackendTest import BaseTestChangeBusAffectRightBus, BaseTestShuntAction
from grid2op.tests.BaseBackendTest import BaseTestResetEqualsLoadGrid, BaseTestVoltageOWhenDisco, BaseTestChangeBusSlack
from grid2op.tests.BaseBackendTest import BaseIssuesTest, BaseStatusActions
from grid2op.tests.test_Environment import TestLoadingBackendPandaPower, TestResetOk
from grid2op.tests.test_Environment import TestResetAfterCascadingFailure, TestCascadingFailure

if __has_storage:
    from grid2op.tests.BaseBackendTest import BaseTestStorageAction

PATH_DATA_TEST_INIT = PATH_DATA_TEST
PATH_DATA_TEST = PATH_DATA_TEST_PP
from lightsim2grid.lightSimBackend import LightSimBackend


class TestNames(HelperTests, BaseTestNames):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk

    def get_path(self):
        return PATH_DATA_TEST_INIT


class TestLoadingCase(HelperTests, BaseTestLoadingCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk

    def get_path(self):
        return PATH_DATA_TEST

    def get_casefile(self):
        return "test_case14.json"


class TestLoadingBackendFunc(HelperTests, BaseTestLoadingBackendFunc):
    def setUp(self):
        # TODO find something more elegant
        BaseTestLoadingBackendFunc.setUp(self)
        self.tests_skipped = set()

        # lightsim does not support DC powerflow at the moment
        # self.tests_skipped.add("test_pf_ac_dc")
        # self.tests_skipped.add("test_apply_action_active_value")
        # self.tests_skipped.add("test_runpf_dc")
        # Now (version >= 0.5.5) it does

    def tearDown(self):
        # TODO find something more elegant
        BaseTestLoadingBackendFunc.tearDown(self)

    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk

    def get_path(self):
        return PATH_DATA_TEST

    def get_casefile(self):
        return "test_case14.json"


class TestTopoAction(HelperTests, BaseTestTopoAction):
    def setUp(self):
        BaseTestTopoAction.setUp(self)

    def tearDown(self):
        # TODO find something more elegant
        BaseTestTopoAction.tearDown(self)

    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk

    def get_path(self):
        return PATH_DATA_TEST

    def get_casefile(self):
        return "test_case14.json"


class TestEnvPerformsCorrectCascadingFailures(HelperTests, BaseTestEnvPerformsCorrectCascadingFailures):
    def setUp(self):
        BaseTestEnvPerformsCorrectCascadingFailures.setUp(self)

    def tearDown(self):
        # TODO find something more elegant
        BaseTestEnvPerformsCorrectCascadingFailures.tearDown(self)

    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk

    def get_casefile(self):
        return "test_case14.json"

    def get_path(self):
        return PATH_DATA_TEST


class TestChangeBusAffectRightBus(HelperTests, BaseTestChangeBusAffectRightBus):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestShuntAction(HelperTests, BaseTestShuntAction):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestResetEqualsLoadGrid(HelperTests, BaseTestResetEqualsLoadGrid):
    def setUp(self):
        BaseTestResetEqualsLoadGrid.setUp(self)

    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestVoltageOWhenDisco(HelperTests, BaseTestVoltageOWhenDisco):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestChangeBusSlack(HelperTests, BaseTestChangeBusSlack):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestIssuesTest(HelperTests, BaseIssuesTest):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestStatusAction(HelperTests, BaseStatusActions):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


if __has_storage:
    class TestStorageAction(HelperTests, BaseTestStorageAction):
        def setUp(self):
            self.tests_skipped = ["test_storage_action_topo"]  # TODO this test is super weird ! It's like we impose
            # TODO a behaviour from pandapower (weird one) to all backends...

        def make_backend(self, detailed_infos_for_cascading_failures=False):
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
            return bk


class TestLoadingBackendLightSim(TestLoadingBackendPandaPower):
    def get_backend(self):
        return LightSimBackend()


class TestResetOkLS(TestResetOk):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        return LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)


class TestResetAfterCascadingFailureLS(TestResetAfterCascadingFailure):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        return LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)


class TestCascadingFailureLS(TestCascadingFailure):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        return LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)


class TestTheta(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env_pp = grid2op.make("l2rpn_case14_sandbox", test=True)
            self.env_ls = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())

    def _check_obs(self, obs_ls, obs_pp):
        for attr_nm in ["theta_or", "theta_ex", "load_theta", "gen_theta", "storage_theta"]:
            theta_ls = getattr(obs_ls, attr_nm)
            theta_pp = getattr(obs_pp, attr_nm)
            assert theta_ls.shape == theta_pp.shape
            if theta_ls.shape[0]:
                assert np.max(np.abs(theta_ls - theta_pp)) <= 1e-5, f"error for theta for {attr_nm}"

    def test_theta(self):
        # test in regular circumstances
        obs_ls = self.env_ls.reset()
        obs_pp = self.env_pp.reset()
        self._check_obs(obs_ls, obs_pp)

        # test after a topological change
        sub_id = 4
        act = self.env_pp.action_space({"set_bus": {"substations_id": [(sub_id, [1, 2, 2, 1, 1])]}})
        obs_ls, *_ = self.env_ls.step(act)
        obs_pp, *_ = self.env_pp.step(act)
        self._check_obs(obs_ls, obs_pp)


if __name__ == "__main__":
    unittest.main()
