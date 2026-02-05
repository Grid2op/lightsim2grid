# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import pdb
import sys
import unittest
import warnings
import grid2op
import numpy as np
from packaging import version

from grid2op.tests.helper_path_test import PATH_DATA_TEST_PP, PATH_DATA_TEST

from grid2op.Space import GridObjects  # lazy import
__has_storage = hasattr(GridObjects, "n_storage")

try:
    # new way of doing, does not need to inherit from HelperTests but from unittest.TestCase
    from grid2op._create_test_suite import create_test_suite
    from grid2op.tests.helper_path_test import HelperTests as DEPRECATEDHelper
    
    class _Garbage:
        def setUp(self):
            pass
        
    class _SuperGarbage(DEPRECATEDHelper, _Garbage):
        pass
    
    _garbage = _SuperGarbage()
    _garbage.setUp()
    
    class HelperTests(unittest.TestCase):
        def setUp(self) -> None:
            self.tol_one = _garbage.tol_one
            self.tolvect = _garbage.tolvect
            return super().setUp()
        
        def tearDown(self) -> None:
            return super().tearDown()
    
except ImportError as exc_:
    # old way of doing, need to inherit from that
    from grid2op.tests.helper_path_test import HelperTests  # noqa: F401
    
from grid2op.tests.BaseBackendTest import (BaseTestNames,  # noqa: E402
                                           BaseTestLoadingCase,
                                           BaseTestLoadingBackendFunc, 
                                           BaseTestTopoAction,
                                           BaseTestEnvPerformsCorrectCascadingFailures,
                                           BaseTestChangeBusAffectRightBus,
                                           BaseTestShuntAction,
                                           BaseTestResetEqualsLoadGrid,
                                           BaseTestVoltageOWhenDisco,
                                           BaseTestChangeBusSlack,
                                           BaseIssuesTest,
                                           BaseStatusActions,
                                        #    MakeBackend, AlwaysLegal, CompleteAction, AmbiguousAction
                                           )

from grid2op.tests.test_Environment import (BaseTestLoadingBackendPandaPower,  # noqa: E402
                                            BaseTestResetOk,
                                            BaseTestResetAfterCascadingFailure,
                                            BaseTestCascadingFailure)

if __has_storage:
    from grid2op.tests.BaseBackendTest import BaseTestStorageAction

PATH_DATA_TEST_INIT = PATH_DATA_TEST
PATH_DATA_TEST = PATH_DATA_TEST_PP
from lightsim2grid.lightSimBackend import LightSimBackend
from lightsim2grid.solver import SolverType
from grid2op.Runner import Runner

        
class TestNames(BaseTestNames, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestLoadingCase(BaseTestLoadingCase, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestLoadingBackendFunc(BaseTestLoadingBackendFunc, unittest.TestCase):
    def setUp(self):
        # TODO find something more elegant
        BaseTestLoadingBackendFunc.setUp(self)
        self.tests_skipped = set()

    def tearDown(self):
        # TODO find something more elegant
        BaseTestLoadingBackendFunc.tearDown(self)

    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestTopoAction(BaseTestTopoAction, unittest.TestCase):
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


class TestEnvPerformsCorrectCascadingFailures(BaseTestEnvPerformsCorrectCascadingFailures, unittest.TestCase):
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
    

class TestChangeBusAffectRightBus(BaseTestChangeBusAffectRightBus, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestShuntAction(BaseTestShuntAction, unittest.TestCase):
    # tests_skipped = ["test_shunt_effect"]  # "bug" in grid2op that performs an equality check for shunt_p real values which does not pass (shunt_p = 1e-16)
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk
    

class TestResetEqualsLoadGrid(BaseTestResetEqualsLoadGrid, unittest.TestCase):
    def setUp(self):
        BaseTestResetEqualsLoadGrid.setUp(self)

    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestVoltageOWhenDisco(BaseTestVoltageOWhenDisco, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestChangeBusSlack(BaseTestChangeBusSlack, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestIssuesTest(BaseIssuesTest, unittest.TestCase):
    tests_skipped = ["test_issue_125"]  if version.parse(grid2op.__version__) < version.parse("1.9.2") else []
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestStatusAction(BaseStatusActions, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


if __has_storage:
    class TestStorageAction(BaseTestStorageAction, unittest.TestCase):
        def setUp(self):
            super().setUp()
            self.tests_skipped = ["test_storage_action_topo"]  # TODO this test is super weird ! It's like we impose
            # TODO a behaviour from pandapower (weird one) to all backends...

        def make_backend(self, detailed_infos_for_cascading_failures=False):
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
            return bk


class TestLoadingBackendLightSim(BaseTestLoadingBackendPandaPower, unittest.TestCase):
    def get_backend(self, detailed_infos_for_cascading_failures=True):
        return LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)


class TestResetOkLS(BaseTestResetOk, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        return LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)


class TestResetAfterCascadingFailureLS(BaseTestResetAfterCascadingFailure, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        return LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)


class TestCascadingFailureLS(BaseTestCascadingFailure, unittest.TestCase):
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


class TestBackendArgument(unittest.TestCase):
    def setUp(self) -> None:
        if version.parse(grid2op.__version__) < version.parse("1.7.1"):
            self.skipTest(f"grid2op version too old for the feature. Expecting "
                          f"grid2op >= 1.7.1 found {grid2op.__version__}")
            
    def test_non_default_argument(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make("l2rpn_case14_sandbox",
                               test=True,
                               backend=LightSimBackend(solver_type=SolverType.SparseLUSingleSlack,
                                                       max_iter=100,
                                                       tol=1e-6))
        assert env.backend._grid.get_solver_type() == SolverType.SparseLUSingleSlack
        assert env.backend.max_it == 100
        assert env.backend.tol == 1e-6
        runner = Runner(**env.get_params_for_runner())
        env_cpy = runner.init_env()
        assert env_cpy.backend._grid.get_solver_type() == SolverType.SparseLUSingleSlack
        assert env_cpy.backend.max_it == 100
        assert env_cpy.backend.tol == 1e-6
            
    def test_default_arguments(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make("l2rpn_case14_sandbox",
                               test=True,
                               backend=LightSimBackend())
        if SolverType.KLUSingleSlack in env.backend.available_solvers:
            # in the test KLU is not available on windows
            correct_type = SolverType.KLUSingleSlack
        else:
            correct_type = SolverType.SparseLUSingleSlack
            
        assert env.backend._grid.get_solver_type() == correct_type
        assert env.backend.max_it == 10
        assert env.backend.tol == 1e-8
        runner = Runner(**env.get_params_for_runner())
        env_cpy = runner.init_env()            
        assert env_cpy.backend._grid.get_solver_type() == correct_type
        assert env_cpy.backend.max_it == 10
        assert env_cpy.backend.tol == 1e-8
       
            
if __name__ == "__main__":
    unittest.main()
