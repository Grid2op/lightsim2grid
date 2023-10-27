# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings

from grid2op.tests.BaseRedispTest import (BaseTestRedispatch,
                                          BaseTestRedispatchChangeNothingEnvironment,
                                          BaseTestRedispTooLowHigh,
                                          BaseTestDispatchRampingIllegalETC,
                                          BaseTestLoadingAcceptAlmostZeroSumRedisp)

from lightsim2grid.lightSimBackend import LightSimBackend


class TestRedispatch(BaseTestRedispatch, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestRedispatchChangeNothingEnvironment(BaseTestRedispatchChangeNothingEnvironment, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestRedispTooLowHigh(BaseTestRedispTooLowHigh, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestDispatchRampingIllegalETC(BaseTestDispatchRampingIllegalETC, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


class TestLoadingAcceptAlmostZeroSumRedisp(BaseTestLoadingAcceptAlmostZeroSumRedisp, unittest.TestCase):
    def make_backend(self, detailed_infos_for_cascading_failures=False):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bk = LightSimBackend(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
        return bk


if __name__ == "__main__":
    unittest.main()
