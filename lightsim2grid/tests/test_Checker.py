# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import pandas as pd
import unittest
import grid2op
import warnings
import pdb
import numpy as np
from lightsim2grid import LightSimBackend

SparseLUSolver_AVAILBLE = False
try:
    from lightsim2grid_cpp import SparseLUSolver
    SparseLUSolver_AVAILBLE = True
except ImportError:
    # KLU solver is not available, these tests cannot be carried out
    pass


class MakeTestsCase14(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make('rte_case14_realistic',
                                    backend=LightSimBackend(),
                                    test=True,
                                    _add_to_name=f"_{type(self).__name__}")
        self.nb_bus = 14

    def tearDown(self) -> None:
        self.env.close()

    def aux_test_check(self, qlim=False):
        tol = self.env.backend.tol

        # it properly detects when something is working
        mismatch = self.env.backend._grid.check_solution(self.env.backend.V, qlim)
        assert mismatch.shape == (2*self.nb_bus, )
        assert np.all(np.abs(mismatch) <= tol)

        # it detects mismatch in non working cases
        mismatch = self.env.backend._grid.check_solution(2*self.env.backend.V, qlim)
        assert mismatch.shape == (2*self.nb_bus, )
        assert np.any(np.abs(mismatch) > tol)

        mismatch = self.env.backend._grid.check_solution(np.ones(self.env.backend.V.shape[0],
                                                                 dtype=complex), qlim)
        assert mismatch.shape == (2*self.nb_bus, )
        assert np.any(np.abs(mismatch) > tol)

    def test_check_withoutqlims(self):
        self.aux_test_check(False)

    def test_check_withqlims(self):
        with self.assertRaises(AssertionError):
            # this should raise because Q lim are not met
            self.aux_test_check(True)

        # and now i check that there are mismatch only due to qlim
        tol = self.env.backend.tol
        tol_gen_q = 1e-5  # because of the "algo" to compute the gen_q i am using in ls
        qlim = True
        mismatch = self.env.backend._grid.check_solution(self.env.backend.V, qlim)
        _, q_gen_tmp, _ = self.env.backend.generators_info()
        # gen_qmin = self.env.backend.init_pp_backend._grid.gen["min_q_mvar"].values
        # gen_qmax = self.env.backend.init_pp_backend._grid.gen["max_q_mvar"].values
        gen_qmin = self.env.backend.init_pp_backend._grid.gen[["bus", "min_q_mvar"]].groupby("bus").sum().values[:, 0]
        gen_qmax = self.env.backend.init_pp_backend._grid.gen[["bus", "max_q_mvar"]].groupby("bus").sum().values[:, 0]
        gen_bus = self.env.backend.init_pp_backend._grid.gen[["bus", "min_q_mvar"]].groupby("bus").sum().index.values

        dt = pd.DataFrame({"bus":  self.env.gen_to_subid, "gen_q": q_gen_tmp})
        q_gen = dt.groupby("bus").sum().values[:, 0]
        assert np.all(dt.groupby("bus").sum().index.values == gen_bus)
        # check the the stuff that "breaks" the limits are indeed due to issue with qmin / qmax
        assert np.all(np.real(mismatch) <= tol), "some real value are not correct, so the error is not due to qlims"
        mismatch_q = np.imag(mismatch)

        # I lack q when i should have generated too much
        not_enough = q_gen > gen_qmax
        assert np.all(mismatch_q[gen_bus][not_enough] > 0.)
        assert np.all(np.abs(mismatch_q[gen_bus][not_enough] + (gen_qmax[not_enough] - q_gen[not_enough])) <= tol_gen_q)

        # I have too much q when i could not absorb enough
        too_much = q_gen < gen_qmin
        assert np.all(mismatch_q[gen_bus][too_much] < 0.)
        assert np.all(np.abs(mismatch_q[gen_bus][too_much] + (gen_qmin[too_much] - q_gen[too_much])) <= tol_gen_q)
        # I have no mismatch when generator
        gen_ok = (q_gen >= gen_qmin) & (q_gen <= gen_qmax)
        assert np.all(mismatch_q[gen_bus][gen_ok] <= tol)
        # no problem when there are no generators
        assert np.all(mismatch_q[~np.isin(np.arange(mismatch.shape[0]), gen_bus)] <= tol)


class MakeTestsCase118(MakeTestsCase14):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("rte_case118_example", backend=LightSimBackend(), test=True,
                                    _add_to_name=f"_{type(self).__name__}")
        self.nb_bus = type(self.env).n_sub


class MakeTestsNeurips2020Track1(MakeTestsCase14):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_neurips_2020_track1", backend=LightSimBackend(), test=True,
                                    _add_to_name=f"_{type(self).__name__}")
        self.nb_bus = type(self.env).n_sub


class MakeTestsNeurips2020Track2(MakeTestsCase14):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_neurips_2020_track2", backend=LightSimBackend(), test=True,
                                    _add_to_name=f"_{type(self).__name__}")
        self.nb_bus = type(self.env).n_sub


if __name__ == "__main__":
    unittest.main()
