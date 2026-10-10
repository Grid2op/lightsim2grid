# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The loaders leave a grid on KLU when this build has it, and on SparseLU otherwise."""

import unittest
import warnings

import numpy as np
import pandapower.networks as pn

from lightsim2grid.network import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import ScenarioSweepCPP
from lightsim2grid.algorithm import AlgorithmType


def _load():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pandapower(pn.case118())


class TestLoaderDefaultAlgorithm(unittest.TestCase):
    def test_fastest_available(self):
        grid = _load()
        available = grid.available_default_algorithms()
        expected_ac = AlgorithmType.NR_KLU if AlgorithmType.NR_KLU in available else AlgorithmType.NR_SparseLU
        expected_dc = AlgorithmType.DC_KLU if AlgorithmType.DC_KLU in available else AlgorithmType.DC_SparseLU
        assert grid.get_algo_type() == expected_ac
        assert grid.get_dc_algo_type() == expected_dc
        # a batch class built from it starts from the grid's AC algorithm
        assert ScenarioSweepCPP(grid).get_algo_type() == expected_ac

    def test_same_answer_as_sparselu(self):
        V0 = np.ones(_load().total_bus(), dtype=complex)
        for ac, sparselu in ((True, AlgorithmType.NR_SparseLU), (False, AlgorithmType.DC_SparseLU)):
            default = _load()
            reference = _load()
            reference.change_algorithm(sparselu)
            solve = (lambda g: g.ac_pf(V0, 10, 1e-10)) if ac else (lambda g: g.dc_pf(V0, 10, 1e-10))
            V = solve(default)
            V_ref = solve(reference)
            assert V.shape[0] and V_ref.shape[0]
            assert np.abs(V - V_ref).max() <= 1e-9

    def test_still_selectable(self):
        grid = _load()
        grid.change_algorithm(AlgorithmType.NR_SparseLU)
        assert grid.get_algo_type() == AlgorithmType.NR_SparseLU


if __name__ == "__main__":
    unittest.main()
