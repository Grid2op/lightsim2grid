# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Faithfulness check for the generated, dependency-free exotic-elements case
(`_exotic_elements_case_ls.build_exotic_elements_case_grid`, produced by
`_gen_exotic_elements_case.py`): its converged AC powerflow must match, to (near)
machine precision, the grid built the normal way through `init_from_pypowsybl`
(`_exotic_elements_fixture.build_exotic_elements_grid`). Requires pypowsybl only for
that reference build -- `build_exotic_elements_case_grid` itself does not.
"""

import unittest

import numpy as np

try:
    import pypowsybl  # noqa: F401
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

from lightsim2grid.tests._exotic_elements_case_ls import build_exotic_elements_case_grid


class TestExoticElementsCaseLS(unittest.TestCase):
    def test_case_builds_and_converges_without_pypowsybl(self):
        grid = build_exotic_elements_case_grid()
        v0 = np.ones(grid.total_bus(), dtype=complex)
        conv = grid.ac_pf(v0, 20, 1e-7)
        self.assertNotEqual(len(conv), 0, "the generated case grid did not converge")

    @unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
    def test_case_matches_pypowsybl_build(self):
        from lightsim2grid.tests._exotic_elements_fixture import build_exotic_elements_grid

        grid_ref = build_exotic_elements_grid(solve=True)
        grid_case = build_exotic_elements_case_grid()
        v0 = np.ones(grid_case.total_bus(), dtype=complex)
        conv = grid_case.ac_pf(v0, 20, 1e-7)
        self.assertNotEqual(len(conv), 0, "the generated case grid did not converge")

        np.testing.assert_allclose(grid_case.get_V(), grid_ref.get_V(), rtol=0, atol=1e-10)


if __name__ == "__main__":
    unittest.main()
