# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Malformed reactive capability curve data (min_q/max_q effectively "swapped" at
one of the curve points) was found on a real grid snapshot: pypowsybl's
"at target p" interpolation then yields min_q > max_q, which OpenLoadFlow
tolerates silently but which used to make lightsim2grid's
``GeneratorContainer::init`` hard-crash with ``RuntimeError: ... min_q being
above max_q``. `from_pypowsybl.init()` now sorts the pair instead of crashing.
"""

import unittest
import warnings

try:
    import pypowsybl.network as pn
    from lightsim2grid.network import init_from_pypowsybl
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False


class TestGenReactiveCurveSwap(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is required")
        self.net = pn.create_ieee14()
        self.gid = self.net.get_generators().index[0]
        # reproduces the real-world curve found in the wild: the p=0.99 point
        # has min_q (0.053) > max_q (0.04), a data-entry error.
        try:
            self.net.create_curve_reactive_limits(
                id=[self.gid, self.gid], p=[-0.15, 0.99],
                min_q=[0.0, 0.053], max_q=[0.0, 0.04])
        except Exception as exc:  # noqa: BLE001 - pypowsybl's error type is not public
            # pypowsybl 1.17 validates the curve: such a network can no longer be built,
            # nor imported from a file (its `iidm.import.repair-invalid-reactive-curve-limits`
            # sorts the pair on the way in), so the converter cannot meet it any more.
            if "maximum reactive power is expected to be greater" not in str(exc):
                raise
            self.skipTest("this pypowsybl refuses a reactive curve point with min_q > max_q")

    def test_does_not_crash_and_sorts_the_pair(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pypowsybl(self.net, sort_index=False)
        gen = next(g for g in model.get_generators() if g.name == self.gid)
        self.assertLessEqual(gen.min_q_mvar, gen.max_q_mvar)


if __name__ == "__main__":
    unittest.main()
