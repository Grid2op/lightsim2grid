# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The pypowsybl converter under pandas' copy-on-write (always on from pandas 3).

With copy-on-write, ``Series.to_numpy`` returns a read-only view of the frame whenever no
conversion is needed, so a converter that modifies such an array in place fails with
"output array is read-only" / "assignment destination is read-only". The networks here are
pypowsybl's built-in ones: they carry a regulating on-load ratio tap changer and hvdc lines,
which is where that happened."""

import contextlib
import unittest
import warnings

import numpy as np
import pandas as pd

try:
    import pypowsybl as pp
    from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl, OlfLoadingParameters
    from lightsim2grid.network.from_pypowsybl._aux_add_storage import _aux_battery_q_limits
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

#: pandas 3 has no way to turn copy-on-write off: there, both loads below run with it
_COW_ALWAYS_ON = int(pd.__version__.split(".")[0]) >= 3


def _copy_on_write(enabled):
    if _COW_ALWAYS_ON:
        return contextlib.nullcontext()
    return pd.option_context("mode.copy_on_write", enabled)


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestPypowsyblCopyOnWrite(unittest.TestCase):
    NETWORKS = ("create_four_substations_node_breaker_network",
                "create_eurostag_tutorial_example1_network",
                "create_ieee9",
                "create_ieee14")

    def _solved_v(self, factory, olf_rules, cow):
        net = getattr(pp.network, factory)()
        with _copy_on_write(cow), warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pypowsybl(net, gen_slack_id=net.get_generators().index[0], sort_index=False,
                                       buses_for_sub=False, olf_rules=olf_rules)
        grid.change_algorithm("NRSing_KLU")
        return grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)

    def test_loads_as_without_copy_on_write(self):
        for factory in self.NETWORKS:
            for olf_rules in (False, OlfLoadingParameters()):
                with self.subTest(network=factory, olf_rules=bool(olf_rules)):
                    v_cow = self._solved_v(factory, olf_rules, cow=True)
                    v_ref = self._solved_v(factory, olf_rules, cow=False)
                    self.assertGreater(v_cow.shape[0], 0)
                    np.testing.assert_array_equal(v_cow, v_ref)

    def test_battery_swapped_q_limits(self):
        # IIDM refuses min_q > max_q on a battery's box, a malformed curve at the target P is
        # where it comes from: the swap is done on the arrays read off the frame
        df_batt = pd.DataFrame({"min_q": [10., -5.], "max_q": [-10., 5.]}, index=["B0", "B1"])
        with _copy_on_write(True):
            min_q, max_q = _aux_battery_q_limits(df_batt)
        np.testing.assert_array_equal(min_q, [-10., -5.])
        np.testing.assert_array_equal(max_q, [10., 5.])
        # the frame itself is untouched
        np.testing.assert_array_equal(df_batt["min_q"].to_numpy(), [10., -5.])


if __name__ == "__main__":
    unittest.main()
