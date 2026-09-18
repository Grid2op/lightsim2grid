# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Tests for the *distributed slack* in DC power flow.

lightsim2grid distributes the global active power imbalance ``-sum(Pbus)`` across the
participating slack buses proportionally to their (normalized) ``slack_weight`` - the
same convention as the AC distributed slack. This is deliberately different from
pandapower's DC power flow, which treats every slack bus as an independent angle
reference (and ignores the weights). The tests below therefore validate the physical
invariants of the weight-based convention rather than equality with pandapower's DC.
"""

import unittest
import warnings

import numpy as np
import pandapower as pp
import pandapower.networks as pn

from lightsim2grid.gridmodel import init_from_pandapower
from test_multiple_slack import make_grid_multiple_slack, int_or_null


class TestDCDistributedSlack(unittest.TestCase):
    def setUp(self) -> None:
        if [int_or_null(el) for el in pp.__version__.split(".")] < [2, 8, 0]:
            self.skipTest("installed pandapower version does not support multiple slack")

        self.max_it = 10
        self.tol = 1e-8
        self.nb_bus_total = 14
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.case = pn.case14()
            self.net = make_grid_multiple_slack(self.case)
        # id of the (single) reference slack generator added by make_grid_multiple_slack
        self.id_ref_slack = self.net.gen.shape[0] - 1
        if "slack_weight" not in self.net.gen:
            self.net.gen["slack_weight"] = 0.
        self.net.gen["slack_weight"][[self.id_ref_slack]] = 1.

    def _build_and_solve(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            gridmodel = init_from_pandapower(self.net)
        setpoints = np.array([g.target_p_mw for g in gridmodel.get_generators()])
        V = gridmodel.dc_pf(np.ones(self.nb_bus_total, dtype=complex), self.max_it, self.tol)
        assert len(V), "lightsim2grid DC powerflow diverged"
        gen_p = np.array([g.res_p_mw for g in gridmodel.get_generators()])
        load_p = sum(l.res_p_mw for l in gridmodel.get_loads())
        return gridmodel, V, setpoints, gen_p, load_p

    def _assert_balanced(self, gen_p, load_p):
        # DC is lossless: total generation == total load (case14 has no sgen / shunt P here)
        assert abs(gen_p.sum() - load_p) <= 1e-6, f"power not balanced: {gen_p.sum() - load_p}"

    def test_single_slack_matches_pandapower(self):
        """with a single slack, the slack generator absorbs the whole imbalance, as pandapower does.

        (DC voltage angles vs pandapower are already validated on the reference grids in
        test_DCPF.py; here we focus on the single-slack power balance / slack power.)
        """
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.rundcpp(self.net)
        gridmodel, V, setpoints, gen_p, load_p = self._build_and_solve()
        self._assert_balanced(gen_p, load_p)
        assert np.max(np.abs(gen_p - self.net.res_gen["p_mw"].values)) <= 1e-6, "gen P differs from pandapower"

    def test_two_slacks_equal_weights(self):
        """two slacks with equal weights: balanced, and the imbalance is split 50/50."""
        gen_added = 1
        self.net.gen["slack"][[gen_added]] = True
        self.net.gen["slack_weight"][[gen_added, self.id_ref_slack]] = 0.5, 0.5
        gridmodel, V, setpoints, gen_p, load_p = self._build_and_solve()
        self._assert_balanced(gen_p, load_p)
        # each slack absorbs an equal share of the global imbalance
        adj = gen_p - setpoints
        assert abs(adj[gen_added] - adj[self.id_ref_slack]) <= 1e-6, f"imbalance not split equally: {adj}"
        # the non-slack generators are not affected
        non_slack = [i for i in range(len(gen_p)) if i not in (gen_added, self.id_ref_slack)]
        assert np.max(np.abs(adj[non_slack])) <= 1e-6, "non-slack generators should not absorb the imbalance"

    def test_two_slacks_weight_proportional(self):
        """two slacks with different weights: the adjustment is proportional to the weight."""
        gen_added = 1
        w_added, w_ref = 0.75, 0.25
        self.net.gen["slack"][[gen_added]] = True
        self.net.gen["slack_weight"][[gen_added, self.id_ref_slack]] = w_added, w_ref
        gridmodel, V, setpoints, gen_p, load_p = self._build_and_solve()
        self._assert_balanced(gen_p, load_p)
        adj = gen_p - setpoints
        # adjustment_i / weight_i must be identical for every participating slack bus
        assert abs(adj[gen_added] / w_added - adj[self.id_ref_slack] / w_ref) <= 1e-6, \
            f"slack adjustment not proportional to weight: {adj}"
        # and the shares sum to the full imbalance
        assert abs((adj[gen_added] + adj[self.id_ref_slack]) - (load_p - setpoints.sum())) <= 1e-6

    def test_weights_change_the_split_not_the_total(self):
        """changing the weights redistributes the imbalance but keeps the total slack power."""
        gen_added = 1
        self.net.gen["slack"][[gen_added]] = True

        self.net.gen["slack_weight"][[gen_added, self.id_ref_slack]] = 0.5, 0.5
        _, _, setp_a, gen_p_a, load_a = self._build_and_solve()

        self.net.gen["slack_weight"][[gen_added, self.id_ref_slack]] = 0.9, 0.1
        _, _, setp_b, gen_p_b, load_b = self._build_and_solve()

        # total slack power identical, individual split different
        total_a = gen_p_a[gen_added] + gen_p_a[self.id_ref_slack]
        total_b = gen_p_b[gen_added] + gen_p_b[self.id_ref_slack]
        assert abs(total_a - total_b) <= 1e-6, "total slack power should not depend on the weights"
        assert abs(gen_p_a[gen_added] - gen_p_b[gen_added]) > 1e-3, "the split should depend on the weights"


if __name__ == "__main__":
    unittest.main()
