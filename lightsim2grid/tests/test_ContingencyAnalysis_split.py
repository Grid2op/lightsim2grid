# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Tests for the `handle_disconnected_grid` mode of the contingency analysis.

When this mode is ON, a contingency that splits the grid into several connected
components is no longer skipped: the largest component is solved while the buses
of the other component(s) are "masked" (their NR equations are forced to identity
and their voltage reported as 0), without any extra symbolic factorization.
"""

import unittest
import warnings
import numpy as np
import pandapower as pp

from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP
from lightsim2grid.algorithm import AlgorithmType


def _line(net, f, t):
    pp.create_line_from_parameters(net, f, t, length_km=1.,
                                   r_ohm_per_km=0.1, x_ohm_per_km=0.3,
                                   c_nf_per_km=0., max_i_ka=1.)


def _build_radial(n_bus, with_island=True):
    """A simple radial grid 0-1-2(-3). Disconnecting the last line (id 2)
    isolates bus 3. `with_island=False` builds the reference grid (no bus 3)."""
    net = pp.create_empty_network(sn_mva=1.)
    nb = n_bus if with_island else 3
    for _ in range(nb):
        pp.create_bus(net, vn_kv=20.)
    pp.create_ext_grid(net, 0, vm_pu=1.0)
    _line(net, 0, 1)               # line 0
    _line(net, 1, 2)               # line 1
    pp.create_load(net, 1, p_mw=1.0, q_mvar=0.2)
    pp.create_load(net, 2, p_mw=1.0, q_mvar=0.2)
    if with_island:
        _line(net, 2, 3)           # line 2 (splits the grid -> isolates bus 3)
        pp.create_load(net, 3, p_mw=1.0, q_mvar=0.2)
    return net


class TestContingencySplitMode(unittest.TestCase):
    def setUp(self):
        self.max_it = 30
        self.tol = 1e-8
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(_build_radial(4))
            self.ref = init_from_pandapower(_build_radial(4, with_island=False))
        self.V0 = np.ones(self.grid.get_bus_vn_kv().shape[0], dtype=complex)
        # reference voltages for the surviving component {0, 1, 2}
        self.Vref = self.ref.ac_pf(np.ones(self.ref.get_bus_vn_kv().shape[0], dtype=complex),
                                   self.max_it, self.tol)
        assert self.Vref.size == 3, "reference powerflow did not converge"

    def _run(self, handle_disconnected):
        SA = ContingencyAnalysisCPP(self.grid)
        SA.add_n1(2)  # disconnect line 2 -> isolates bus 3
        SA.handle_disconnected_grid = handle_disconnected
        SA.compute(1. * self.V0, self.max_it, self.tol)
        return SA.get_voltages().copy()

    def test_default_is_off(self):
        SA = ContingencyAnalysisCPP(self.grid)
        assert SA.handle_disconnected_grid is False

    def test_flag_off_skips_split(self):
        # legacy behaviour: a splitting contingency is not simulated -> all 0.
        v = self._run(False)
        assert np.max(np.abs(v[0])) == 0., f"split contingency should be skipped, got {v[0]}"

    def test_flag_on_solves_main_component(self):
        v = self._run(True)
        # the surviving component {0, 1, 2} matches the reference powerflow
        assert np.max(np.abs(v[0, :3] - self.Vref)) <= 1e-6, \
            f"main component mismatch: {np.max(np.abs(v[0, :3] - self.Vref)):.2e}"
        # the isolated bus 3 is masked -> reported as exactly 0
        assert v[0, 3] == 0., f"masked bus should be 0, got {v[0, 3]}"

    def test_masked_setter_validates_bool(self):
        SA = ContingencyAnalysisCPP(self.grid)
        # the C++ property only accepts a bool; a clean round-trip is enough here
        SA.handle_disconnected_grid = True
        assert SA.handle_disconnected_grid is True
        SA.handle_disconnected_grid = False
        assert SA.handle_disconnected_grid is False

    def test_error_on_non_nr_ac_algorithm(self):
        # the mode requires an NR algorithm: an AC non-NR solver must be rejected
        SA = ContingencyAnalysisCPP(self.grid)
        SA.change_algorithm(AlgorithmType.GaussSeidel)
        SA.add_n1(2)
        SA.handle_disconnected_grid = True
        with self.assertRaises(RuntimeError):
            SA.compute(1. * self.V0, self.max_it, self.tol)

    def test_dc_flag_on_is_a_noop(self):
        # DC handles connectivity internally: enabling the flag must not raise and
        # must keep the legacy DC behaviour (no exception, computation runs).
        SA = ContingencyAnalysisCPP(self.grid)
        SA.change_algorithm(AlgorithmType.DC_SparseLU)
        SA.add_n1(2)
        SA.handle_disconnected_grid = True
        SA.compute(1. * self.V0, self.max_it, self.tol)  # must not raise

    def test_is_grid_connected_no_segfault(self):
        # regression: is_grid_connected_after_contingency() used to segfault because
        # it relied on the (empty) grid-model Ybus. It must now work both standalone
        # (before compute) and after compute, and agree.
        SA = ContingencyAnalysisCPP(self.grid)
        SA.add_n1(2)  # splits the grid (isolates bus 3)
        SA.add_n1(0)  # also splits (isolates {1, 2, 3} from the slack)
        standalone = np.asarray(SA.is_grid_connected_after_contingency())
        assert standalone.shape == (2,)
        assert np.all(standalone == 0), f"both radial cuts disconnect the grid, got {standalone}"
        SA.compute(1. * self.V0, self.max_it, self.tol)
        after = np.asarray(SA.is_grid_connected_after_contingency())
        assert np.array_equal(standalone, after)


class TestContingencySplitMultiSlack(unittest.TestCase):
    """Two slacks at the two ends of a line; a contingency strands one of them.

    Topology: ext_grid@0 - 1 - 2 - 3 - ext_grid@4. Disconnecting line 2 (2-3)
    splits it into {0,1,2} (with slack 0) and {3,4} (with slack 4). The largest
    component {0,1,2} is solved; the stranded slack 4 has its weight zeroed.
    """
    def setUp(self):
        self.max_it = 30
        self.tol = 1e-8
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(self._net(with_island=True))
            self.ref = init_from_pandapower(self._net(with_island=False))
        self.V0 = np.ones(self.grid.get_bus_vn_kv().shape[0], dtype=complex)
        self.Vref = self.ref.ac_pf(np.ones(self.ref.get_bus_vn_kv().shape[0], dtype=complex),
                                   self.max_it, self.tol)
        assert self.Vref.size == 3

    @staticmethod
    def _net(with_island):
        net = pp.create_empty_network(sn_mva=1.)
        nb = 5 if with_island else 3
        for _ in range(nb):
            pp.create_bus(net, vn_kv=20.)
        pp.create_ext_grid(net, 0, vm_pu=1.0)
        _line(net, 0, 1)   # line 0
        _line(net, 1, 2)   # line 1
        pp.create_load(net, 1, p_mw=1.0, q_mvar=0.2)
        pp.create_load(net, 2, p_mw=1.0, q_mvar=0.2)
        if with_island:
            _line(net, 2, 3)   # line 2 (splits the grid)
            _line(net, 3, 4)   # line 3
            pp.create_ext_grid(net, 4, vm_pu=1.0)   # second slack, stranded by the split
            pp.create_load(net, 3, p_mw=1.0, q_mvar=0.2)
        return net

    def test_stranded_slack_is_masked(self):
        SA = ContingencyAnalysisCPP(self.grid)
        SA.add_n1(2)  # isolates {3, 4} (which holds the second slack)
        SA.handle_disconnected_grid = True
        SA.compute(1. * self.V0, self.max_it, self.tol)
        v = SA.get_voltages()
        # surviving component {0,1,2} solved with slack 0 absorbing everything;
        # the reference grid keeps only slack 0, so the voltages must match.
        assert np.max(np.abs(v[0, :3] - self.Vref)) <= 1e-6, \
            f"main component mismatch: {np.max(np.abs(v[0, :3] - self.Vref)):.2e}"
        # the stranded island {3, 4} is masked -> reported as 0
        assert np.max(np.abs(v[0, 3:5])) == 0., f"masked island should be 0, got {v[0, 3:5]}"


class TestContingencySplitCase14(unittest.TestCase):
    """Integration test on l2rpn_case14_sandbox: enabling the mode must leave the
    connected contingencies bit-identical and additionally solve at least one of
    the contingencies that split the grid (skipped when the mode is off)."""
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})
        self.nb_sub = self.env.n_sub

    def tearDown(self):
        self.env.close()

    def _voltages(self, handle_disconnected):
        SA = ContingencyAnalysisCPP(self.env.backend._grid)
        SA.add_all_n1()
        SA.handle_disconnected_grid = handle_disconnected
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        return SA.get_voltages().copy()

    def test_regression_and_new_contingencies(self):
        voff = self._voltages(False)
        von = self._voltages(True)
        nb_sub = self.nb_sub

        # contingencies solved when the mode is OFF (connected, converged) must be
        # bit-identical when the mode is ON (their masked set is empty)
        off_solved = np.array([np.any(np.abs(voff[i, :nb_sub]) > 1e-9) for i in range(voff.shape[0])])
        assert off_solved.any(), "sanity: some contingencies should be solved with the mode off"
        assert np.max(np.abs(von[off_solved] - voff[off_solved])) == 0., \
            "connected contingencies must be unchanged when the mode is enabled"

        # at least one contingency skipped when OFF is now solved on its largest
        # connected component when ON
        newly_solved = [i for i in range(voff.shape[0])
                        if (not off_solved[i]) and np.any(np.abs(von[i, :nb_sub]) > 1e-9)]
        assert len(newly_solved) >= 1, "the mode should solve at least one split contingency"

        for i in newly_solved:
            vm = np.abs(von[i, :nb_sub])
            live = vm[vm > 1e-9]
            assert np.all(np.isfinite(von[i, :nb_sub])), f"cont {i}: non-finite voltage"
            assert live.min() > 0.5 and live.max() < 1.5, \
                f"cont {i}: unrealistic voltages on the live component [{live.min()}, {live.max()}]"
            assert np.sum(vm <= 1e-9) >= 1, f"cont {i}: expected at least one masked (0) bus"


if __name__ == "__main__":
    unittest.main()
