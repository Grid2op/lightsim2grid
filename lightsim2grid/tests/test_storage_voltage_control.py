# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A storage unit regulating the voltage of its own bus (``init_storages_full``):
it is a PV bus exactly like a local voltage-regulating generator would be, its
reactive output is solved for (within its limits for the residual split), the
setting survives pickle / binary round trips, and the pypowsybl converter reads it
off an IIDM battery's ``voltageRegulation`` extension."""

import os
import pickle
import tempfile
import unittest
import warnings

import numpy as np

from lightsim2grid.network import init_from_pandapower, LSGrid

try:
    import pandapower as pp
    import pandapower.networks as pn
    PP_AVAILABLE = True
except ImportError:
    PP_AVAILABLE = False

try:
    import pypowsybl as pypo
    import pypowsybl.loadflow as pypo_lf
    from lightsim2grid.network import init_from_pypowsybl
    PYPO_AVAILABLE = True
except ImportError:
    PYPO_AVAILABLE = False


BUS = 9          # pandapower bus of case14 hosting the unit under test (a load bus, no generator)
P_MW = 12.0      # load convention: the unit charges 12 MW
VM = 1.035
MIN_Q, MAX_Q = -50.0, 50.0
TOL = 1e-10


def _solve(grid):
    n_bus = grid.get_bus_vn_kv().shape[0]
    v0 = grid.dc_pf(np.ones(n_bus, dtype=complex), 1, 1e-6)
    V = grid.ac_pf(v0.copy(), 30, TOL)
    assert V.shape[0] > 0, "powerflow did not converge"
    return V


@unittest.skipIf(not PP_AVAILABLE, "pandapower is not installed")
class TestStorageVoltageControl(unittest.TestCase):
    def _grid_with_storage(self, regulating=True):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grid = init_from_pandapower(pn.case14())
        grid.init_storages_full(np.array([P_MW]), np.array([0.0]),
                                [bool(regulating)], np.array([VM]),
                                np.array([MIN_Q]), np.array([MAX_Q]),
                                np.array([BUS], dtype=np.int32))
        grid.tell_solver_need_reset()
        return grid

    def _grid_with_generator(self):
        """The same grid with a generator playing the storage unit's role."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            net = pn.case14()
            pp.create_gen(net, bus=BUS, p_mw=-P_MW, vm_pu=VM, min_q_mvar=MIN_Q, max_q_mvar=MAX_Q)
            return init_from_pandapower(net)

    def test_info(self):
        grid = self._grid_with_storage()
        sto = grid.get_storages()
        self.assertEqual(len(sto), 1)
        self.assertTrue(sto[0].voltage_regulator_on)
        self.assertAlmostEqual(sto[0].target_vm_pu, VM)
        self.assertAlmostEqual(sto[0].min_q_mvar, MIN_Q)
        self.assertAlmostEqual(sto[0].max_q_mvar, MAX_Q)
        self.assertEqual(sto[0].regulated_bus_id, BUS)
        self.assertEqual(sto[0].target_p_mw, P_MW)
        # the plain init keeps every unit PQ
        grid.init_storages(np.array([P_MW]), np.array([1.0]), np.array([BUS], dtype=np.int32))
        self.assertFalse(grid.get_storages()[0].voltage_regulator_on)

    def test_regulating_storage_is_pv(self):
        grid = self._grid_with_storage()
        V = _solve(grid)
        bus_solver = int(grid.id_me_to_ac_solver()[BUS])
        self.assertIn(bus_solver, list(grid.get_ac_pv_solver()))
        self.assertAlmostEqual(abs(V[BUS]), VM, places=8)
        sto = grid.get_storages()[0]
        self.assertAlmostEqual(sto.res_p_mw, P_MW)
        # the reactive output is solved for: NOT the (zero) setpoint, and inside the range
        self.assertGreater(abs(sto.res_q_mvar), 1e-3)
        self.assertLessEqual(-sto.res_q_mvar, MAX_Q + 1e-6)     # load convention -> generator
        self.assertGreaterEqual(-sto.res_q_mvar, MIN_Q - 1e-6)

    def test_matches_equivalent_generator(self):
        V_sto = _solve(self._grid_with_storage())
        grid_gen = self._grid_with_generator()
        V_gen = _solve(grid_gen)
        np.testing.assert_allclose(V_sto, V_gen, atol=1e-8, rtol=0)
        gen = [g for g in grid_gen.get_generators() if g.bus_id == BUS][0]
        sto = self._grid_with_storage()
        _solve(sto)
        # same reactive output, opposite conventions
        self.assertAlmostEqual(sto.get_storages()[0].res_q_mvar, -gen.res_q_mvar, places=6)

    def test_non_regulating_storage_is_pq(self):
        grid = self._grid_with_storage(regulating=False)
        V = _solve(grid)
        bus_solver = int(grid.id_me_to_ac_solver()[BUS])
        self.assertNotIn(bus_solver, list(grid.get_ac_pv_solver()))
        self.assertNotAlmostEqual(abs(V[BUS]), VM, places=4)
        self.assertAlmostEqual(grid.get_storages()[0].res_q_mvar, 0.0)

    def test_change_v_storage(self):
        grid = self._grid_with_storage()
        _solve(grid)
        grid.change_v_storage(0, 1.02)
        V = _solve(grid)
        self.assertAlmostEqual(abs(V[BUS]), 1.02, places=8)
        self.assertAlmostEqual(grid.get_storages()[0].target_vm_pu, 1.02)

    def test_pickle_and_binary_round_trip(self):
        grid = self._grid_with_storage()
        V_ref = _solve(grid)
        for other in (pickle.loads(pickle.dumps(grid)), grid.copy()):
            sto = other.get_storages()[0]
            self.assertTrue(sto.voltage_regulator_on)
            self.assertAlmostEqual(sto.target_vm_pu, VM)
            self.assertAlmostEqual(sto.max_q_mvar, MAX_Q)
            np.testing.assert_allclose(_solve(other), V_ref, atol=1e-10, rtol=0)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "grid.l2gb")
            grid.save_binary(path)
            other = LSGrid.load_binary(path)
        sto = other.get_storages()[0]
        self.assertTrue(sto.voltage_regulator_on)
        self.assertAlmostEqual(sto.target_vm_pu, VM)
        self.assertAlmostEqual(sto.min_q_mvar, MIN_Q)
        np.testing.assert_allclose(_solve(other), V_ref, atol=1e-10, rtol=0)


@unittest.skipIf(not PYPO_AVAILABLE, "pypowsybl is not installed")
class TestBatteryVoltageRegulationPypowsybl(unittest.TestCase):
    """An IIDM battery with a ``voltageRegulation`` extension becomes a
    voltage-regulating storage unit, and the solved voltages match OpenLoadFlow."""

    def _network(self, regulating=True):
        net = pypo.network.create_ieee14()
        vls = net.get_voltage_levels()
        vl = vls.index[4]
        bus = net.get_bus_breaker_topology(vl).buses.index[0]
        net.create_batteries(id="BAT1", voltage_level_id=vl, bus_id=bus, min_p=-20., max_p=20.,
                             target_p=3.0, target_q=0.0)
        net.create_minmax_reactive_limits(id="BAT1", min_q=-100., max_q=100.)
        net.create_extensions("voltageRegulation", id="BAT1", voltage_regulator_on=bool(regulating),
                              target_v=VM * vls.at[vl, "nominal_v"])
        return net

    def test_converter_reads_the_extension(self):
        net = self._network()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grid = init_from_pypowsybl(net, sort_index=False, buses_for_sub=False)
        sto = grid.get_storages()
        self.assertEqual(len(sto), 1)
        self.assertTrue(sto[0].voltage_regulator_on)
        self.assertAlmostEqual(sto[0].target_vm_pu, VM, places=9)
        self.assertAlmostEqual(sto[0].target_p_mw, -3.0)        # IIDM generator -> load convention
        self.assertAlmostEqual(sto[0].min_q_mvar, -100.)
        self.assertAlmostEqual(sto[0].max_q_mvar, 100.)
        # and off, a plain PQ unit
        net_off = self._network(regulating=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grid_off = init_from_pypowsybl(net_off, sort_index=False, buses_for_sub=False)
        self.assertFalse(grid_off.get_storages()[0].voltage_regulator_on)

    def test_matches_openloadflow(self):
        net = self._network()
        res = pypo_lf.run_ac(net, pypo_lf.Parameters(distributed_slack=False, use_reactive_limits=False))
        self.assertEqual(res[0].status, pypo_lf.ComponentStatus.CONVERGED)
        buses = net.get_buses()
        nominal = net.get_voltage_levels()["nominal_v"].reindex(buses["voltage_level_id"]).to_numpy()
        vm_olf = buses["v_mag"].to_numpy() / nominal
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grid = init_from_pypowsybl(net, sort_index=False, buses_for_sub=False)
        V = _solve(grid)
        vm_ls = np.abs(V)[grid._orig_to_ls]
        bat_bus = list(buses.index).index(net.get_batteries().at["BAT1", "bus_id"])
        self.assertAlmostEqual(vm_olf[bat_bus], VM, places=5)
        self.assertAlmostEqual(vm_ls[bat_bus], VM, places=8)
        np.testing.assert_allclose(vm_ls, vm_olf, atol=2e-4, rtol=0)


if __name__ == "__main__":
    unittest.main()
