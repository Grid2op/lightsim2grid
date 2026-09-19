# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Storage units taking part in the distributed slack (``add_storage_slackbus``), as
OpenLoadFlow distributes it on batteries.

The oracle for the C++ side is a generator: a storage unit carrying a slack weight must
land on exactly the solution of a generator carrying the same weight on the same bus,
its share of the slack written back to its active result in the load convention. The
oracle for the pypowsybl side is OpenLoadFlow itself, with its distributed slack outer
loop on (``get_pypowsybl_loopfree_distributed_slack_parameters``)."""

import copy
import os
import pickle
import re
import tempfile
import unittest
import warnings

import numpy as np

from lightsim2grid.network import LSGrid

try:
    import pypowsybl as pypo
    import pypowsybl.loadflow as pypo_lf
    from lightsim2grid.network import init_from_pypowsybl
    from lightsim2grid.network.from_pypowsybl import (
        bake_outer_loops,
        get_pypowsybl_loopfree_distributed_slack_parameters,
    )
    PYPO_AVAILABLE = True
except ImportError:
    PYPO_AVAILABLE = False


TOL = 1e-10
BAT_P = 10.         # IIDM generator convention: the battery discharges 10 MW
EXTRA_LOAD = 50.    # MW added to a load, so that the slack has something to distribute


def _network(bat_p=BAT_P, min_p=-300., max_p=300., bus="B3", regulating=False, target_v_pu=1.03):
    """IEEE 14 with a battery ``BAT`` on ``bus`` (B3 also holds generator B3-G, dispatched
    at 0 MW; B4 holds no generator)."""
    net = pypo.network.create_ieee14()
    vl = "VL" + bus[1:]
    net.create_batteries(id="BAT", voltage_level_id=vl, bus_id=bus, target_p=bat_p, target_q=0.,
                         min_p=min_p, max_p=max_p)
    if regulating:
        net.create_minmax_reactive_limits(id="BAT", min_q=-100., max_q=100.)
        nominal_v = net.get_voltage_levels().at[vl, "nominal_v"]
        net.create_extensions("voltageRegulation", id="BAT", voltage_regulator_on=True,
                              target_v=target_v_pu * nominal_v)
    net.update_loads(id="B3-L", p0=net.get_loads().at["B3-L", "p0"] + EXTRA_LOAD)
    return net


def _grid(net, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return init_from_pypowsybl(net, sort_index=False, buses_for_sub=False, **kwargs)


def _solve(grid, ac=True):
    n_bus = grid.get_bus_vn_kv().shape[0]
    V = grid.dc_pf(np.ones(n_bus, dtype=complex), 10, 1e-8)
    assert V.shape[0] > 0, "the DC powerflow did not converge"
    if ac:
        V = grid.ac_pf(V.copy(), 30, TOL)
        assert V.shape[0] > 0, "the AC powerflow did not converge"
    return V


def _gen_id(grid, name):
    return [g.name for g in grid.get_generators()].index(name)


def _olf_params():
    params = get_pypowsybl_loopfree_distributed_slack_parameters()
    # distribute (nearly) all of the mismatch: OLF stops at 1 MW by default
    provider = dict(params.provider_parameters or {})
    provider.update({"slackBusPMaxMismatch": "1e-6", "newtonRaphsonConvEpsPerEq": "1e-11"})
    params.provider_parameters = provider
    return params


@unittest.skipIf(not PYPO_AVAILABLE, "pypowsybl is not installed")
class TestStorageSlack(unittest.TestCase):
    """The C++ side: weights, the result split, the invalidation flags, the state."""

    def setUp(self):
        self.net = _network()

    def test_share_equals_a_generator_share_on_the_same_bus(self):
        for ac in (True, False):
            with self.subTest(ac=ac):
                # gen weights normalised to (0.5, 0.5), the storage adds 0.25 raw: 0.4 / 0.4 / 0.2
                grid_sto = _grid(self.net, gen_slack_id={"B1-G": 1., "B2-G": 1.})
                grid_sto.add_storage_slackbus(0, 0.25)
                grid_gen = _grid(self.net, gen_slack_id={"B1-G": 1., "B2-G": 1., "B3-G": 0.5})
                sto = grid_sto.get_storages()[0]
                self.assertTrue(sto.is_slack)
                self.assertAlmostEqual(sto.slack_weight, 0.25)
                V_sto = _solve(grid_sto, ac)
                V_gen = _solve(grid_gen, ac)
                np.testing.assert_allclose(V_sto, V_gen, atol=1e-8, rtol=0)
                np.testing.assert_allclose(grid_sto.get_slack_weights_solver(),
                                           grid_gen.get_slack_weights_solver(), atol=1e-12, rtol=0)

                g3 = _gen_id(grid_gen, "B3-G")
                share_gen = grid_gen.get_gen_res()[0][g3] - grid_gen.get_generators()[g3].target_p_mw
                # load convention: the unit draws -BAT_P, minus the share it produces
                share_sto = -(grid_sto.get_storages_res()[0][0] - (-BAT_P))
                self.assertGreater(abs(share_gen), 1., "the test needs a real share of the slack")
                self.assertAlmostEqual(share_sto, share_gen, delta=1e-6)
                for name in ("B1-G", "B2-G"):
                    self.assertAlmostEqual(grid_sto.get_gen_res()[0][_gen_id(grid_sto, name)],
                                           grid_gen.get_gen_res()[0][_gen_id(grid_gen, name)], delta=1e-6)

    def test_split_between_a_generator_and_a_storage_unit_on_one_bus(self):
        for ac in (True, False):
            with self.subTest(ac=ac):
                # raw (1, 1, 0.25) / 2.25 for the generators and 0.25 / 2.25 for the storage:
                # bus B3 carries 0.5 / 2.25, the same bus weights as (1, 1, 0.5) / 2.5
                grid_mix = _grid(self.net, gen_slack_id={"B1-G": 1., "B2-G": 1., "B3-G": 0.25})
                grid_mix.add_storage_slackbus(0, 0.25 / 2.25)
                grid_gen = _grid(self.net, gen_slack_id={"B1-G": 1., "B2-G": 1., "B3-G": 0.5})
                np.testing.assert_allclose(_solve(grid_mix, ac), _solve(grid_gen, ac), atol=1e-8, rtol=0)

                g3 = _gen_id(grid_gen, "B3-G")
                bus_share = grid_gen.get_gen_res()[0][g3]
                share_gen = grid_mix.get_gen_res()[0][g3]
                share_sto = -(grid_mix.get_storages_res()[0][0] + BAT_P)
                self.assertAlmostEqual(share_gen, share_sto, delta=1e-6)
                self.assertAlmostEqual(share_gen + share_sto, bus_share, delta=1e-6)

    def test_disconnecting_a_storage_participant_reweights_the_slack(self):
        grid = _grid(self.net, gen_slack_id={"B1-G": 1., "B2-G": 1.})
        grid.add_storage_slackbus(0, 0.25)
        V_on = _solve(grid)
        gen_p_on = grid.get_gen_res()[0].copy()

        grid.deactivate_storage(0)
        V_off = _solve(grid)
        ref = _grid(self.net, gen_slack_id={"B1-G": 1., "B2-G": 1.})
        ref.deactivate_storage(0)
        V_ref = _solve(ref)
        np.testing.assert_allclose(V_off, V_ref, atol=1e-8, rtol=0)
        np.testing.assert_allclose(grid.get_gen_res()[0], ref.get_gen_res()[0], atol=1e-6, rtol=0)

        grid.reactivate_storage(0)
        np.testing.assert_allclose(_solve(grid), V_on, atol=1e-8, rtol=0)
        np.testing.assert_allclose(grid.get_gen_res()[0], gen_p_on, atol=1e-6, rtol=0)

    def test_storage_as_the_only_participant(self):
        grid = _grid(self.net, gen_slack_id={"B1-G": 1.})
        grid.remove_gen_slackbus(_gen_id(grid, "B1-G"))
        grid.add_storage_slackbus(0, 1.)
        grid.check_grid()
        _solve(grid)
        gens = grid.get_generators()
        # the storage unit absorbs everything: every generator stays at its setpoint
        np.testing.assert_allclose(grid.get_gen_res()[0], [g.target_p_mw for g in gens], atol=1e-6, rtol=0)
        self.assertGreater(abs(grid.get_storages_res()[0][0] + BAT_P), 1.)

        # and a grid whose only participant is disconnected is refused
        grid.deactivate_storage(0)
        with self.assertRaises(RuntimeError):
            grid.check_grid()

    def test_remove_storage_slackbus(self):
        grid = _grid(self.net, gen_slack_id={"B1-G": 1., "B2-G": 1.})
        V_ref = _solve(grid)
        grid.add_storage_slackbus(0, 0.25)
        _solve(grid)
        grid.remove_storage_slackbus(0)
        self.assertFalse(grid.get_storages()[0].is_slack)
        np.testing.assert_allclose(_solve(grid), V_ref, atol=1e-8, rtol=0)
        self.assertAlmostEqual(grid.get_storages_res()[0][0], -BAT_P, delta=1e-8)

        with self.assertRaises(RuntimeError):
            grid.add_storage_slackbus(1, 1.)
        with self.assertRaises(RuntimeError):
            grid.add_storage_slackbus(0, 0.)

    def test_state_round_trips(self):
        grid = _grid(self.net, gen_slack_id={"B1-G": 1., "B2-G": 1.})
        grid.add_storage_slackbus(0, 0.25)
        V = _solve(grid)

        grid_pickle = pickle.loads(pickle.dumps(grid))
        self.assertTrue(grid_pickle.get_storages()[0].is_slack)
        self.assertAlmostEqual(grid_pickle.get_storages()[0].slack_weight, 0.25)
        np.testing.assert_allclose(_solve(grid_pickle), V, atol=1e-8, rtol=0)

        with tempfile.TemporaryDirectory() as tmp_dir:
            path = os.path.join(tmp_dir, "grid.lsb")
            grid.save_binary(path)
            grid_binary = LSGrid.load_binary(path)
        self.assertTrue(grid_binary.get_storages()[0].is_slack)
        self.assertAlmostEqual(grid_binary.get_storages()[0].slack_weight, 0.25)
        np.testing.assert_allclose(_solve(grid_binary), V, atol=1e-8, rtol=0)

    def test_scenario_sweep_row_with_a_participating_generator_off(self):
        """A sweep row taking a slack generator out re-weights the slack over what is left,
        the storage participant included (LSGrid::get_slack_weights_solver_without)."""
        from lightsim2grid.scenarioSweep import ScenarioSweepCPP
        grid = _grid(self.net, gen_slack_id={"B1-G": 1., "B2-G": 1.})
        grid.add_storage_slackbus(0, 0.25)
        V0 = _solve(grid)
        b2 = _gen_id(grid, "B2-G")
        gen_mask = np.zeros((2, len(grid.get_generators())), dtype=bool)
        gen_mask[1, b2] = True
        sweep = ScenarioSweepCPP(grid)
        sweep.set_contingency_gens(gen_mask)
        sweep.compute(1.0 * V0, 30, TOL)
        self.assertTrue(np.all(sweep.converged_mask()))

        ref = copy.deepcopy(grid)
        ref.deactivate_gen(b2)
        V_ref = ref.ac_pf(1.0 * V0, 30, TOL)
        self.assertGreater(V_ref.shape[0], 0)
        buses = np.asarray(grid.id_ac_solver_to_me(), dtype=int)
        np.testing.assert_allclose(sweep.get_voltages()[0][buses], V0[buses], atol=1e-8, rtol=1e-8)
        np.testing.assert_allclose(sweep.get_voltages()[1][buses], V_ref[buses], atol=1e-8, rtol=1e-8)

    def test_slack_bus_held_by_a_regulating_storage_keeps_its_voltage(self):
        """A slack bus whose magnitude a battery holds is Vm-fixed (VoltageControlPlan::
        build_free_vm_slack): it used to get a free Vm, the setpoint being ignored."""
        net = _network(bus="B4", regulating=True, target_v_pu=1.03)
        grid = _grid(net, gen_slack_id={"B1-G": 1., "B2-G": 1.})
        grid.add_storage_slackbus(0, 0.25)
        V = _solve(grid)
        sto = grid.get_storages()[0]
        self.assertTrue(sto.voltage_regulator_on)
        self.assertAlmostEqual(abs(V[sto.bus_id]), 1.03, places=8)


@unittest.skipIf(not PYPO_AVAILABLE, "pypowsybl is not installed")
class TestStorageSlackOpenLoadFlow(unittest.TestCase):
    """The pypowsybl side: the default distributed slack includes the batteries OpenLoadFlow
    distributes on, with its weights."""

    def _assert_matches_olf(self, net, grid):
        res = pypo_lf.run_ac(net, _olf_params())
        self.assertEqual(res[0].status, pypo_lf.ComponentStatus.CONVERGED)
        V = _solve(grid)
        # results: IIDM p is in the load convention, lightsim2grid's generators are not
        self.assertAlmostEqual(grid.get_storages_res()[0][0], net.get_batteries().at["BAT", "p"], delta=1e-3)
        np.testing.assert_allclose(grid.get_gen_res()[0], -net.get_generators()["p"].to_numpy(), atol=1e-3, rtol=0)
        buses = net.get_buses()
        nominal_v = net.get_voltage_levels()["nominal_v"].reindex(buses["voltage_level_id"]).to_numpy()
        np.testing.assert_allclose(np.abs(V)[grid._orig_to_ls], buses["v_mag"].to_numpy() / nominal_v,
                                   atol=1e-5, rtol=0)
        return net.get_batteries().at["BAT", "p"]

    def test_matches_openloadflow(self):
        cases = {
            "discharging": (10., -300., 300., True),
            "charging": (-10., -300., 300., True),   # a charging battery participates too
            "idle": (0., -300., 300., False),
            "above max_p": (10., -5., 5., False),
            "below min_p": (-10., -5., 300., False),
        }
        for label, (bat_p, min_p, max_p, participates) in cases.items():
            with self.subTest(label):
                net = _network(bat_p=bat_p, min_p=min_p, max_p=max_p, bus="B4")
                grid = _grid(net)
                sto = grid.get_storages()[0]
                self.assertEqual(sto.is_slack, participates)
                olf_p = self._assert_matches_olf(net, grid)
                share = -olf_p - bat_p
                if participates:
                    self.assertGreater(abs(share), 0.1)
                    # OLF's key, max_p / droop (default droop 4), in the unit of the generators'
                    b1 = grid.get_generators()[_gen_id(grid, "B1-G")]
                    self.assertAlmostEqual(sto.slack_weight / b1.slack_weight, (max_p / 4.) / (9999. / 4.), places=9)
                else:
                    self.assertAlmostEqual(share, 0., delta=1e-6)

    @staticmethod
    def _with_battery_apc(net, **attrs):
        """``net`` round-tripped through XIIDM with an ``activePowerControl`` extension on the
        battery (pypowsybl <= 1.16.1 cannot create it on a battery)."""
        net.create_extensions("activePowerControl", id="B2-G", participate=True, droop=4.)  # OLF's default
        xml = net.save_to_string("XIIDM")
        block = re.search(r'<iidm:extension id="B2-G">.*?</iidm:extension>', xml, re.S).group(0)
        attr_txt = " ".join(f'{key}="{val}"' for key, val in attrs.items())
        bat_block = re.sub(r"<(\w+):activePowerControl[^>]*/>", rf"<\1:activePowerControl {attr_txt}/>", block)
        bat_block = bat_block.replace('id="B2-G"', 'id="BAT"')
        return pypo.network.load_from_string("net.xiidm", xml.replace("</iidm:network>", bat_block + "</iidm:network>"))

    def test_battery_extension(self):
        cases = {
            "not participating": ({"participate": "false", "droop": "4.0"}, None),
            "droop 2": ({"participate": "true", "droop": "2.0"}, 2.),
            "droop 0": ({"participate": "true", "droop": "0.0"}, None),
        }
        for label, (attrs, droop) in cases.items():
            with self.subTest(label):
                net = self._with_battery_apc(_network(bus="B4"), **attrs)
                grid = _grid(net)
                sto = grid.get_storages()[0]
                self.assertEqual(sto.is_slack, droop is not None)
                if droop is not None:
                    b1 = grid.get_generators()[_gen_id(grid, "B1-G")]
                    self.assertAlmostEqual(sto.slack_weight / b1.slack_weight, (300. / droop) / (9999. / 4.), places=9)
                olf_p = self._assert_matches_olf(net, grid)
                if droop is None:
                    self.assertAlmostEqual(olf_p, -BAT_P, delta=1e-6)

                # without reading the extension, the battery gets OLF's defaults
                grid_default = _grid(net, battery_active_power_control="default")
                b1 = grid_default.get_generators()[_gen_id(grid_default, "B1-G")]
                self.assertTrue(grid_default.get_storages()[0].is_slack)
                self.assertAlmostEqual(grid_default.get_storages()[0].slack_weight / b1.slack_weight,
                                       300. / 9999., places=9)

    def test_explicit_slack_leaves_the_batteries_out(self):
        grid = _grid(_network(bus="B4"), gen_slack_id={"B1-G": 1., "B2-G": 1.})
        self.assertFalse(grid.get_storages()[0].is_slack)
        with self.assertRaises(RuntimeError):
            _grid(_network(bus="B4"), battery_active_power_control="nope")

    def test_bake_excludes_a_capped_charging_battery(self):
        """A charging battery is capped at 0 MW by a positive mismatch (OLF never pushes a unit
        across 0): the bake takes it out of the slack, and lightsim2grid then reproduces OLF."""
        net = _network(bat_p=-10., min_p=-9000., max_p=9000., bus="B4")
        res = pypo_lf.run_ac(net, _olf_params())
        self.assertEqual(res[0].status, pypo_lf.ComponentStatus.CONVERGED)
        self.assertAlmostEqual(net.get_batteries().at["BAT", "p"], 0., delta=1e-6)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bake_outer_loops(net)
        grid = _grid(net)
        self.assertFalse(grid.get_storages()[0].is_slack)
        _solve(grid)
        np.testing.assert_allclose(grid.get_gen_res()[0], net.get_generators()["target_p"].to_numpy(),
                                   atol=1e-3, rtol=0)
        self.assertAlmostEqual(grid.get_storages_res()[0][0], 0., delta=1e-6)


if __name__ == "__main__":
    unittest.main()
