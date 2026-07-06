# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import json
import os
import unittest

import numpy as np

from lightsim2grid.network import init_from_pf_delta, init_from_powermodels
from lightsim2grid.network.from_powermodels._aux_add_branch import classify_branches

_FIXTURE_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "pf_delta_case14.json")


def validate_against_row(model, row, tol=1e-5):
    """Compares `model`'s CURRENT (already solved) AC state against the PFΔ row's own
    solved `solution`. Does not run any powerflow itself. Returns a dict with the max
    absolute error for each of "vm", "va", "pf", "qf", "pt", "qt" -- restricted to the
    entries actually present in the solution (inactive elements are absent, not zero)."""
    network = row["network"]
    solution = row["solution"]["solution"]

    bus_keys = sorted(network["bus"], key=int)
    line_keys, trafo_keys = classify_branches(network)

    vm = model.get_Vm()
    va = model.get_Va()
    errors = {"vm": 0.0, "va": 0.0}
    for ls_id, bus_key in enumerate(bus_keys):
        bus_sol = solution["bus"].get(bus_key)
        if bus_sol is None:
            continue
        errors["vm"] = max(errors["vm"], abs(vm[ls_id] - bus_sol["vm"]))
        errors["va"] = max(errors["va"], abs(va[ls_id] - bus_sol["va"]))

    por, qor, _, _ = model.get_line_res1()
    pex, qex, _, _ = model.get_line_res2()
    phv, qhv, _, _ = model.get_trafo_res1()
    plv, qlv, _, _ = model.get_trafo_res2()

    errors.update({"pf": 0.0, "qf": 0.0, "pt": 0.0, "qt": 0.0})
    for ls_id, branch_key in enumerate(line_keys):
        branch_sol = solution["branch"].get(branch_key)
        if branch_sol is None:
            continue
        errors["pf"] = max(errors["pf"], abs(por[ls_id] - branch_sol["pf"]))
        errors["qf"] = max(errors["qf"], abs(qor[ls_id] - branch_sol["qf"]))
        errors["pt"] = max(errors["pt"], abs(pex[ls_id] - branch_sol["pt"]))
        errors["qt"] = max(errors["qt"], abs(qex[ls_id] - branch_sol["qt"]))
    for ls_id, branch_key in enumerate(trafo_keys):
        branch_sol = solution["branch"].get(branch_key)
        if branch_sol is None:
            continue
        errors["pf"] = max(errors["pf"], abs(phv[ls_id] - branch_sol["pf"]))
        errors["qf"] = max(errors["qf"], abs(qhv[ls_id] - branch_sol["qf"]))
        errors["pt"] = max(errors["pt"], abs(plv[ls_id] - branch_sol["pt"]))
        errors["qt"] = max(errors["qt"], abs(qlv[ls_id] - branch_sol["qt"]))

    return errors


class BaseTests:
    def setUp(self):
        with open(_FIXTURE_PATH, "r") as f:
            self.row = json.load(f)
        self.model = init_from_pf_delta(self.row)

    def test_counts(self):
        network = self.row["network"]
        line_keys, trafo_keys = classify_branches(network)
        assert len(self.model.get_lines()) == len(line_keys)
        assert len(self.model.get_trafos()) == len(trafo_keys)
        assert len(self.model.get_generators()) == len(network["gen"])
        assert len(self.model.get_loads()) == len(network["load"])
        assert len(self.model.get_shunts()) == len(network.get("shunt", {}))

    def test_fixture_exercises_taps_and_shunt(self):
        # sanity check on the fixture itself: if these were all trivial (tap==1,
        # shunt absent) the AC powerflow test below would not actually be exercising
        # the tap-ratio / shunt-sign conversions
        network = self.row["network"]
        _, trafo_keys = classify_branches(network)
        assert len(trafo_keys) > 0, "fixture has no transformers, tap conversion untested"
        assert any(network["branch"][k]["tap"] != 1.0 for k in trafo_keys)
        assert len(network.get("shunt", {})) > 0, "fixture has no shunt, sign convention untested"

    def test_ac_pf_matches_solution(self):
        nb_bus = len(self.row["network"]["bus"])
        V0 = np.full(nb_bus, 1.0, dtype=complex)
        Vfinal = self.model.ac_pf(V0, 10, 1e-8)
        assert Vfinal.shape[0] > 0, "powerflow diverged"

        errors = validate_against_row(self.model, self.row, tol=1e-5)
        for name, err in errors.items():
            assert err <= 1e-5, f"error too large for {name}: {err}"

    def test_accepts_path(self):
        model = init_from_pf_delta(_FIXTURE_PATH)
        nb_bus = len(self.row["network"]["bus"])
        Vfinal = model.ac_pf(np.full(nb_bus, 1.0, dtype=complex), 10, 1e-8)
        assert Vfinal.shape[0] > 0, "powerflow diverged"

    def test_init_from_powermodels_on_bare_network_dict(self):
        # init_from_pf_delta is a thin wrapper unwrapping row["network"] and calling
        # init_from_powermodels directly -- both should behave identically
        model = init_from_powermodels(self.row["network"])
        assert len(model.get_lines()) == len(self.model.get_lines())
        assert len(model.get_trafos()) == len(self.model.get_trafos())
        nb_bus = len(self.row["network"]["bus"])
        Vfinal = model.ac_pf(np.full(nb_bus, 1.0, dtype=complex), 10, 1e-8)
        assert Vfinal.shape[0] > 0, "powerflow diverged"

    def test_missing_network_key_raises(self):
        with self.assertRaises(RuntimeError):
            init_from_pf_delta(self.row["network"])


class TestLSGridPFDelta(BaseTests, unittest.TestCase):
    pass


class TestPFDeltaBranchConductance(unittest.TestCase):
    """Unlike `init_from_matpower` (which routes through matpower's raw `mpc.branch`,
    whose single `BR_B` column cannot represent line conductance at all),
    `init_from_powermodels` builds lines directly from `"g_fr"`/`"g_to"`, so an
    asymmetric, non-zero line conductance should convert and solve without issue."""

    @staticmethod
    def _network(g_fr=0.0, g_to=0.0):
        return {
            "baseMVA": 100.0,
            "bus": {"1": {"bus_i": 1, "bus_type": 3, "vmax": 1.1, "vmin": 0.9},
                    "2": {"bus_i": 2, "bus_type": 1, "vmax": 1.1, "vmin": 0.9}},
            "gen": {"1": {"gen_bus": 1, "pg": 0.0, "qg": 0.0, "qmax": 10.0, "qmin": -10.0,
                          "vg": 1.0, "pmax": 100.0, "pmin": 0.0, "gen_status": 1}},
            "load": {"1": {"load_bus": 2, "pd": 10.0, "qd": 2.0, "status": 1}},
            "shunt": {},
            "branch": {"1": {"f_bus": 1, "t_bus": 2, "br_r": 0.01, "br_x": 0.1,
                             "g_fr": g_fr, "g_to": g_to, "br_status": 1}},
        }

    def test_asymmetric_line_conductance_converts_and_solves(self):
        model = init_from_powermodels(self._network(g_fr=0.001, g_to=0.002))
        assert len(model.get_lines()) == 1
        Vfinal = model.ac_pf(np.full(2, 1.0, dtype=complex), 10, 1e-8)
        assert Vfinal.shape[0] > 0, "powerflow diverged"


class TestPFDeltaShiftConversion(unittest.TestCase):
    """case14 has no phase shifters (`shift == 0` everywhere in the real pglib case),
    so this exercises the radians -> degrees -> radians round trip through the
    `init_trafo` API boundary directly on a small synthetic network instead."""

    def test_shift_round_trips_through_init_trafo(self):
        shift_rad = np.pi / 6
        network = {
            "baseMVA": 100.0,
            "bus": {"1": {"bus_i": 1, "bus_type": 3, "vmax": 1.1, "vmin": 0.9},
                    "2": {"bus_i": 2, "bus_type": 1, "vmax": 1.1, "vmin": 0.9}},
            "gen": {"1": {"gen_bus": 1, "pg": 0.0, "qg": 0.0, "qmax": 10.0, "qmin": -10.0,
                          "vg": 1.0, "pmax": 100.0, "pmin": 0.0, "gen_status": 1}},
            "load": {}, "shunt": {},
            "branch": {"1": {"f_bus": 1, "t_bus": 2, "br_r": 0.01, "br_x": 0.1,
                             "tap": 1.05, "shift": shift_rad, "br_status": 1}},
        }
        model = init_from_powermodels(network)
        trafos = model.get_trafos()
        assert len(trafos) == 1
        assert abs(trafos[0].shift_rad - shift_rad) <= 1e-8


class TestPFDeltaStorage(unittest.TestCase):
    """PFΔ's own pglib-derived cases never contain a `"storage"` entry, but
    PowerModels' schema supports one and lightsim2grid itself supports storage
    (unlike `from_matpower`'s raw mpc format, which has no storage table at all), so
    this exercises that conversion path directly."""

    @staticmethod
    def _network_with_storage(status=1):
        return {
            "baseMVA": 100.0,
            "bus": {"1": {"bus_i": 1, "bus_type": 3, "vmax": 1.1, "vmin": 0.9},
                    "2": {"bus_i": 2, "bus_type": 1, "vmax": 1.1, "vmin": 0.9}},
            "gen": {"1": {"gen_bus": 1, "pg": 0.0, "qg": 0.0, "qmax": 100.0, "qmin": -100.0,
                          "vg": 1.0, "pmax": 200.0, "pmin": 0.0, "gen_status": 1}},
            "load": {"1": {"load_bus": 2, "pd": 10.0, "qd": 2.0, "status": 1}},
            "shunt": {},
            "branch": {"1": {"f_bus": 1, "t_bus": 2, "br_r": 0.01, "br_x": 0.1, "br_status": 1}},
            "storage": {"1": {"storage_bus": 2, "ps": 5.0, "qs": 1.0, "status": status}},
        }

    def test_no_storage_key_at_all(self):
        network = self._network_with_storage()
        del network["storage"]
        model = init_from_powermodels(network)
        assert len(model.get_storages()) == 0

    def test_storage_is_converted(self):
        model = init_from_powermodels(self._network_with_storage())
        storages = model.get_storages()
        assert len(storages) == 1
        # PowerModels' "ps"/"qs" are load-convention (withdrawn=positive), matching
        # lightsim2grid's own init_storages convention directly -- no sign flip
        Vfinal = model.ac_pf(np.full(2, 1.0, dtype=complex), 10, 1e-8)
        assert Vfinal.shape[0] > 0, "powerflow diverged"

    def test_deactivated_storage(self):
        model = init_from_powermodels(self._network_with_storage(status=0))
        Vfinal = model.ac_pf(np.full(2, 1.0, dtype=complex), 10, 1e-8)
        assert Vfinal.shape[0] > 0, "powerflow diverged"


class TestPFDeltaDCLine(unittest.TestCase):
    """PFΔ's own pglib-derived cases (case14/30/57/118/500/2000) never contain a
    `"dcline"` entry, but PowerModels' schema supports one and lightsim2grid now
    converts it (`model.init_dclines`), so this exercises that translation path
    directly rather than leaving it silently unsupported."""

    @staticmethod
    def _network_with_dcline():
        return {
            "baseMVA": 100.0,
            "bus": {"1": {"bus_i": 1, "bus_type": 3, "vmax": 1.1, "vmin": 0.9},
                    "2": {"bus_i": 2, "bus_type": 1, "vmax": 1.1, "vmin": 0.9},
                    "3": {"bus_i": 3, "bus_type": 1, "vmax": 1.1, "vmin": 0.9}},
            "gen": {"1": {"gen_bus": 1, "pg": 0.0, "qg": 0.0, "qmax": 100.0, "qmin": -100.0,
                          "vg": 1.0, "pmax": 200.0, "pmin": 0.0, "gen_status": 1}},
            "load": {"1": {"load_bus": 2, "pd": 10.0, "qd": 2.0, "status": 1},
                     "2": {"load_bus": 3, "pd": 5.0, "qd": 1.0, "status": 1}},
            "shunt": {},
            "branch": {"1": {"f_bus": 1, "t_bus": 2, "br_r": 0.01, "br_x": 0.1, "br_status": 1},
                       "2": {"f_bus": 2, "t_bus": 3, "br_r": 0.01, "br_x": 0.1, "br_status": 1}},
            # PowerModels keeps "pf" at matpower's own sign, but negates "pt"/"qf"/"qt";
            # only f_bus/t_bus/br_status/pf/loss0/loss1/vf/vt/qmin*/qmax* are actually
            # consumed by the converter (see from_powermodels/_aux_add_dc_line.py)
            "dcline": {"1": {"f_bus": 2, "t_bus": 3, "pf": 10.0, "pt": -9.0, "qf": 0.0, "qt": 0.0,
                             "vf": 1.0, "vt": 1.0, "qminf": -10.0, "qmaxf": 10.0,
                             "qmint": -10.0, "qmaxt": 10.0, "loss0": 0.5, "loss1": 0.05,
                             "br_status": 1}},
        }

    def test_no_dcline_key_at_all(self):
        network = self._network_with_dcline()
        del network["dcline"]
        model = init_from_pf_delta({"network": network})
        assert len(model.get_dclines()) == 0

    def test_dcline_is_converted_and_matches_matpower_loss_formula(self):
        model = init_from_pf_delta({"network": self._network_with_dcline()})
        dclines = model.get_dclines()
        assert len(dclines) == 1
        # matpower's PF (10 MW, "from" -> "to") maps to lightsim2grid's "power received
        # at side 1" convention, so it must be negated (same as the pandapower loader)
        assert abs(dclines[0].target_p1_mw - (-10.0)) <= 1e-8

        Vfinal = model.ac_pf(np.full(3, 1.0, dtype=complex), 10, 1e-8)
        assert Vfinal.shape[0] > 0, "powerflow diverged"
        dclines = model.get_dclines()
        # PT = PF - (LOSS0 + LOSS1 * PF) = 10 - (0.5 + 0.05 * 10) = 9.0
        assert abs(dclines[0].res_p1_mw - (-10.0)) <= 1e-4
        assert abs(dclines[0].res_p2_mw - 9.0) <= 1e-4

    def test_deactivated_dcline(self):
        network = self._network_with_dcline()
        network["dcline"]["1"]["br_status"] = 0
        model = init_from_pf_delta({"network": network})
        Vfinal = model.ac_pf(np.full(3, 1.0, dtype=complex), 10, 1e-8)
        assert Vfinal.shape[0] > 0, "powerflow diverged"


if __name__ == "__main__":
    unittest.main()
