# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import copy
import json
import os
import unittest

import numpy as np

from lightsim2grid.network import init_from_pf_delta
from lightsim2grid.network.from_pf_delta._pf_delta_to_mpc import network_to_mpc
from lightsim2grid.network.from_matpower._aux_add_branch import get_branch_split

_FIXTURE_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "pf_delta_case14.json")


def _branch_split_keys(network):
    """Reproduces, from the row's `network` dict alone, which branch keys end up as
    lines vs. transformers and in what order -- by reusing `network_to_mpc` (the exact
    function `init_from_pf_delta` uses) and `from_matpower`'s own classification, so
    this never drifts out of sync with the loader itself."""
    branch_keys = sorted(network["branch"], key=int)
    is_trafo = get_branch_split(network_to_mpc(network)["branch"])
    line_keys = [k for k, t in zip(branch_keys, is_trafo) if not t]
    trafo_keys = [k for k, t in zip(branch_keys, is_trafo) if t]
    return line_keys, trafo_keys


def validate_against_row(model, row, tol=1e-5):
    """Compares `model`'s CURRENT (already solved) AC state against the PFΔ row's own
    solved `solution`. Does not run any powerflow itself. Returns a dict with the max
    absolute error for each of "vm", "va", "pf", "qf", "pt", "qt" -- restricted to the
    entries actually present in the solution (inactive elements are absent, not zero)."""
    network = row["network"]
    solution = row["solution"]["solution"]

    bus_keys = sorted(network["bus"], key=int)
    line_keys, trafo_keys = _branch_split_keys(network)

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
        line_keys, trafo_keys = _branch_split_keys(network)
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
        _, trafo_keys = _branch_split_keys(network)
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

    def test_missing_network_key_raises(self):
        with self.assertRaises(RuntimeError):
            init_from_pf_delta(self.row["network"])

    def test_nonzero_g_fr_raises(self):
        row = copy.deepcopy(self.row)
        first_branch_key = next(iter(row["network"]["branch"]))
        row["network"]["branch"][first_branch_key]["g_fr"] = 0.01
        with self.assertRaises(RuntimeError):
            init_from_pf_delta(row)

    def test_shift_is_converted_from_radians_to_degrees(self):
        # case14 has no phase shifters (shift == 0 everywhere), so this exercises the
        # rad -> deg conversion directly on a small synthetic network instead
        network = {
            "baseMVA": 100.0,
            "bus": {"1": {"bus_i": 1, "bus_type": 3, "vmax": 1.1, "vmin": 0.9},
                    "2": {"bus_i": 2, "bus_type": 1, "vmax": 1.1, "vmin": 0.9}},
            "gen": {"1": {"gen_bus": 1, "pg": 0.0, "qg": 0.0, "qmax": 10.0, "qmin": -10.0,
                          "vg": 1.0, "pmax": 100.0, "pmin": 0.0, "gen_status": 1}},
            "load": {}, "shunt": {},
            "branch": {"1": {"f_bus": 1, "t_bus": 2, "br_r": 0.01, "br_x": 0.1,
                             "tap": 1.05, "shift": np.pi / 6, "br_status": 1}},
        }
        mpc = network_to_mpc(network)
        np.testing.assert_allclose(mpc["branch"][0, 9], 30.0)  # SHIFT column, degrees


class TestPFDeltaDCLine(unittest.TestCase):
    """PFΔ's own pglib-derived cases (case14/30/57/118/500/2000) never contain a
    `"dcline"` entry, but PowerModels' schema supports one and `from_matpower` now
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
            # consumed by from_matpower's converter (see _aux_add_dc_line.py)
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


class TestLSGridPFDelta(BaseTests, unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
