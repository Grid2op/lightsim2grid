# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Storage-unit (battery) parity against PowSyBl OpenLoadFlow through the pypowsybl
converter.

IIDM batteries use the *generator* convention for their `target_p` / `target_q`
setpoints (positive = produced), while lightsim2grid stores storage as PQ in the
*load* convention (positive = drawn from the grid). The converter negates the
setpoints; this test pins that sign down with a *non-zero* battery (the
`case_14_storage_iidm` fixture only has `targetP=0`, so it cannot catch a sign
regression). Skipped when pypowsybl is unavailable.
"""

import unittest
import numpy as np

try:
    import pypowsybl as pp
    from lightsim2grid.network import init_from_pypowsybl
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False


def _build_net(target_p=50.0, target_q=10.0, connected=True):
    """2-bus network: slack gen + line + load, and a battery on bus 2."""
    n = pp.network.create_empty()
    n.create_substations(id=["S1", "S2"])
    n.create_voltage_levels(id=["VL1", "VL2"], substation_id=["S1", "S2"],
                            topology_kind=["BUS_BREAKER"] * 2, nominal_v=[400.0, 400.0])
    n.create_buses(id=["B1", "B2"], voltage_level_id=["VL1", "VL2"])
    n.create_lines(id="L", voltage_level1_id="VL1", voltage_level2_id="VL2",
                   bus1_id="B1", bus2_id="B2", r=1.0, x=20.0, g1=0.0, b1=0.0, g2=0.0, b2=0.0)
    n.create_generators(id="G", voltage_level_id="VL1", bus_id="B1",
                        target_p=100.0, target_q=0.0, target_v=400.0,
                        voltage_regulator_on=True, max_p=1000.0, min_p=0.0)
    n.create_loads(id="LD", voltage_level_id="VL2", bus_id="B2", p0=80.0, q0=20.0)
    kwargs = dict(id="BATT", voltage_level_id="VL2", bus_id="B2",
                  max_p=1000.0, min_p=-1000.0, target_p=target_p, target_q=target_q)
    n.create_batteries(**kwargs)
    if not connected:
        n.update_batteries(id="BATT", connected=False)
    return n


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestStoragePypowsybl(unittest.TestCase):
    def setUp(self):
        self.params = pp.loadflow.Parameters(distributed_slack=False, use_reactive_limits=False)
        self.tol_v = 5e-6
        self.tol_a = 5e-6
        self.tol_pq = 1e-4  # MW / MVar

    def _run_ls(self, n):
        n.per_unit = False
        model = init_from_pypowsybl(n, gen_slack_id="G", sort_index=True)
        nb_bus = len(model.get_bus_status())
        V = model.ac_pf(np.ones(nb_bus, dtype=np.complex128), 30, 1e-10)
        self.assertGreater(V.shape[0], 0, "lightsim2grid diverged")
        return model, V

    def test_one_storage_loaded(self):
        model, _ = self._run_ls(_build_net(target_p=50.0, target_q=10.0))
        self.assertEqual(len(model.get_storages()), 1)
        sto = model.get_storages()[0]
        # IIDM generator convention target_p=+50 -> load convention -50 (discharging)
        self.assertAlmostEqual(sto.target_p_mw, -50.0, places=4)
        self.assertAlmostEqual(sto.target_q_mvar, -10.0, places=4)

    def test_parity_voltage(self):
        n = _build_net(target_p=50.0, target_q=10.0)
        ref = pp.loadflow.run_ac(n, parameters=self.params)
        self.assertEqual(ref[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        n.per_unit = False
        buses_si = n.get_buses()
        vlevels = n.get_voltage_levels()
        model, V = self._run_ls(n)
        for i, (bus_id, row) in enumerate(buses_si.iterrows()):
            nominal_v = vlevels.loc[row["voltage_level_id"], "nominal_v"]
            v_olf = row["v_mag"] / nominal_v
            a_olf = np.deg2rad(row["v_angle"])
            self.assertLess(abs(abs(V[i]) - v_olf), self.tol_v, f"V at {bus_id}")
            self.assertLess(abs(np.angle(V[i]) - a_olf), self.tol_a, f"angle at {bus_id}")

    def test_parity_storage_pq_sign(self):
        """lightsim2grid storage result p/q (load convention) must match the
        pypowsybl battery output p/q columns (also load convention)."""
        n = _build_net(target_p=50.0, target_q=10.0)
        ref = pp.loadflow.run_ac(n, parameters=self.params)
        self.assertEqual(ref[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        n.per_unit = False
        batt = n.get_batteries()
        p_olf = batt.loc["BATT", "p"]
        q_olf = batt.loc["BATT", "q"]
        model, _ = self._run_ls(n)
        res_p, res_q, _ = model.get_storages_res()
        self.assertLess(abs(res_p[0] - p_olf), self.tol_pq, f"storage p: ls={res_p[0]} olf={p_olf}")
        self.assertLess(abs(res_q[0] - q_olf), self.tol_pq, f"storage q: ls={res_q[0]} olf={q_olf}")

    def test_disconnected_storage(self):
        n = _build_net(target_p=50.0, target_q=10.0, connected=False)
        model, V = self._run_ls(n)
        self.assertEqual(len(model.get_storages()), 1)
        self.assertFalse(model.get_storages_status()[0])


if __name__ == "__main__":
    unittest.main()
