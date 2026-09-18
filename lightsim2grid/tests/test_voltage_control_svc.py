# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Static Var Compensators (SVC) through the bordered VoltageControl NRSystem
extension (configured directly on a case14 model, no pypowsybl needed)."""

import warnings
import unittest
import numpy as np
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower

# SvcContainer.RegulationMode
OFF, VOLTAGE, REACTIVE_POWER = 0, 1, 2


def _model():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        net = pn.case14()
        model = init_from_pandapower(net)
    return net, model


def _add_svc(model, *, mode, bus, reg_bus=None, target_vm=1.03, target_q=0.0,
             slope_pu=0.0, b_min=-100.0, b_max=100.0):
    if reg_bus is None:
        reg_bus = bus
    model.init_svcs([int(mode)], np.array([target_vm]), np.array([target_q]),
                    np.array([slope_pu]), np.array([b_min]), np.array([b_max]),
                    np.array([reg_bus], dtype=np.int32), np.array([bus], dtype=np.int32))
    model.tell_solver_need_reset()


class TestSvcVoltageControl(unittest.TestCase):
    def setUp(self):
        self.sn_mva = 100.0
        self.tol = 1e-6
        self.local_bus = 9    # a PQ load bus
        self.remote_bus = 10  # another PQ load bus

    def _solve(self, model, nb_bus, tol=1e-11):
        V = model.ac_pf(np.ones(nb_bus, dtype=np.complex128), 30, tol)
        self.assertGreater(V.shape[0], 0, "ac_pf diverged with an SVC")
        return V

    def test_local_voltage_no_slope(self):
        net, model = _model()
        _add_svc(model, mode=VOLTAGE, bus=self.local_bus, target_vm=1.03)
        V = self._solve(model, net.bus.shape[0])
        self.assertLess(abs(abs(V[self.local_bus]) - 1.03), self.tol)
        svc = model.get_svcs()[0]
        self.assertAlmostEqual(svc.res_p_mw, 0.0, places=10)  # SVC has no active power

    def test_remote_voltage_no_slope(self):
        net, model = _model()
        _add_svc(model, mode=VOLTAGE, bus=self.local_bus, reg_bus=self.remote_bus, target_vm=1.04)
        V = self._solve(model, net.bus.shape[0])
        self.assertLess(abs(abs(V[self.remote_bus]) - 1.04), self.tol)

    def test_local_voltage_with_slope(self):
        # the sloped voltage constraint is  Vm(reg) = vset - s_pu * Q_pu  (Q = injection)
        net, model = _model()
        slope = 0.02
        _add_svc(model, mode=VOLTAGE, bus=self.local_bus, target_vm=1.03, slope_pu=slope)
        V = self._solve(model, net.bus.shape[0])
        svc = model.get_svcs()[0]
        q_pu = svc.res_q_mvar / self.sn_mva
        expected = 1.03 - slope * q_pu
        self.assertLess(abs(abs(V[self.local_bus]) - expected), self.tol,
                        f"slope law broken: |V|={abs(V[self.local_bus]):.8f} expected {expected:.8f}")
        # with a slope the regulated magnitude is NOT exactly the setpoint
        self.assertGreater(abs(abs(V[self.local_bus]) - 1.03), 1e-9)

    def test_reactive_power_mode_injects_fixed_q(self):
        net, model = _model()
        _add_svc(model, mode=REACTIVE_POWER, bus=self.local_bus, target_q=25.0)
        self._solve(model, net.bus.shape[0])
        svc = model.get_svcs()[0]
        self.assertAlmostEqual(svc.res_q_mvar, 25.0, places=6)
        self.assertAlmostEqual(svc.res_p_mw, 0.0, places=10)

    def test_off_mode_is_noop(self):
        # an OFF svc must not change the solution at all
        net, plain = _model()
        _, withsvc = _model()
        _add_svc(withsvc, mode=OFF, bus=self.local_bus, target_vm=1.03)
        nb = net.bus.shape[0]
        Vp = self._solve(plain, nb)
        Vs = self._solve(withsvc, nb)
        self.assertEqual(np.abs(Vp - Vs).max(), 0.0)

    def test_no_svc_backward_compatible(self):
        # declaring an empty SVC container must be bit-identical to never touching it
        net, plain = _model()
        _, empty = _model()
        empty.init_svcs([], np.array([]), np.array([]), np.array([]),
                        np.array([]), np.array([]),
                        np.array([], dtype=np.int32), np.array([], dtype=np.int32))
        empty.tell_solver_need_reset()
        nb = net.bus.shape[0]
        Vp = self._solve(plain, nb)
        Ve = self._solve(empty, nb)
        self.assertEqual(np.abs(Vp - Ve).max(), 0.0)
        self.assertEqual(plain.get_solver().get_J().nnz, empty.get_solver().get_J().nnz)

    def test_kirchhoff_balance(self):
        net, model = _model()
        _add_svc(model, mode=VOLTAGE, bus=self.local_bus, reg_bus=self.remote_bus,
                 target_vm=1.03, slope_pu=0.01)
        V = self._solve(model, net.bus.shape[0])
        mismatch = model.check_solution(V, False)
        self.assertLess(np.abs(np.real(mismatch)).max(), 1e-6)
        self.assertLess(np.abs(np.imag(mismatch)).max(), 1e-6)

    def test_svc_sharing_bus_with_gen_raises(self):
        # v1: an SVC may only be the sole controller of its regulated bus
        net, model = _model()
        # gen 3 (bus 7) regulates bus 9; an SVC at bus 9 also regulates bus 9
        model.set_gen_regulated_bus(3, self.local_bus)
        _add_svc(model, mode=VOLTAGE, bus=self.local_bus, target_vm=1.03)
        with self.assertRaises(Exception):
            model.ac_pf(np.ones(net.bus.shape[0], dtype=np.complex128), 30, 1e-10)


if __name__ == "__main__":
    unittest.main()
