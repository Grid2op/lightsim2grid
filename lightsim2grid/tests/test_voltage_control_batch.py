# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The batch algorithms (``TimeSeriesCPP`` / ``ContingencyAnalysisCPP``) must apply
remote voltage control -- a remote-regulating generator or a voltage-mode SVC --
exactly as ``ac_pf`` does.

They used not to. The ``VoltageControl`` NR extension is not built from the solver's
arguments: it reaches back into the ``LSGrid`` through the algorithm's ``lsgrid_ptr``
and reads that grid's member solver-side bus labelling. ``ac_pf`` passes its own
members to ``pre_process_solver`` so they are always current, but the batch classes
work on a private copy of the grid and pass their own local mapping vectors, and
only the FORWARD map was written back. With the reverse map left empty,
``fill_voltage_control_solver_data`` returned on its ``nb_bus_solver == 0`` guard --
so the batch solved a grid with no regulation at all, converged, and returned
plausible wrong voltages instead of raising.

Companion C++ suite: ``src/tests/test_batch_voltage_control.cpp``.
"""

import warnings
import unittest
import numpy as np
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import TimeSeriesCPP, ContingencyAnalysisCPP  # noqa: E402

MAX_IT = 30
TOL = 1e-11
# SvcContainer.RegulationMode
SVC_VOLTAGE_MODE = 1
# gen 3 sits on bus 7; bus 9 is a PQ load bus no generator regulates locally.
GEN_ID = 3
REG_BUS = 9
SVC_BUS = 4       # a PQ bus, so the SVC's own bus carries a reactive equation
SVC_V_SET = 1.06


def _case14():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        net = pn.case14()
        model = init_from_pandapower(net)
    return net, model


def _make(kind):
    """case14 plus one remote voltage controller; ``kind`` selects which."""
    net, model = _case14()
    if kind == "remote_gen":
        model.set_gen_regulated_bus(GEN_ID, REG_BUS)
    elif kind in ("svc_flat", "svc_sloped"):
        slope_pu = 0.02 if kind == "svc_sloped" else 0.
        model.init_svcs([SVC_VOLTAGE_MODE], np.array([SVC_V_SET]), np.array([0.]),
                        np.array([slope_pu]), np.array([-100.]), np.array([100.]),
                        np.array([REG_BUS], dtype=np.int32),
                        np.array([SVC_BUS], dtype=np.int32))
    else:
        raise ValueError(kind)
    model.tell_solver_need_reset()
    return net, model


KINDS = ["remote_gen", "svc_flat", "svc_sloped"]


def _base_injections(model):
    gen_p = np.array(model.get_gen_target_p(), dtype=np.float64)
    load_p = np.array(model.get_load_target_p(), dtype=np.float64)
    load_q = np.array([ld.target_q_mvar for ld in model.get_loads()], dtype=np.float64)
    return gen_p, load_p, load_q


def _single_shot(model, gen_p, load_p, load_q, nb_bus):
    for i, p in enumerate(gen_p):
        model.change_p_gen(i, float(p))
    for i, p in enumerate(load_p):
        model.change_p_load(i, float(p))
    for i, q in enumerate(load_q):
        model.change_q_load(i, float(q))
    return model.ac_pf(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)


class TestVoltageControlIsActuallyOn(unittest.TestCase):
    """Guard the guard: if these setups did not regulate anything even single-shot,
    every comparison below would pass vacuously."""

    def test_single_shot_regulates(self):
        for kind in KINDS:
            with self.subTest(kind=kind):
                net, model = _make(kind)
                nb_bus = net.bus.shape[0]
                V = model.ac_pf(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
                self.assertGreater(V.shape[0], 0, f"ac_pf diverged ({kind})")

                _, plain = _case14()
                Vplain = plain.ac_pf(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
                self.assertGreater(np.abs(abs(V[REG_BUS]) - abs(Vplain[REG_BUS])), 1e-4,
                                   f"{kind} does not move the regulated bus at all")

                if kind in ("remote_gen", "svc_flat"):
                    vset = (model.get_generators()[GEN_ID].target_vm_pu
                            if kind == "remote_gen" else SVC_V_SET)
                    self.assertLess(abs(abs(V[REG_BUS]) - vset), 1e-6,
                                    f"{kind}: |V[{REG_BUS}]| != setpoint")


class TestVoltageControlTimeSeries(unittest.TestCase):
    def test_timeseries_matches_single_shot(self):
        for kind in KINDS:
            with self.subTest(kind=kind):
                net, model = _make(kind)
                nb_bus = net.bus.shape[0]
                gen_p, load_p, load_q = _base_injections(model)

                # two distinct steps (base and a +8% load / +8% gen scaling)
                scales = np.array([1.0, 1.08])
                gen_mat = np.outer(scales, gen_p)
                lp_mat = np.outer(scales, load_p)
                lq_mat = np.outer(scales, load_q)
                sgen_mat = np.zeros((len(scales), 0))

                ts = TimeSeriesCPP(model)
                status = ts.compute_Vs(gen_mat, sgen_mat, lp_mat, lq_mat,
                                       np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
                self.assertEqual(status, 1, f"time series diverged ({kind})")
                Vs = np.array(ts.get_voltages())

                for step in range(len(scales)):
                    _, ref_model = _make(kind)
                    Vref = _single_shot(ref_model, gen_mat[step], lp_mat[step],
                                        lq_mat[step], nb_bus)
                    self.assertGreater(Vref.shape[0], 0, f"{kind} ref step {step} diverged")
                    self.assertLess(np.abs(Vs[step] - Vref).max(), 1e-7,
                                    f"{kind} step {step}")


class TestVoltageControlContingencyAnalysis(unittest.TestCase):
    def test_contingency_matches_single_shot(self):
        # powerlines away from the controller / regulated buses, so every
        # contingency leaves the control problem itself well posed
        cont_lines = [0, 5, 10]
        for kind in KINDS:
            with self.subTest(kind=kind):
                net, model = _make(kind)
                nb_bus = net.bus.shape[0]

                ca = ContingencyAnalysisCPP(model)
                for line_id in cont_lines:
                    ca.add_n1(line_id)
                ca.compute(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
                Vs = np.array(ca.get_voltages())

                for k, line_id in enumerate(cont_lines):
                    if np.abs(Vs[k]).max() == 0.0:
                        # this contingency diverged / islanded: nothing to compare
                        continue
                    _, ref_model = _make(kind)
                    ref_model.deactivate_powerline(line_id)
                    Vref = ref_model.ac_pf(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
                    self.assertGreater(Vref.shape[0], 0,
                                       f"{kind} ref cont {line_id} diverged")
                    mask = np.abs(Vs[k]) > 0.0
                    self.assertLess(np.abs(Vs[k][mask] - Vref[mask]).max(), 1e-6,
                                    f"{kind} contingency line {line_id}")

    def test_regulated_bus_is_held_in_every_contingency(self):
        # the substantive claim, independent of the single-shot reference: with a
        # flat (slope-free) controller the regulated bus sits on the setpoint in
        # every contingency the analysis actually solves.
        for kind, vset in (("remote_gen", None), ("svc_flat", SVC_V_SET)):
            with self.subTest(kind=kind):
                net, model = _make(kind)
                nb_bus = net.bus.shape[0]
                if vset is None:
                    vset = model.get_generators()[GEN_ID].target_vm_pu

                ca = ContingencyAnalysisCPP(model)
                for line_id in [0, 5, 10]:
                    ca.add_n1(line_id)
                ca.compute(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
                Vs = np.array(ca.get_voltages())

                nb_checked = 0
                for k in range(Vs.shape[0]):
                    if np.abs(Vs[k]).max() == 0.0:
                        continue
                    self.assertLess(abs(abs(Vs[k][REG_BUS]) - vset), 1e-6,
                                    f"{kind}: contingency {k} left bus {REG_BUS} at "
                                    f"{abs(Vs[k][REG_BUS]):.6f} instead of {vset:.6f}")
                    nb_checked += 1
                self.assertGreater(nb_checked, 0, "no contingency was actually solved")


if __name__ == "__main__":
    unittest.main()
