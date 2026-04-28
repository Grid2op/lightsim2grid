# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Tests for NRAlgo scaling and refactor policies:
  - ScalingPolicyType / RefactorPolicyType enum binding
  - Per-policy parameter getters/setters
  - AlgoConfig round-trip (get_config -> set_config -> get_config)
  - Convergence with every (scaling policy, refactor policy) combination
  - GridModel.get_ac_algo_config / set_ac_algo_config
"""

import os
import unittest
import numpy as np
import zipfile
from scipy import sparse

from lightsim2grid.lightsim2grid_cpp import (
    ScalingPolicyType,
    RefactorPolicyType,
    AlgoConfig,
    SparseLUSolver,
)

try:
    import lightsim2grid
    from lightsim2grid.gridmodel import init as init_gridmodel
    _GRIDMODEL_AVAILABLE = True
except Exception:
    _GRIDMODEL_AVAILABLE = False

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_CASE14_ZIP = os.path.join(_THIS_DIR, "case14.zip")


def _load_case(zip_path):
    """Return (Ybus, Sbus, V_init, pv, pq, ref, slack_weights)."""
    def _arr(z, name):
        path = z.extract(name)
        arr = np.load(path)
        os.remove(path)
        return arr

    with zipfile.ZipFile(zip_path) as z:
        V_init = _arr(z, "V0.npy")
        pq     = _arr(z, "pq.npy")
        pv     = _arr(z, "pv.npy")
        Sbus   = _arr(z, "Sbus.npy")
        Ybus   = sparse.csc_matrix(_arr(z, "Ybus.npy"))

    ref = np.array(sorted(set(range(Sbus.shape[0])) - set(pv) - set(pq)), dtype=np.int32)
    slack_weights = np.zeros(Sbus.shape[0])
    slack_weights[ref] = 1.0 / ref.shape[0]
    return Ybus, Sbus, V_init, pv, pq, ref, slack_weights


class TestEnumBinding(unittest.TestCase):
    """ScalingPolicyType and RefactorPolicyType are importable with expected values."""

    def test_scaling_policy_values(self):
        self.assertEqual(int(ScalingPolicyType.NoScaling),        0)
        self.assertEqual(int(ScalingPolicyType.MaxVoltageChange), 1)
        self.assertEqual(int(ScalingPolicyType.LineSearch),       2)
        self.assertEqual(int(ScalingPolicyType.Iwamoto),          3)

    def test_refactor_policy_values(self):
        self.assertEqual(int(RefactorPolicyType.AlwaysRefactor), 0)
        self.assertEqual(int(RefactorPolicyType.EveryN),         1)
        self.assertEqual(int(RefactorPolicyType.Chord),          2)

    def test_scaling_members_complete(self):
        expected = {"NoScaling", "MaxVoltageChange", "LineSearch", "Iwamoto"}
        self.assertEqual(set(ScalingPolicyType.__members__), expected)

    def test_refactor_members_complete(self):
        expected = {"AlwaysRefactor", "EveryN", "Chord"}
        self.assertEqual(set(RefactorPolicyType.__members__), expected)


class TestDefaultPolicyState(unittest.TestCase):
    """Fresh SparseLUSolver has expected default policy settings."""

    def setUp(self):
        self.s = SparseLUSolver()

    def test_default_scaling_policy(self):
        self.assertEqual(self.s.get_scaling_policy_type(), ScalingPolicyType.NoScaling)

    def test_default_refactor_policy(self):
        self.assertEqual(self.s.get_refactor_policy(), RefactorPolicyType.AlwaysRefactor)

    def test_default_max_dVa(self):
        self.assertAlmostEqual(self.s.get_max_dVa(), 0.5)

    def test_default_max_dVm(self):
        self.assertAlmostEqual(self.s.get_max_dVm(), 0.1)

    def test_default_ls_c(self):
        self.assertAlmostEqual(self.s.get_ls_c(), 1e-4, places=10)

    def test_default_ls_rho(self):
        self.assertAlmostEqual(self.s.get_ls_rho(), 0.5)

    def test_default_ls_max_iter(self):
        self.assertEqual(self.s.get_ls_max_iter(), 10)

    def test_default_iw_mu_min(self):
        self.assertAlmostEqual(self.s.get_iw_mu_min(), 1e-4, places=10)

    def test_default_iw_mu_max(self):
        self.assertAlmostEqual(self.s.get_iw_mu_max(), 1.0)

    def test_default_refactor_every_n(self):
        self.assertEqual(self.s.get_refactor_every_n(), 4)


class TestPolicySetters(unittest.TestCase):
    """All parameter setters round-trip correctly."""

    def setUp(self):
        self.s = SparseLUSolver()

    def test_set_scaling_policy_all_values(self):
        for pol in ScalingPolicyType.__members__.values():
            self.s.set_scaling_policy(pol)
            self.assertEqual(self.s.get_scaling_policy_type(), pol)

    def test_set_refactor_policy_all_values(self):
        for pol in RefactorPolicyType.__members__.values():
            self.s.set_refactor_policy(pol)
            self.assertEqual(self.s.get_refactor_policy(), pol)

    def test_max_dVa_setter(self):
        self.s.set_max_dVa(0.3)
        self.assertAlmostEqual(self.s.get_max_dVa(), 0.3)

    def test_max_dVm_setter(self):
        self.s.set_max_dVm(0.05)
        self.assertAlmostEqual(self.s.get_max_dVm(), 0.05)

    def test_ls_c_setter(self):
        self.s.set_ls_c(1e-3)
        self.assertAlmostEqual(self.s.get_ls_c(), 1e-3, places=10)

    def test_ls_rho_setter(self):
        self.s.set_ls_rho(0.7)
        self.assertAlmostEqual(self.s.get_ls_rho(), 0.7)

    def test_ls_max_iter_setter(self):
        self.s.set_ls_max_iter(30)
        self.assertEqual(self.s.get_ls_max_iter(), 30)

    def test_iw_mu_min_setter(self):
        self.s.set_iw_mu_min(1e-3)
        self.assertAlmostEqual(self.s.get_iw_mu_min(), 1e-3, places=10)

    def test_iw_mu_max_setter(self):
        self.s.set_iw_mu_max(0.9)
        self.assertAlmostEqual(self.s.get_iw_mu_max(), 0.9)

    def test_refactor_every_n_setter(self):
        self.s.set_refactor_every_n(3)
        self.assertEqual(self.s.get_refactor_every_n(), 3)

    def test_set_scaling_policy_preserves_params(self):
        # Params stored on the solver survive a policy switch
        self.s.set_iw_mu_min(1e-3)
        self.s.set_scaling_policy(ScalingPolicyType.Iwamoto)
        self.assertAlmostEqual(self.s.get_iw_mu_min(), 1e-3, places=10)


class TestAlgoConfigRoundTrip(unittest.TestCase):
    """AlgoConfig get_config/set_config round-trips all parameters."""

    def setUp(self):
        self.s = SparseLUSolver()

    def _configure(self):
        self.s.set_scaling_policy(ScalingPolicyType.LineSearch)
        self.s.set_refactor_policy(RefactorPolicyType.EveryN)
        self.s.set_max_dVa(0.3)
        self.s.set_max_dVm(0.05)
        self.s.set_ls_c(2e-4)
        self.s.set_ls_rho(0.6)
        self.s.set_ls_max_iter(15)
        self.s.set_iw_mu_min(5e-4)
        self.s.set_iw_mu_max(0.95)
        self.s.set_refactor_every_n(3)

    def test_get_config_returns_AlgoConfig(self):
        cfg = self.s.get_config()
        self.assertIsInstance(cfg, AlgoConfig)

    def test_config_int_params_length(self):
        cfg = self.s.get_config()
        self.assertEqual(len(cfg.int_params), 4)

    def test_config_real_params_length(self):
        cfg = self.s.get_config()
        self.assertEqual(len(cfg.real_params), 6)

    def test_round_trip_scaling_policy(self):
        self._configure()
        cfg = self.s.get_config()
        s2 = SparseLUSolver()
        s2.set_config(cfg)
        self.assertEqual(s2.get_scaling_policy_type(), ScalingPolicyType.LineSearch)

    def test_round_trip_refactor_policy(self):
        self._configure()
        cfg = self.s.get_config()
        s2 = SparseLUSolver()
        s2.set_config(cfg)
        self.assertEqual(s2.get_refactor_policy(), RefactorPolicyType.EveryN)

    def test_round_trip_max_dVa(self):
        self._configure()
        cfg = self.s.get_config()
        s2 = SparseLUSolver()
        s2.set_config(cfg)
        self.assertAlmostEqual(s2.get_max_dVa(), 0.3)

    def test_round_trip_max_dVm(self):
        self._configure()
        cfg = self.s.get_config()
        s2 = SparseLUSolver()
        s2.set_config(cfg)
        self.assertAlmostEqual(s2.get_max_dVm(), 0.05)

    def test_round_trip_ls_c(self):
        self._configure()
        cfg = self.s.get_config()
        s2 = SparseLUSolver()
        s2.set_config(cfg)
        self.assertAlmostEqual(s2.get_ls_c(), 2e-4, places=10)

    def test_round_trip_ls_rho(self):
        self._configure()
        cfg = self.s.get_config()
        s2 = SparseLUSolver()
        s2.set_config(cfg)
        self.assertAlmostEqual(s2.get_ls_rho(), 0.6)

    def test_round_trip_ls_max_iter(self):
        self._configure()
        cfg = self.s.get_config()
        s2 = SparseLUSolver()
        s2.set_config(cfg)
        self.assertEqual(s2.get_ls_max_iter(), 15)

    def test_round_trip_iw_mu_min(self):
        self._configure()
        cfg = self.s.get_config()
        s2 = SparseLUSolver()
        s2.set_config(cfg)
        self.assertAlmostEqual(s2.get_iw_mu_min(), 5e-4, places=10)

    def test_round_trip_iw_mu_max(self):
        self._configure()
        cfg = self.s.get_config()
        s2 = SparseLUSolver()
        s2.set_config(cfg)
        self.assertAlmostEqual(s2.get_iw_mu_max(), 0.95)

    def test_round_trip_refactor_every_n(self):
        self._configure()
        cfg = self.s.get_config()
        s2 = SparseLUSolver()
        s2.set_config(cfg)
        self.assertEqual(s2.get_refactor_every_n(), 3)

    def test_config_int_params_encoding(self):
        self._configure()
        cfg = self.s.get_config()
        # int_params: [ScalingPolicyType, RefactorPolicyType, ls_max_iter, refactor_every_n]
        self.assertEqual(cfg.int_params[0], int(ScalingPolicyType.LineSearch))
        self.assertEqual(cfg.int_params[1], int(RefactorPolicyType.EveryN))
        self.assertEqual(cfg.int_params[2], 15)
        self.assertEqual(cfg.int_params[3], 3)

    def test_config_real_params_encoding(self):
        self._configure()
        cfg = self.s.get_config()
        # real_params: [max_dVa, max_dVm, ls_c, ls_rho, iw_mu_min, iw_mu_max]
        self.assertAlmostEqual(cfg.real_params[0], 0.3)
        self.assertAlmostEqual(cfg.real_params[1], 0.05)
        self.assertAlmostEqual(cfg.real_params[2], 2e-4, places=10)
        self.assertAlmostEqual(cfg.real_params[3], 0.6)
        self.assertAlmostEqual(cfg.real_params[4], 5e-4, places=10)
        self.assertAlmostEqual(cfg.real_params[5], 0.95)


class TestConvergenceWithPolicies(unittest.TestCase):
    """
    Verify that all scaling policies × all refactor policies converge on case14.
    The reference solution (Va/Vm) must match the default NoScaling/AlwaysRefactor result.
    """

    TOL_SOLVER = 1e-8
    TOL_RESULT = 1e-5
    MAX_IT     = 30

    @classmethod
    def setUpClass(cls):
        if not os.path.isfile(_CASE14_ZIP):
            raise unittest.SkipTest("case14.zip not found")
        cls.Ybus, cls.Sbus, cls.V_init, cls.pv, cls.pq, cls.ref, cls.sw = \
            _load_case(_CASE14_ZIP)

        # reference: default policy
        ref_solver = SparseLUSolver()
        ok = ref_solver.compute_pf(
            cls.Ybus, cls.V_init.copy(), cls.Sbus,
            cls.ref, cls.sw, cls.pv, cls.pq,
            cls.MAX_IT, cls.TOL_SOLVER)
        assert ok, "reference solve failed"
        cls.Va_ref = ref_solver.get_Va().copy()
        cls.Vm_ref = ref_solver.get_Vm().copy()

    def _run(self, scaling, refactor, refactor_every_n=4):
        s = SparseLUSolver()
        s.set_scaling_policy(scaling)
        s.set_refactor_policy(refactor)
        s.set_refactor_every_n(refactor_every_n)
        ok = s.compute_pf(
            self.Ybus, self.V_init.copy(), self.Sbus,
            self.ref, self.sw, self.pv, self.pq,
            self.MAX_IT, self.TOL_SOLVER)
        return ok, s

    def _assert_converged_and_matches_ref(self, scaling, refactor, **kw):
        ok, s = self._run(scaling, refactor, **kw)
        label = f"{scaling.name}/{refactor.name}"
        self.assertTrue(ok, f"{label}: solver did not converge")
        np.testing.assert_allclose(
            s.get_Va(), self.Va_ref, atol=self.TOL_RESULT,
            err_msg=f"{label}: Va differs from reference")
        np.testing.assert_allclose(
            s.get_Vm(), self.Vm_ref, atol=self.TOL_RESULT,
            err_msg=f"{label}: Vm differs from reference")

    # --- NoScaling ---
    def test_noscaling_always(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.NoScaling, RefactorPolicyType.AlwaysRefactor)

    def test_noscaling_everyn(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.NoScaling, RefactorPolicyType.EveryN, refactor_every_n=2)

    def test_noscaling_chord(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.NoScaling, RefactorPolicyType.Chord)

    # --- MaxVoltageChange ---
    def test_maxvc_always(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.MaxVoltageChange, RefactorPolicyType.AlwaysRefactor)

    def test_maxvc_everyn(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.MaxVoltageChange, RefactorPolicyType.EveryN, refactor_every_n=2)

    def test_maxvc_chord(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.MaxVoltageChange, RefactorPolicyType.Chord)

    # --- LineSearch ---
    def test_linesearch_always(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.LineSearch, RefactorPolicyType.AlwaysRefactor)

    def test_linesearch_everyn(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.LineSearch, RefactorPolicyType.EveryN, refactor_every_n=2)

    def test_linesearch_chord(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.LineSearch, RefactorPolicyType.Chord)

    # --- Iwamoto ---
    def test_iwamoto_always(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.Iwamoto, RefactorPolicyType.AlwaysRefactor)

    def test_iwamoto_everyn(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.Iwamoto, RefactorPolicyType.EveryN, refactor_every_n=2)

    def test_iwamoto_chord(self):
        self._assert_converged_and_matches_ref(
            ScalingPolicyType.Iwamoto, RefactorPolicyType.Chord)


class TestAlgoConfigViaGridModel(unittest.TestCase):
    """GridModel.get_ac_algo_config / set_ac_algo_config round-trips."""

    @unittest.skipUnless(_GRIDMODEL_AVAILABLE, "lightsim2grid gridmodel not available")
    def _make_gridmodel(self):
        import pandapower as pp
        net = pp.networks.case14()
        from lightsim2grid.gridmodel import init as gm_init
        return gm_init(net)

    def _try_import_gridmodel_direct(self):
        try:
            from lightsim2grid.lightsim2grid_cpp import GridModel
            return GridModel()
        except Exception:
            return None

    def test_get_ac_algo_config_returns_AlgoConfig(self):
        gm = self._try_import_gridmodel_direct()
        if gm is None:
            self.skipTest("GridModel C++ class not directly instantiable without a grid")
        cfg = gm.get_ac_algo_config()
        self.assertIsInstance(cfg, AlgoConfig)

    def test_get_dc_algo_config_returns_AlgoConfig(self):
        gm = self._try_import_gridmodel_direct()
        if gm is None:
            self.skipTest("GridModel C++ class not directly instantiable without a grid")
        cfg = gm.get_dc_algo_config()
        self.assertIsInstance(cfg, AlgoConfig)

    def test_set_ac_algo_config_round_trip(self):
        gm = self._try_import_gridmodel_direct()
        if gm is None:
            self.skipTest("GridModel C++ class not directly instantiable without a grid")

        # Build a custom config
        cfg = AlgoConfig()
        cfg.int_params  = [int(ScalingPolicyType.Iwamoto),
                           int(RefactorPolicyType.EveryN),
                           20, 3]
        cfg.real_params = [0.4, 0.08, 1e-3, 0.7, 5e-4, 0.9]

        gm.set_ac_algo_config(cfg)
        out = gm.get_ac_algo_config()

        self.assertEqual(out.int_params[0], int(ScalingPolicyType.Iwamoto))
        self.assertEqual(out.int_params[1], int(RefactorPolicyType.EveryN))
        self.assertEqual(out.int_params[2], 20)
        self.assertEqual(out.int_params[3], 3)
        self.assertAlmostEqual(out.real_params[0], 0.4)
        self.assertAlmostEqual(out.real_params[4], 5e-4, places=10)


if __name__ == "__main__":
    unittest.main()
