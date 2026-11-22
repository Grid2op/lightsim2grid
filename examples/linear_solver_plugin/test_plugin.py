#!/usr/bin/env python3
"""Smoke-test for the linear_solver_plugin example: a linear solver (the mylu stand-in C
library, see mylu/mylu.h) brought into lightsim2grid through an algorithm plugin.

Build the plugin first (see CMakeLists.txt or MyLU_plugin.cpp), then run:
    python test_plugin.py

Requires pandapower. Checks that:
- the plugin registers NR_MyLU, NRSing_MyLU, NROuter_MyLU and DC_MyLU;
- each gives what its SparseLU counterpart gives, on a single powerflow;
- an N-1 contingency analysis with NR_MyLU, some rows of which cut buses off, gives what
  NR_SparseLU gives, while analyzing the matrix once for the whole batch (every row after that
  is a refactorization).
"""
import pathlib
import platform
import warnings

import numpy as np
import pandapower.networks as pn

import lightsim2grid
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP
from lightsim2grid.network import init_from_pandapower

TOL = 1e-10
# the pandapower conversion warns about the test cases' missing data, not about the plugin
warnings.filterwarnings("ignore", module="lightsim2grid.network")


def find_plugin():
    build = pathlib.Path(__file__).parent / "build"
    if platform.system() == "Windows":
        candidates = [build / "Release" / "linear_solver_plugin.dll", build / "linear_solver_plugin.dll"]
    else:
        candidates = [build / "liblinear_solver_plugin.so"]
    for p in candidates:
        if p.exists():
            return str(p)
    raise FileNotFoundError(f"Plugin not found (tried {[str(c) for c in candidates]}). "
                            "Build it first (see CMakeLists.txt).")


def same(a, b, what):
    a = np.asarray(a)
    b = np.asarray(b)
    assert a.shape == b.shape, f"{what}: shapes {a.shape} vs {b.shape}"
    # a bus a contingency cut off is NaN in both
    assert np.array_equal(np.isnan(a), np.isnan(b)), f"{what}: not the same NaN entries"
    ok = ~np.isnan(a)
    err = np.max(np.abs(a[ok] - b[ok])) if ok.any() else 0.
    assert err < 1e-8, f"{what}: max difference {err}"


lightsim2grid.load_algorithm_plugin(find_plugin())
print("Plugin loaded successfully.")

grid = init_from_pandapower(pn.case118())
names = grid.available_algorithm_names()
for name in ("NR_MyLU", "NRSing_MyLU", "NROuter_MyLU", "DC_MyLU"):
    assert name in names, f"{name} not in {names}"
print("Registered: NR_MyLU, NRSing_MyLU, NROuter_MyLU, DC_MyLU")

# ------------------------------------------------------------------ single powerflows
for mine, ref in (("NR_MyLU", "NR_SparseLU"), ("NRSing_MyLU", "NRSing_SparseLU"),
                  ("NROuter_MyLU", "NROuter_SparseLU")):
    results = []
    for name in (mine, ref):
        g = init_from_pandapower(pn.case118())
        g.change_algorithm(name)
        V = g.ac_pf(np.ones(g.total_bus(), dtype=complex), 30, TOL)
        assert V.shape[0] > 0, f"{name} did not converge: {g.get_algo().get_error()}"
        results.append(V)
        if name == mine:
            stats = g.get_algo().get_linear_solver_stats()
            assert stats.nb_analyze == 1, f"{name}: {stats.nb_analyze} analyze"
    same(results[0], results[1], f"{mine} vs {ref}")
    print(f"{mine}: same voltages as {ref}")

results = []
for name in ("DC_MyLU", "DC_SparseLU"):
    g = init_from_pandapower(pn.case118())
    g.change_algorithm(name)
    V = g.dc_pf(np.ones(g.total_bus(), dtype=complex), 10, TOL)
    assert V.shape[0] > 0, f"{name} failed"
    results.append(V)
same(results[0], results[1], "DC_MyLU vs DC_SparseLU")
print("DC_MyLU: same voltages as DC_SparseLU")

# ------------------------------------------------------------------ N-1
# Some contingencies of case118 cut buses off: with handle_disconnected_grid, those buses are
# masked (their rows of the Jacobian rewritten by value) and the matrix keeps its pattern.
results = []
for name in ("NR_MyLU", "NR_SparseLU"):
    g = init_from_pandapower(pn.case118())
    ca = ContingencyAnalysisCPP(g)
    ca.change_algorithm(name)
    ca.handle_disconnected_grid = True
    ca.add_all_n1()
    ca.compute(np.ones(g.total_bus(), dtype=complex), 30, TOL)
    assert all(ca.converged_mask()), f"{name}: a contingency did not converge"
    V = np.array(ca.get_voltages())
    results.append(V)
    if name == "NR_MyLU":
        assert (np.abs(V) == 0.).any(axis=1).any(), "no contingency cut a bus off: nothing was masked"
        stats = ca.get_linear_solver_stats()
        print(f"NR_MyLU N-1: {V.shape[0]} contingencies, {stats.nb_analyze} analyze, "
              f"{stats.nb_factorize} factorize, {stats.nb_refactorize} refactorize")
        assert stats.nb_analyze == 1, f"the pattern changed during the batch: {stats.nb_analyze} analyze"
        assert stats.nb_refactorize > 0
same(results[0], results[1], "N-1 NR_MyLU vs NR_SparseLU")
print("NR_MyLU N-1: same voltages as NR_SparseLU")
print("All checks passed.")
