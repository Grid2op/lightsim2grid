# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The backward of a powerflow, from the loss's sensitivity to V to each input's gradient.

Shared by :class:`BatchCPUPowerFlow` (a batch, one row per scenario) and
:class:`GridPowerFlow` (one grid, one solve, seen here as a batch of one row): what
differs between them is only who solves ``J^T lambda = xbar`` and who knows the solve,
which is the ``solver`` argument below. Every array is ``(n_rows, ...)``.
"""

import numpy as np


def cotangent(torch, V, grad_V, theta_col, vm_col, dim_J):
    """Project dL/dV onto the Newton-Raphson unknowns: an angle unknown takes
    ``Im(conj(V) . gV)``, a magnitude one ``Re(conj(V/|V|) . gV)`` -- the two halves of
    dL/dV in polar coordinates.

    Returns ``(xbar, gV, Vs)``: the cotangent, and the loss's sensitivity and the
    voltages with the rows that were never solved (NaN) neutralised -- such a row
    carries no cotangent, and must not put a NaN into the system either.
    """
    n_rows = V.shape[0]
    finite = torch.isfinite(V.real) & torch.isfinite(V.imag)
    gV = torch.where(finite, grad_V, torch.zeros_like(grad_V))
    Vs = torch.where(finite, V, torch.ones_like(V))

    xbar = torch.zeros(n_rows, dim_J, dtype=torch.float64)
    th_bus = np.flatnonzero(theta_col >= 0)
    vm_bus = np.flatnonzero(vm_col >= 0)
    if th_bus.size:
        xbar[:, torch.as_tensor(theta_col[th_bus])] = (torch.conj(Vs[:, th_bus]) * gV[:, th_bus]).imag
    if vm_bus.size:
        v_hat = Vs[:, vm_bus] / Vs[:, vm_bus].abs()
        xbar[:, torch.as_tensor(vm_col[vm_bus])] = (torch.conj(v_hat) * gV[:, vm_bus]).real
    return xbar, gV, Vs


def element_grad(torch, lam, bus_of_element, row_of_bus, sign, sn_mva):
    """lambda is the gradient with respect to the per-unit bus injection; an element's
    own is that, times how it enters its bus (``sign``: -1 for a load, +1 for a source)."""
    n_rows = lam.shape[0]
    grad = torch.zeros(n_rows, bus_of_element.shape[0], dtype=torch.float64)
    rows = np.where(bus_of_element >= 0, row_of_bus[np.clip(bus_of_element, 0, None)], -1)
    sel = np.flatnonzero(rows >= 0)   # an element off the grid, or on a bus
    if sel.size:                      # with no such equation, keeps its zero
        grad[:, torch.as_tensor(sel)] = sign * lam[:, torch.as_tensor(rows[sel])] / sn_mva
    return grad


def gen_v_grad(torch, indirect, target, share, Vs, gV):
    """The whole ``gen_v`` gradient from its two halves.

    gen_v is not an injection: it FIXES the magnitude of the bus it regulates -- through
    the starting voltage, or through the identity row of a bus pinned PV -- so lambda is
    not its gradient, and the two halves are put together here:

      indirect  ``-lambda^T dF/dv``, from the C++ that knows the solve (the dS/d|V|
                column of that bus, or lambda at a voltage-control group's row);
      direct    the loss's own dependence on ``V_k = v_k . exp(j.theta_k)`` at fixed
                unknowns, ``Re(conj(V_k/|V_k|) . gV_k)`` -- the very quantity a PQ bus
                contributes to xbar. A bus whose magnitude is an unknown hands it to
                the adjoint; a bus whose magnitude is a set-point hands it to that
                set-point.

    ``share`` is ``(n_rows, n_gen)``: where several generators regulate one bus their
    set-points are tied, so each holds a share of the one derivative that exists, and
    a generator whose set-point reaches nothing in a row (its bus switched PQ there)
    holds none. The indirect half already carries it.
    """
    grad = indirect.clone()
    live = np.flatnonzero(target >= 0)   # a generator whose set-point never
    if live.size:                        # reaches the solve keeps its zero
        bus = torch.as_tensor(target[live])
        v_hat = Vs[:, bus] / Vs[:, bus].abs()
        grad[:, torch.as_tensor(live)] += (torch.conj(v_hat) * gV[:, bus]).real * share[:, torch.as_tensor(live)]
    return grad


def input_grads(torch, solver, V, grad_V, needs, sn_mva, load_bus, gen_bus, sgen_bus):
    """Every input's gradient, ``(load_p, load_q, gen_p, sgen_p, gen_v)``, None where
    ``needs`` says the forward was not given it.

    ``solver`` knows the solve(s): ``dim_J``, the four bus maps ``theta_col``,
    ``vm_col``, ``p_row``, ``q_row`` (numpy, grid buses), ``solve_JT(xbar)``,
    ``gen_v_indirect_grad(lam)``, ``gen_p_indirect_grad(lam)``, ``gen_v_target_bus``
    (numpy) and ``gen_v_share()`` (``(n_rows, n_gen)``), all on ``(n_rows, ...)``
    numpy arrays.
    """
    xbar, gV, Vs = cotangent(torch, V, grad_V, solver.theta_col, solver.vm_col, solver.dim_J)
    lam_np = solver.solve_JT(np.ascontiguousarray(xbar.numpy()))
    lam = torch.as_tensor(lam_np)

    wants_lp, wants_lq, wants_gp, wants_sp, wants_gv = needs
    grad_load_p = element_grad(torch, lam, load_bus, solver.p_row, -1., sn_mva) if wants_lp else None
    grad_load_q = element_grad(torch, lam, load_bus, solver.q_row, -1., sn_mva) if wants_lq else None
    grad_gen_p = None
    if wants_gp:
        # plus what the outer loops change in how a unit's target P reaches the solve (a
        # unit frozen on its capability curve, one the distributed slack holds at a bound)
        grad_gen_p = (element_grad(torch, lam, gen_bus, solver.p_row, +1., sn_mva)
                      + torch.as_tensor(solver.gen_p_indirect_grad(lam_np)))
    grad_sgen_p = element_grad(torch, lam, sgen_bus, solver.p_row, +1., sn_mva) if wants_sp else None
    grad_gen_v = None
    if wants_gv:
        grad_gen_v = gen_v_grad(torch, torch.as_tensor(solver.gen_v_indirect_grad(lam_np)),
                                solver.gen_v_target_bus, torch.as_tensor(solver.gen_v_share()), Vs, gV)
    return grad_load_p, grad_load_q, grad_gen_p, grad_sgen_p, grad_gen_v
