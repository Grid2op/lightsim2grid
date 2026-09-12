# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Branch flows from the bus voltages, in pure torch.

The batch algorithms already compute flows in C++ (``compute_flows`` /
``compute_power_flows``), but those read the voltages the solver stored and
return numbers: nothing differentiates through them. The functions here are the
same formulas written as torch operations, so a loss on a current -- an overload
penalty, the usual reason to differentiate a powerflow at all -- reaches the
injections through :class:`BatchCPUPowerFlow`'s backward.

They are kept deliberately equal, coefficient for coefficient, to what
``BaseBatchSolverSynch::_flows_of_row`` does, and a test pins that agreement
against the C++ on the same voltages.
"""

from collections import namedtuple

__all__ = ["BranchFlows", "compute_branch_flows"]


#: what :func:`compute_branch_flows` returns, all of shape ``(n_scenarios, n_branch)``
#: with the powerlines first and the transformers after, in grid order. "or" / "ex"
#: are grid2op's names for the two ends: the origin side (side 1, the high voltage
#: side of a transformer) and the extremity side.
BranchFlows = namedtuple("BranchFlows",
                         ["p_or_mw", "q_or_mvar", "p_ex_mw", "q_ex_mvar", "i_or_a", "i_ex_a"])


def _one_side(v_here, v_there, y_here, y_there, vn_kv_here, sn_mva, connected):
    """Flows measured at one end of a branch: ``S = V . conj(y_here V + y_there V')``.

    ``connected`` is the boolean "this branch carries something in this scenario";
    where it is False the end reads exactly zero, which is what the C++ does by
    zeroing the voltage of an open side rather than by masking afterwards.
    """
    import torch

    current = y_here * v_here + y_there * v_there           # per unit
    s = v_here * torch.conj(current)                        # per unit
    zero = torch.zeros((), dtype=s.real.dtype)
    p_mw = torch.where(connected, s.real * sn_mva, zero)
    q_mvar = torch.where(connected, s.imag * sn_mva, zero)

    # |S| / (sqrt(3) . |V| . vn_kv) in kA, then A. The magnitude of the voltage at
    # the measuring end is the base; a disconnected end would put 0 there, hence the
    # guard -- the numerator is zero there anyway, and this only avoids the 0/0.
    v_kv = v_here.abs() * vn_kv_here
    v_kv = torch.where(connected & (v_kv > 0.), v_kv, torch.ones((), dtype=v_kv.dtype))
    i_a = torch.where(connected, 1000. * s.abs() * sn_mva / ((3. ** 0.5) * v_kv), zero)
    return p_mw, q_mvar, i_a


def compute_branch_flows(V, y_11, y_12, y_21, y_22, bus_or, bus_ex,
                         bus_vn_kv, sn_mva, connected=None):
    """Flows at both ends of every branch, differentiable with respect to ``V``.

    V : complex ``(n_scenarios, n_bus)``, per unit, in GRID bus numbering
    y_* : complex ``(n_branch,)``, the branch's Kron-reduced admittances
    bus_or / bus_ex : int ``(n_branch,)``, the grid bus at each end
    bus_vn_kv : real ``(n_bus,)``, nominal voltage of each bus
    connected : bool ``(n_scenarios, n_branch)`` or None (everything connected)
    """
    import torch

    v_or = V[:, bus_or]
    v_ex = V[:, bus_ex]
    if connected is None:
        connected = torch.ones(v_or.shape, dtype=torch.bool)

    p_or, q_or, i_or = _one_side(v_or, v_ex, y_11, y_12, bus_vn_kv[bus_or], sn_mva, connected)
    p_ex, q_ex, i_ex = _one_side(v_ex, v_or, y_22, y_21, bus_vn_kv[bus_ex], sn_mva, connected)
    return BranchFlows(p_or_mw=p_or, q_or_mvar=q_or, p_ex_mw=p_ex, q_ex_mvar=q_ex,
                       i_or_a=i_or, i_ex_a=i_ex)
