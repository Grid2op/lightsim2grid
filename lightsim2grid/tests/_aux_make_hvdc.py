# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Shared helpers to build small grids carrying an HVDC line for the HVDC test suite.

This module is intentionally NOT named ``test_*`` so that ``unittest`` discovery
ignores it.  The test files add the tests directory to ``sys.path`` and import it
directly (see the header of each ``test_hvdc_*.py``).
"""

import warnings
import numpy as np
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower

# the angle-droop slope is stored as MW / radian internally, while the user gives MW / deg
DEG2RAD = np.pi / 180.0
MW_PER_DEG_TO_MW_PER_RAD = 180.0 / np.pi

# a "very large" reactive limit so regulating VSC never hits its bounds in the tests
_BIG_Q = 1.0e6


def make_case14_hvdc(b1=3, b2=9, *,
                     type1=0, type2=0,
                     loss_factor1=0.0, loss_factor2=0.0,
                     vreg1=False, vreg2=False,
                     vm1=1.0, vm2=1.0,
                     q1=0.0, q2=0.0,
                     power_factor1=1.0, power_factor2=1.0,
                     converters_mode=0,
                     p_setpoint=20.0,
                     r_ohm=0.0, nominal_v=0.0,
                     droop_enabled=False, p0=10.0, droop_mw_per_deg=5.0,
                     pmax12=300.0, pmax21=300.0):
    """Build a ``case14`` grid with a single HVDC line wired through ``init_hvdc_lines``.

    Returns the pandapower ``net`` (untouched, for bus count / references) and the
    lightsim2grid ``model`` with the HVDC line already declared.  ``type1/type2`` are
    ``0`` for VSC and ``1`` for LCC.
    """
    net = pn.case14()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        model = init_from_pandapower(net)
    model.init_hvdc_lines(
        np.array([b1], dtype=np.int32), np.array([b2], dtype=np.int32),
        [int(type1)], [int(type2)],
        np.array([loss_factor1]), np.array([loss_factor2]),
        [bool(vreg1)], [bool(vreg2)],
        np.array([vm1]), np.array([vm2]),
        np.array([q1]), np.array([q2]),
        np.array([-_BIG_Q]), np.array([_BIG_Q]),
        np.array([-_BIG_Q]), np.array([_BIG_Q]),
        np.array([power_factor1]), np.array([power_factor2]),
        [int(converters_mode)],
        np.array([float(p_setpoint)]),
        np.array([float(r_ohm)]), np.array([float(nominal_v)]),
        [bool(droop_enabled)],
        np.array([float(p0)]), np.array([float(droop_mw_per_deg)]),
        np.array([float(pmax12)]), np.array([float(pmax21)]),
    )
    model.tell_solver_need_reset()
    return net, model


def recv_mw(p_abs, lf_ctrl, lf_recv, r_pu=0.0, fixed_loss_mw=0.0):
    """Received power of the HVDC loss chain (mirrors ``HvdcLineContainer::recv_mw``)."""
    x = (1.0 - lf_ctrl) * p_abs
    return (1.0 - lf_recv) * (x - r_pu * x * x) - fixed_loss_mw


def droop_raw_mw(p0_mw, droop_mw_per_deg, theta1_deg, theta2_deg):
    """``raw = p0 + k * (theta1 - theta2)`` with ``k`` in MW / radian."""
    k = droop_mw_per_deg * MW_PER_DEG_TO_MW_PER_RAD
    return p0_mw + k * np.deg2rad(theta1_deg - theta2_deg)
