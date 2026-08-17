# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Phase 0 probe #2 (documented, throwaway-but-kept).

Determines OLF's reactive-power sharing rule when N>1 generators regulate the SAME
(remote) bus, so the Python importer can compute the bordered sharing keys ``w_i``.

The bordered sharing rows are homogeneous: ``w_1.Q_{k+1} - w_{k+1}.Q_1 = 0`` -> they
encode pure proportional sharing ``Q_i / w_i = const`` with NO offset. Q_EQUAL_PROPORTION
should match with ``w_i = qmax_i - qmin_i``; K_EQUAL_PROPORTION carries a (qmin+qmax)/2
offset the homogeneous form cannot represent. This probe confirms which mode to pin.

Run:  venv_debug_pypo/bin/python probe_olf_qsharing.py
"""

import numpy as np
import pypowsybl as pp


def build(ranges):
    # two gens at DIFFERENT buses (B1, B3) both regulating remote bus B2, with
    # different reactive ranges; a load at B2 forces nonzero reactive injection.
    n = pp.network.create_empty()
    n.create_substations(id=["S1", "S2", "S3"])
    n.create_voltage_levels(id=["VL1", "VL2", "VL3"], substation_id=["S1", "S2", "S3"],
                            topology_kind=["BUS_BREAKER"] * 3, nominal_v=[400.0] * 3)
    n.create_buses(id=["B1", "B2", "B3"], voltage_level_id=["VL1", "VL2", "VL3"])
    n.create_lines(id=["L12", "L32"],
                   voltage_level1_id=["VL1", "VL3"], voltage_level2_id=["VL2", "VL2"],
                   bus1_id=["B1", "B3"], bus2_id=["B2", "B2"],
                   r=[1.0, 1.0], x=[20.0, 20.0], g1=[0, 0], b1=[0, 0], g2=[0, 0], b2=[0, 0])
    n.create_loads(id="LD", voltage_level_id="VL2", bus_id="B2", p0=120.0, q0=60.0)
    # busbar section B2 must be the regulated element. Use it as regulated_element_id.
    (qmin1, qmax1), (qmin2, qmax2) = ranges
    # slack generator G1 at B1 also regulates B2 (remote)
    n.create_generators(id="G1", voltage_level_id="VL1", bus_id="B1", target_p=60.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True,
                        max_p=1000.0, min_p=0.0)
    n.create_generators(id="G2", voltage_level_id="VL3", bus_id="B3", target_p=60.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True,
                        max_p=1000.0, min_p=0.0)
    n.create_minmax_reactive_limits(id=["G1", "G2"],
                                    min_q=[qmin1, qmin2], max_q=[qmax1, qmax2])
    # point both gens' regulating terminal at B2 (remote control)
    n.update_generators(id=["G1", "G2"], regulated_element_id=["LD", "LD"])
    return n, (qmax1 - qmin1, qmax2 - qmin2)


def run(mode, ranges):
    n, (r1, r2) = build(ranges)
    params = pp.loadflow.Parameters(
        distributed_slack=False, use_reactive_limits=False,
        provider_parameters={"generatorReactivePowerRemoteControl": "true",
                            "reactivePowerDispatchMode": mode})
    res = pp.loadflow.run_ac(n, parameters=params)
    status = res[0].status
    n.per_unit = False
    gens = n.get_generators()
    q1 = gens.loc["G1", "q"]   # pypo terminal q (load sign); injection = -q
    q2 = gens.loc["G2", "q"]
    qi1, qi2 = -q1, -q2        # generator injection convention
    print(f"\n=== mode={mode}  ranges r1={r1}, r2={r2}  status={status} ===")
    print(f"  Q_inj G1 = {qi1:.6f} MVar, G2 = {qi2:.6f} MVar")
    print(f"  Q1/Q2            = {qi1/qi2:.6f}")
    print(f"  r1/r2            = {r1/r2:.6f}   (proportional-to-range hypothesis)")
    print(f"  Q1/r1 = {qi1/r1:.6f}, Q2/r2 = {qi2/r2:.6f}   (equal => pure proportional)")


if __name__ == "__main__":
    # asymmetric ranges so a proportional rule is distinguishable from equal split
    ranges = [(-50.0, 50.0), (-100.0, 200.0)]   # r1=100, r2=300
    run("Q_EQUAL_PROPORTION", ranges)
    run("K_EQUAL_PROPORTION", ranges)
