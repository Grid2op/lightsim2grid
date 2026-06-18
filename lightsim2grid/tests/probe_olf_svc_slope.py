# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Phase 0 probe #1 (documented, throwaway-but-kept).

Determines, empirically against PowSyBl OpenLoadFlow (OLF), the SVC voltage/slope
("droop") convention so the C++ ``VoltageControl`` extension uses the right sign:

    is OLF satisfying   Vm(reg) = Vset - s.Q   or   Vm(reg) = Vset + s.Q   ?
    in which Q sign convention (generator-injection vs load-absorption)?
    and is   s_pu = slope_kV_per_MVar * sn_mva / vn_kv(reg bus)   numerically right?

Run:  venv_debug_pypo/bin/python probe_olf_svc_slope.py
"""

import numpy as np
import pypowsybl as pp


SN_MVA = 100.0  # the per-unit base lightsim2grid will use


def build(slope_kv_per_mvar, target_v_kv=398.0, nominal_v=400.0, load_q=40.0):
    n = pp.network.create_empty()
    n.create_substations(id=["S1", "S2"])
    n.create_voltage_levels(id=["VL1", "VL2"], substation_id=["S1", "S2"],
                            topology_kind=["BUS_BREAKER"] * 2, nominal_v=[nominal_v, nominal_v])
    n.create_buses(id=["B1", "B2"], voltage_level_id=["VL1", "VL2"])
    n.create_lines(id="L", voltage_level1_id="VL1", voltage_level2_id="VL2",
                   bus1_id="B1", bus2_id="B2", r=1.0, x=20.0, g1=0, b1=0, g2=0, b2=0)
    n.create_generators(id="G", voltage_level_id="VL1", bus_id="B1", target_p=50.0,
                        target_q=0.0, target_v=nominal_v, voltage_regulator_on=True,
                        max_p=1000.0, min_p=0.0)
    n.create_loads(id="LD", voltage_level_id="VL2", bus_id="B2", p0=80.0, q0=load_q)
    n.create_static_var_compensators(id="SVC1", voltage_level_id="VL2", bus_id="B2",
                                     regulation_mode="VOLTAGE", target_v=target_v_kv,
                                     b_min=-1.0, b_max=1.0, regulating=True)
    if slope_kv_per_mvar != 0.0:
        n.create_extensions("voltagePerReactivePowerControl", id="SVC1",
                            slope=slope_kv_per_mvar)
    return n


def run(slope_kv_per_mvar):
    n = build(slope_kv_per_mvar)
    params = pp.loadflow.Parameters(
        distributed_slack=False, use_reactive_limits=False,
        provider_parameters={"voltagePerReactivePowerControl": "true"})
    res = pp.loadflow.run_ac(n, parameters=params)
    assert res[0].status == pp.loadflow.ComponentStatus.CONVERGED, res[0].status
    n.per_unit = False
    svc = n.get_static_var_compensators().loc["SVC1"]
    bus2 = n.get_buses()
    # regulated bus is B2 -> find its v_mag
    vl = n.get_voltage_levels()
    # bus row whose voltage_level_id == VL2
    b2row = bus2[bus2["voltage_level_id"] == "VL2"].iloc[0]
    vn = vl.loc["VL2", "nominal_v"]
    vm_pu = b2row["v_mag"] / vn
    q_mvar = svc["q"]              # pypowsybl terminal reactive power (load sign at terminal)
    target_v_pu = svc["target_v"] / vn
    s_pu = slope_kv_per_mvar * SN_MVA / vn
    q_pu = q_mvar / SN_MVA
    print(f"\n--- slope = {slope_kv_per_mvar} kV/MVar ---")
    print(f"  Vm(reg)      = {vm_pu:.8f} pu")
    print(f"  Vset         = {target_v_pu:.8f} pu")
    print(f"  Q(svc, pypo) = {q_mvar:.6f} MVar  ({q_pu:.8f} pu)")
    print(f"  s_pu         = {s_pu:.8f}")
    print(f"  Vset - s*Q_pypo = {target_v_pu - s_pu * q_pu:.8f}")
    print(f"  Vset + s*Q_pypo = {target_v_pu + s_pu * q_pu:.8f}")
    print(f"  Vset - s*(-Q)   = {target_v_pu - s_pu * (-q_pu):.8f}")
    print("  => the line matching Vm(reg) tells the sign AND the Q convention")
    return vm_pu, target_v_pu, s_pu, q_pu


if __name__ == "__main__":
    print("baseline (no slope): Vm should equal Vset exactly")
    run(0.0)
    run(0.05)
    run(0.10)
