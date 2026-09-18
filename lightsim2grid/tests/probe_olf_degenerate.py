# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Phase 0 probe #3 (documented, throwaway-but-kept).

Checks what OLF does for the configurations that make the bordered Jacobian singular,
so we know which ones lightsim2grid must reject (clear std::runtime_error) and which,
if any, OLF resolves with a cheap deterministic fallback worth mirroring:

  (a) a generator regulating the SLACK bus (regulated bus has no Vm unknown);
  (b) a generator regulating a bus that is already PV via a local generator;
  (c) a sloped SVC sharing its regulated bus with another controller.

Run:  venv_debug_pypo/bin/python probe_olf_degenerate.py
"""

import pypowsybl as pp


def _base():
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
    return n


def _params():
    return pp.loadflow.Parameters(
        distributed_slack=False, use_reactive_limits=False,
        provider_parameters={"generatorReactivePowerRemoteControl": "true",
                            "voltagePerReactivePowerControl": "true"})


def case_a_regulate_slack():
    # G2 (at B3) regulates B1, which carries the slack generator G1 -> reg bus is slack
    n = _base()
    n.create_generators(id="G1", voltage_level_id="VL1", bus_id="B1", target_p=60.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True,
                        max_p=1000.0, min_p=0.0)
    n.create_generators(id="G2", voltage_level_id="VL3", bus_id="B3", target_p=60.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True,
                        max_p=1000.0, min_p=0.0)
    n.update_generators(id="G2", regulated_element_id="G1")  # regulate the slack-bus element
    res = pp.loadflow.run_ac(n, parameters=_params(), )
    return res, n


def case_b_regulate_existing_pv():
    # G2 (remote at B3) regulates B2, which already has a LOCAL voltage-regulating gen GL
    n = _base()
    n.create_generators(id="G1", voltage_level_id="VL1", bus_id="B1", target_p=60.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True,
                        max_p=1000.0, min_p=0.0)
    n.create_generators(id="GL", voltage_level_id="VL2", bus_id="B2", target_p=10.0,
                        target_q=0.0, target_v=399.0, voltage_regulator_on=True,
                        max_p=1000.0, min_p=0.0)
    n.create_generators(id="G2", voltage_level_id="VL3", bus_id="B3", target_p=60.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True,
                        max_p=1000.0, min_p=0.0)
    n.update_generators(id="G2", regulated_element_id="LD")  # regulate B2 remotely
    res = pp.loadflow.run_ac(n, parameters=_params())
    return res, n


def case_c_sloped_svc_shared():
    # sloped SVC at B2 and a remote gen G2 both regulate B2
    n = _base()
    n.create_generators(id="G1", voltage_level_id="VL1", bus_id="B1", target_p=60.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True,
                        max_p=1000.0, min_p=0.0)
    n.create_generators(id="G2", voltage_level_id="VL3", bus_id="B3", target_p=60.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True,
                        max_p=1000.0, min_p=0.0)
    n.update_generators(id="G2", regulated_element_id="LD")
    n.create_static_var_compensators(id="SVC1", voltage_level_id="VL2", bus_id="B2",
                                     regulation_mode="VOLTAGE", target_v=400.0,
                                     b_min=-1.0, b_max=1.0, regulating=True)
    n.create_extensions("voltagePerReactivePowerControl", id="SVC1", slope=0.05)
    res = pp.loadflow.run_ac(n, parameters=_params())
    return res, n


if __name__ == "__main__":
    for name, fn in [("(a) regulate slack bus", case_a_regulate_slack),
                     ("(b) regulate already-PV bus", case_b_regulate_existing_pv),
                     ("(c) sloped SVC sharing reg bus", case_c_sloped_svc_shared)]:
        try:
            res, n = fn()
            st = res[0].status
            n.per_unit = False
            print(f"\n{name}: status={st}")
            gens = n.get_generators()
            print("  gen q (terminal):")
            print(gens["q"].to_string())
            try:
                print("  svc q:", n.get_static_var_compensators()["q"].to_string())
            except Exception:
                pass
        except Exception as e:
            print(f"\n{name}: RAISED {type(e).__name__}: {e}")
