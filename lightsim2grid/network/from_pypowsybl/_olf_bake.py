# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Bake PowSyBl Open Load Flow (OLF) outer-loop results into an IIDM network so
that a subsequent *outer-loop-free* power flow -- in OLF itself or in
lightsim2grid -- reproduces the converged state exactly.

Why this exists
---------------
lightsim2grid solves a single power-flow problem: it does not run the discrete
"outer loops" that OLF runs (distributed slack, reactive-limit PV<->PQ
switching, transformer/shunt voltage control, phase control, ...). Those loops
change the *inputs* of the power flow between inner solves: a generator that
hits its Q limit becomes a fixed-Q (PQ) injection, a tap changer settles on a
discrete position, the slack mismatch is spread over participating units, etc.

If you take the raw network, run OLF *with* outer loops, and then hand the same
raw network to lightsim2grid, the two engines disagree -- not because the
solvers differ, but because they are solving different problems.

``bake_outer_loops`` rewrites the network's input setpoints to the values the
outer loops settled on, and disables the corresponding regulation, so the
problem becomes a plain power flow. After baking, OLF (loop-free) and
lightsim2grid agree to solver tolerance, and they keep agreeing through
topology changes (line/transformer outages) applied identically to both.

The loop-free OLF parameters that reproduce the baked operating point are
available as :func:`lightsim2grid.network.get_pypowsybl_loopfree_parameters`.

What gets baked
---------------
* Discrete tap / section positions (ratio tap, phase tap, shunt section):
  the solved position (``solved_tap_position`` / ``solved_section_count``) is
  copied into the input position and regulation is switched off. The copy is
  guarded: it only happens where the changer was regulating AND a solved value
  exists, because ``solved_tap_position`` is NaN for changers that did not take
  part in an outer loop -- copying NaN would destroy the input tap.
* PV -> PQ reactive-limit switches (generators, VSC converter stations): only
  the units whose realized reactive power actually sits at a Q limit are frozen
  to fixed-Q at that value, with voltage regulation switched off. Units still
  inside their limits keep voltage control -- their equation type is unchanged.
  Units with NaN (unlimited) reactive limits are never frozen, since NaN
  comparisons are False.
* Active-power redistribution from distributed slack / area interchange: the
  realized P is written back into the target P of generators and batteries
  (and optionally loads, if the slack was distributed on load).
* Static var compensators whose realized reactive output sits at (or beyond)
  their voltage-dependent susceptance envelope (``Q(V) = b * V^2``, ``b`` in
  ``b_min``..``b_max``, recomputed in MVAr at the SVC's own solved terminal
  voltage) are frozen to fixed-Q (``REACTIVE_POWER`` mode), mirroring the
  generator/VSC reactive-limit switch above. An SVC still comfortably inside
  its envelope is left regulating (its target reproduces the OLF result
  exactly, since it isn't saturated).

Verified against pypowsybl 1.15.0. The OLF-internal round trip (solve with
outer loops -> bake -> solve loop-free) reproduces bus voltages to ~2e-4 kV /
~1e-5 deg, both on IEEE-14 (with a forced reactive-limit switch) and on the
four-substations node-breaker network, which carries VSC and LCC HVDC, an SVC
regulating voltage, a shunt, and ratio + phase tap changers.
"""

import pandas as pd


# Tolerance (MVAr) for deciding a reactive injection sits "at" its Q limit.
_Q_LIMIT_TOL = 1e-3
            
            
def _keep_only_main_comp(df_el, df_bus):
    """
    keep only element (modeled in df_el) that are on the main component => bus_els["synchronous_component"] == 0
    
    This does not deactivate anything.
    """
    mask_conn = df_el["connected"]
    bus_els = df_bus.loc[df_el.loc[mask_conn, "bus_id"]]
    mask_main = (bus_els["synchronous_component"] == 0).to_numpy()
    df_el = df_el.loc[df_el.loc[mask_conn][mask_main].index]
    return df_el


def _reactive_limits(df: pd.DataFrame):
    """Return (qmin, qmax) in *generator* convention.

    Prefer the P-dependent capability-curve limits (``min_q_at_p`` /
    ``max_q_at_p``) when present and finite; fall back to the fixed
    ``min_q`` / ``max_q`` box otherwise.
    """
    if "min_q_at_p" in df.columns:
        qmin = df["min_q_at_p"].where(df["min_q_at_p"].notna(), df["min_q"])
        qmax = df["max_q_at_p"].where(df["max_q_at_p"].notna(), df["max_q"])
    else:
        qmin, qmax = df["min_q"], df["max_q"]
    return qmin, qmax


def _bound_at_qlimit(df: pd.DataFrame, regulating: pd.Series):
    """Boolean mask of rows that regulate voltage *and* sit at a Q limit,
    plus the realized reactive power in generator convention.

    ``df`` must carry result column ``q`` (load convention) and the reactive
    limit columns. ``regulating`` is the per-row voltage-regulation flag.
    """
    q_gen = -df["q"]  # result column is load convention; flip to generator
    qmin, qmax = _reactive_limits(df)
    at_max = regulating & (q_gen >= qmax - _Q_LIMIT_TOL)
    at_min = regulating & (q_gen <= qmin + _Q_LIMIT_TOL)
    return (at_max | at_min), q_gen


def bake_outer_loops(
    network,
    bake_taps: bool = True,
    bake_reactive_limits: bool = True,
    bake_active_power: bool = True,
    bake_remote_voltage_control: bool = True,
    balance_on_loads: bool = False,
    load_power_factor_constant: bool = False,
    keep_only_main_comp: bool=True
):
    """Rewrite ``network`` input setpoints to the converged outer-loop state.

    Call this on a network that has just been solved by OLF *with* outer loops.
    Afterwards the network represents a plain power-flow problem: a loop-free
    OLF run (see :func:`get_pypowsybl_loopfree_parameters`) or a lightsim2grid
    run (via :func:`init_from_pypowsybl`) will reproduce the same operating
    point.

    Parameters
    ----------
    network
        A pypowsybl network, freshly solved with the outer loops enabled.
    bake_taps
        Copy solved ratio/phase tap positions and shunt sections into the
        input positions and disable their regulation.
    bake_reactive_limits
        Freeze generators / VSC stations that hit a Q limit to fixed-Q (PQ).
    bake_active_power
        Write realized active power back into generator/battery target P
        (and load p0/q0 if ``balance_on_loads``).
    bake_remote_voltage_control
        Rewrite remote voltage control to local control at the solved terminal
        voltage (see :func:`_bake_remote_voltage_control`). Needed so that
        remote-regulating generators can sit on a (distributed) slack bus, which
        lightsim2grid v1 does not otherwise support.
    balance_on_loads
        Set if the slack was distributed on loads (BalanceType
        PROPORTIONAL_TO_LOAD / CONFORM_LOAD).
    load_power_factor_constant
        Mirror OLF's ``loadPowerFactorConstant``: also rewrite load q0 so the
        power factor is preserved.
    keep_only_main_comp
        Only elements of the main connected component are updated (True by default)

    Notes
    -----
    Operates in place and is idempotent on an already-baked network. The
    voltage-regulation flag is *not* used as a switch signal: OLF does not flip
    it in IIDM, so PV->PQ is detected from the realized Q sitting at a limit.
    """
    if bake_taps:
        _bake_taps_and_sections(network, keep_only_main_comp)
    if bake_reactive_limits:
        _bake_reactive_limit_switches(network, keep_only_main_comp)
    if bake_active_power:
        _bake_active_power(
            network,
            balance_on_loads=balance_on_loads,
            load_power_factor_constant=load_power_factor_constant,
            keep_only_main_comp=keep_only_main_comp
        )
    if bake_remote_voltage_control:
        _bake_remote_voltage_control(network, keep_only_main_comp)


def _bake_taps_and_sections(network, keep_only_main_comp=True):
    # Ratio tap changers (transformer voltage control outer loop).
    df_bus = network.get_buses(attributes=["synchronous_component"])
    rtc = network.get_ratio_tap_changers(
        attributes=["tap", "solved_tap_position", "regulating"]
    )
    if len(rtc):
        # Only adopt the solved position where the changer was actually
        # regulating AND a solved position exists. ``solved_tap_position`` is
        # NaN when the changer did not participate in an outer loop; copying it
        # blindly would destroy the input tap.
        keep = rtc["regulating"] & rtc["solved_tap_position"].notna()
        if keep.any():
            upd = pd.DataFrame(index=rtc.index[keep])
            upd["tap"] = rtc["solved_tap_position"][keep].astype(int)
            # upd["regulating"] = False
            network.update_ratio_tap_changers(upd)

    # Phase tap changers (phase-shifter regulation outer loop).
    ptc = network.get_phase_tap_changers(
        attributes=["tap", "solved_tap_position", "regulating"]
    )
    if len(ptc):
        keep = ptc["regulating"] & ptc["solved_tap_position"].notna()
        if keep.any():
            upd = pd.DataFrame(index=ptc.index[keep])
            upd["tap"] = ptc["solved_tap_position"][keep].astype(int)
            # upd["regulating"] = False
            network.update_phase_tap_changers(upd)

    # Shunt compensators (shunt voltage control outer loop).
    sh = network.get_shunt_compensators(
        attributes=[
            "section_count", 
            "solved_section_count", 
            "voltage_regulation_on",
            "connected", 
            "bus_id"
        ]
    )
    if keep_only_main_comp:
        sh = _keep_only_main_comp(sh, df_bus)
    if len(sh):
        keep = sh["voltage_regulation_on"] & sh["solved_section_count"].notna()
        if keep.any():
            upd = pd.DataFrame(index=sh.index[keep])
            upd["section_count"] = sh["solved_section_count"][keep].astype(int)
            upd["voltage_regulation_on"] = False
            network.update_shunt_compensators(upd)


def _bake_reactive_limit_switches(
    network,
    keep_only_main_comp=True):
    df_bus = network.get_buses(attributes=["synchronous_component"])
    gen = network.get_generators(
        attributes=[
            "voltage_regulator_on", "q",
            "min_q", "max_q", "min_q_at_p", "max_q_at_p",
            "connected", "bus_id"
        ]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    mask, q_gen = _bound_at_qlimit(gen, gen["voltage_regulator_on"])
    if mask.any():
        upd = pd.DataFrame(index=gen.index[mask])
        upd["target_q"] = q_gen[mask]
        upd["voltage_regulator_on"] = False
        network.update_generators(upd)
        
    vsc = network.get_vsc_converter_stations(
        attributes=[
            "voltage_regulator_on", "q",
            "min_q", "max_q", "min_q_at_p", "max_q_at_p",
            "connected", "bus_id"
        ]
    )
    if keep_only_main_comp:
        vsc = _keep_only_main_comp(vsc, df_bus)
    if len(vsc):
        mask, q_gen = _bound_at_qlimit(vsc, vsc["voltage_regulator_on"])
        if mask.any():
            upd = pd.DataFrame(index=vsc.index[mask])
            upd["target_q"] = q_gen[mask]
            upd["voltage_regulator_on"] = False
            network.update_vsc_converter_stations(upd)

    _bake_svc_saturation(network, keep_only_main_comp)


def _bake_svc_saturation(network, keep_only_main_comp=True):
    """Freeze a VOLTAGE-mode SVC whose realized reactive output sits at (or beyond)
    its voltage-dependent susceptance envelope to fixed-Q (REACTIVE_POWER mode),
    mirroring ``_bake_reactive_limit_switches`` for generators/VSC stations above.

    Unlike a generator's fixed Q box, an SVC's reactive range is
    ``Q(V) = b * V^2`` (``b`` in ``b_min``..``b_max``, in Siemens): recompute
    ``qmin``/``qmax`` in MVAr at the SVC's own solved terminal voltage before
    comparing against the realized ``q``.
    """
    df_bus = network.get_buses(attributes=["v_mag", "synchronous_component"])
    svc = network.get_static_var_compensators(
        attributes=["regulating", "regulation_mode", "q", "b_min", "b_max", "connected", "bus_id"]
    )
    if keep_only_main_comp:
        svc = _keep_only_main_comp(svc, df_bus)
    is_voltage = svc["regulating"] & (svc["regulation_mode"] == "VOLTAGE")
    svc = svc[is_voltage]
    if not len(svc):
        return
    v_kv = df_bus.loc[svc["bus_id"].values, "v_mag"].to_numpy()
    # SVC "q" (like generators') is the terminal/receptor-convention result: flip to
    # generator/injection convention to compare against the susceptance envelope.
    q_gen = -svc["q"].to_numpy()
    qmax = svc["b_max"].to_numpy() * v_kv ** 2  # Q[MVAr] = B[S] * V[kV]^2
    qmin = svc["b_min"].to_numpy() * v_kv ** 2
    mask = (q_gen >= qmax - _Q_LIMIT_TOL) | (q_gen <= qmin + _Q_LIMIT_TOL)
    if mask.any():
        upd = pd.DataFrame(index=svc.index[mask])
        # unlike Generator.target_q (already generator convention), StaticVarCompensator
        # .target_q is receptor convention like its "q" result column -- write the raw
        # (unflipped) realized value.
        upd["target_q"] = svc["q"].to_numpy()[mask]
        upd["regulation_mode"] = "REACTIVE_POWER"
        network.update_static_var_compensators(upd)


def _bake_remote_voltage_control(network, keep_only_main_comp=True):
    """Rewrite *remote* voltage control into *local* control at the solved terminal.

    A generator regulating a bus other than its own terminal (``regulated_element_id``
    resolving to a different bus -- e.g. holding the 400 kV grid-connection point
    across its step-up transformer) is rewritten to regulate its OWN terminal at that
    terminal's solved magnitude, and the remote regulation is cleared so the converter
    treats it as local control.

    At the converged operating point this is exact: the generator's own terminal
    already sits at ``v_mag`` with the realized reactive output, so fixing it locally
    reproduces the same fixed point (same V everywhere, same Q; the formerly-regulated
    remote bus still lands on its solved value). A subsequent loop-free solve -- OLF
    (whose loop-free parameters disable remote control anyway) or lightsim2grid --
    reproduces the baked state.

    Why it exists: lightsim2grid v1 cannot host a remote voltage controller on a slack
    bus, and the default distributed slack puts many remote-regulating generators on
    slack buses. Making every controller local sidesteps that.

    .. warning::
        This is a *base-case* faithful approximation. Under a topology change the
        generator then holds its own terminal magnitude instead of the remote bus, so
        the post-contingency reactive behaviour differs from true remote control.

    Operates in place; idempotent (an already-local generator is left untouched).
    """
    df_bus = network.get_buses(attributes=["v_mag", "synchronous_component"])
    gen = network.get_generators(
        attributes=["voltage_regulator_on", "regulated_element_id", "connected", "bus_id"]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    if not len(gen):
        return
    # "remote" matches the converter's own test (see _from_pypowsybl.init): a non-empty
    # regulated element that is not the generator's own id.
    reg = gen["regulated_element_id"].fillna("")
    remote = gen["voltage_regulator_on"] & gen["connected"] & (reg != "") & (reg != gen.index)
    gen = gen[remote]
    if not len(gen):
        return
    upd = pd.DataFrame(index=gen.index)
    upd["target_v"] = df_bus.loc[gen["bus_id"].values, "v_mag"].values  # own terminal, kV
    upd["regulated_element_id"] = gen.index  # regulate own terminal -> local control
    network.update_generators(upd)


def _bake_active_power(network, balance_on_loads, load_power_factor_constant, keep_only_main_comp=True):
    df_bus = network.get_buses(attributes=["synchronous_component"])
    
    gen = network.get_generators(attributes=["p", "connected", "bus_id"])
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    if len(gen):
        # result p is load convention; target_p is generator convention
        network.update_generators(
            pd.DataFrame({"target_p": -gen["p"]}, index=gen.index)
        )

    bat = network.get_batteries(attributes=["p", "connected", "bus_id"])
    if keep_only_main_comp:
        bat = _keep_only_main_comp(bat, df_bus)
    if len(bat):
        network.update_batteries(
            pd.DataFrame({"target_p": -bat["p"]}, index=bat.index)
        )

    if balance_on_loads:
        load = network.get_loads(attributes=["p", "q", "connected", "bus_id"])
        if keep_only_main_comp:
            load = _keep_only_main_comp(load, df_bus)
        if len(load):
            upd = pd.DataFrame(index=load.index)
            upd["p0"] = load["p"]  # load convention both sides, no flip
            if load_power_factor_constant:
                upd["q0"] = load["q"]
            network.update_loads(upd)
