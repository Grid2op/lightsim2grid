# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import copy

import numpy as np
import pandas as pd

from ._aux_common import _aux_get_bus, _aux_regulated_bus_view_ids


def _aux_svc_can_be_pv_flags(can_be_pv, svc_index):
    """The static var compensators the ``can_be_pv`` argument of `init` flags, as a
    boolean array aligned on ``svc_index`` (the SVCs in lightsim2grid order): the ids it
    holds (what `bake_outer_loops` returns for the SVCs it froze), or the
    True entries of a boolean Series indexed by id. A boolean array is in the generators'
    order, so it flags no SVC. ``None`` when nothing is flagged."""
    if can_be_pv is None or len(svc_index) == 0:
        return None
    if isinstance(can_be_pv, pd.Series):
        flags = can_be_pv.reindex(svc_index).fillna(False).astype(bool).to_numpy()
    elif isinstance(can_be_pv, np.ndarray) and can_be_pv.dtype == bool:
        return None
    else:
        flags = svc_index.isin(pd.Index([str(el) for el in can_be_pv]))
    flags = np.asarray(flags, dtype=bool)
    return flags if flags.any() else None


def _svc_standby_b0(net, svc_index):
    """The ``b0`` (S) of the ``standbyAutomaton`` extension of each SVC of ``svc_index``, 0
    for an SVC without it, as a ``numpy`` array.

    OpenLoadFlow models an SVC carrying that extension as its ``b0``, a fixed susceptance,
    PLUS the SVC itself -- in standby or not -- and holds the SVC's own susceptance in
    ``[b_min, b_max]``: its ``ReactiveLimits`` loop compares what the SVC part produces with
    those limits. The SVC's total output, the one pypowsybl reports and the one
    lightsim2grid models, can therefore only range over ``[b_min + b0, b_max + b0]``."""
    b0 = np.zeros(len(svc_index))
    try:
        automaton = net.get_extensions("standbyAutomaton")
    except Exception:  # noqa: BLE001 - extension unsupported / absent on old pypowsybl
        return b0
    if automaton is None or not automaton.shape[0] or "b0" not in automaton.columns:
        return b0
    val = automaton["b0"].reindex(pd.Index(svc_index)).to_numpy(float)
    return np.where(np.isfinite(val), val, 0.)


def _aux_add_svc(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl, sn_mva_used,
                 can_be_pv=None, olf_rules=None, olf_vc=None):
    """Add every Static Var Compensator (SVC) of ``net`` to ``model``: VOLTAGE
    (local/remote, optional slope), REACTIVE_POWER (fixed Q) or OFF, all solved
    through the bordered VoltageControl NR extension. A grid with no SVC declares
    no controller and stays byte-identical to before this feature. Returns
    ``df_svc`` (its per-substation ids are not needed downstream, unlike every
    other element type here).

    The SVCs ``can_be_pv`` flags (see `_aux_svc_can_be_pv_flags`) are the ones an
    outer loop froze out of voltage control (what `bake_outer_loops` returns). One
    whose ``standbyAutomaton`` extension says ``standby`` is a standby SVC left
    idle: its thresholds go to ``LSGrid.set_svc_standby`` (the check of its switch
    on). Any other is an SVC frozen at a reactive limit: ``LSGrid.set_svc_can_be_pv``
    (the check of its release, as for a generator).

    With ``olf_rules`` (and ``olf_vc``, :func:`._olf_rules.voltage_controllers`), OpenLoadFlow's
    loading: an SVC whose voltage control it discards is off, one it makes a voltage monitor
    is off too and flagged standby (``LSGrid.set_svc_standby``, the VoltageMonitoring loop
    switches it on), the slope is ignored unless ``voltage_per_reactive_power_control``,
    and without ``svc_voltage_monitoring`` the automaton's ``b0`` does not exist."""
    if sort_index:
        df_svc = net.get_static_var_compensators().sort_index()
    else:
        df_svc = net.get_static_var_compensators()
    nb_svc = df_svc.shape[0]
    svc_bus, svc_disco, svc_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "svc", df_svc)

    # SvcContainer.RegulationMode: OFF=0, VOLTAGE=1, REACTIVE_POWER=2
    OFF_MODE, VOLTAGE_MODE, REACTIVE_POWER_MODE = 0, 1, 2
    svc_mode = np.zeros(nb_svc, dtype=int)
    svc_reg_bus = svc_bus.copy()             # regulated bus (own bus unless remote)
    svc_reg_vn = np.ones(nb_svc)             # nominal v (kV) of the regulated bus
    svc_slope_pu = np.zeros(nb_svc)
    svc_target_q_inject = np.zeros(nb_svc)
    if nb_svc:
        mode_str = df_svc["regulation_mode"].values.astype(str)
        if "regulating" in df_svc.columns:
            regulating = df_svc["regulating"].values.astype(bool)
        else:
            # legacy pypowsybl: "OFF" is encoded directly in the regulation mode
            regulating = np.ones(nb_svc, dtype=bool)
        svc_mode[(mode_str == "VOLTAGE") & regulating] = VOLTAGE_MODE
        svc_mode[(mode_str == "REACTIVE_POWER") & regulating] = REACTIVE_POWER_MODE

        # resolve the regulated bus, mirroring the generator `regulated_element_id`
        # logic above (busbar section in node/breaker grids, or any bus-connected
        # element). Local SVCs keep their own bus.
        svc_vl = copy.deepcopy(df_svc["voltage_level_id"].values)
        if "regulated_element_id" in df_svc.columns:
            reg_id = copy.deepcopy(df_svc["regulated_element_id"].values)
            reg_id = np.where(reg_id == "", df_svc.index, reg_id)
            mask_svc_remote = reg_id != df_svc.index.values
            if mask_svc_remote.any():
                # TODO: resolved once at import; if the regulated element later changes
                # bus inside lightsim2grid this stays frozen and desynchronises from the
                # original grid (see `_aux_regulated_bus_view_ids` warning).
                remote_svc_idx = np.nonzero(mask_svc_remote)[0]
                svc_reg_bus_view = _aux_regulated_bus_view_ids(net, reg_id[mask_svc_remote])
                # same disconnected-remote-target situation as for generators above:
                # fall back to local control rather than crashing on an unresolvable
                # (disconnected) regulated element.
                unresolved_svc = svc_reg_bus_view == ""
                if unresolved_svc.any():
                    mask_svc_remote[remote_svc_idx[unresolved_svc]] = False
                    svc_reg_bus_view = svc_reg_bus_view[~unresolved_svc]
                svc_reg_bus[mask_svc_remote] = bus_df.loc[svc_reg_bus_view, "bus_global_id"].values
                svc_vl[mask_svc_remote] = bus_df.loc[svc_reg_bus_view, "voltage_level_id"].values
        svc_reg_vn = voltage_levels.loc[svc_vl, "nominal_v"].values

        # IIDM gives the SVC reactive setpoint in the receptor (load) convention,
        # whereas lightsim2grid stamps Q with the generator-injection convention
        # (Phase 0 probe: SVC target_q=+30 absorbs 30 MVar) -> negate.
        target_q = df_svc["target_q"].values.astype(float)
        target_q = np.where(np.isfinite(target_q), target_q, 0.)
        svc_target_q_inject = -target_q

        # optional voltage/reactive-power slope ("droop"), in kV/MVar:
        #   s_pu = slope[kV/MVar] * sn_mva / vn_kv(regulated bus)   (Phase 0 probe #1)
        # Read from `net` (not `net_pu`): pypowsybl's per-unit view (native
        # `per_unit=True` and the legacy `PerUnitView`) does not per-unit this
        # extension at all -- `slope` comes back numerically identical whether
        # or not `per_unit` is set (checked empirically) -- so there is no
        # `net_pu` value to defer to here; the conversion has to be done by hand.
        try:
            df_slope = net.get_extensions("voltagePerReactivePowerControl")
        except Exception:
            # extension tables may be unavailable on (very) old pypowsybl versions
            df_slope = None
        if df_slope is not None and df_slope.shape[0]:
            for svc_pos, svc_id in enumerate(df_svc.index):
                if svc_id in df_slope.index:
                    slope_kv_per_mvar = float(df_slope.loc[svc_id, "slope"])
                    svc_slope_pu[svc_pos] = slope_kv_per_mvar * sn_mva_used / svc_reg_vn[svc_pos]

    olf_monitor = np.zeros(nb_svc, dtype=bool)
    olf_b0 = None
    if olf_rules is not None and nb_svc:
        if not olf_rules.voltage_per_reactive_power_control:
            svc_slope_pu[:] = 0.
        if olf_vc is not None:
            vc = olf_vc.reindex(df_svc.index)
            olf_monitor = vc["monitor"].fillna(False).to_numpy(bool)
            # discarded or idle until its automaton switches it on: no reactive output
            svc_mode[vc["discarded"].fillna(False).to_numpy(bool) | olf_monitor] = OFF_MODE

    if nb_svc:
        target_v = df_svc["target_v"].values.astype(float)
        # target_v (kV) -> pu at the regulated bus; NaN (REACTIVE_POWER / OFF) -> 1 pu
        target_vm_pu = np.where(np.isfinite(target_v), target_v, svc_reg_vn) / svc_reg_vn
        # pypowsybl/IIDM gives b_min/b_max in SIEMENS (physical susceptance), while
        # lightsim2grid's SvcContainer expects them in per unit (base sn_mva, at the
        # regulated bus's nominal voltage -- same base as target_vm_pu/slope above).
        # Without this conversion the SVC's modeled Q range is smaller than its real one
        # by a factor of (nominal_v_kv)^2 / sn_mva (eg ~500x on a 225kV/100MVA bus),
        # making a perfectly healthy SVC collapse to a near-zero Q range: it saturates
        # (or, since check_solution never enforces SVC Q limits, silently hits a "hard"
        # voltage pin at its own target_vm_pu regardless of what Q that would truly take)
        # long before the real device would.
        # This CANNOT be replaced by reading `net_pu.get_static_var_compensators()`
        # instead: checked empirically that pypowsybl's per-unit view (native
        # `per_unit=True` and the legacy `PerUnitView`) leaves `b_min`/`b_max`
        # in SIEMENS even under `per_unit=True` -- it only per-units `target_v`
        # (and, generically, elements it fully models) -- so relying on it here
        # would silently reintroduce this exact bug.
        # An SVC carrying a standby automaton: OpenLoadFlow holds the SVC part apart from
        # the automaton's fixed b0, so the total output lightsim2grid models ranges over the
        # shifted interval (see `_svc_standby_b0`).
        b0 = _svc_standby_b0(net, df_svc.index)
        if olf_rules is not None:
            if olf_rules.svc_voltage_monitoring:
                # OpenLoadFlow's: a shunt of its own (LSGrid.set_svc_b0), carried by an SVC
                # that regulates voltage in the network, in standby or not
                regulates_voltage = (mode_str == "VOLTAGE") & regulating
                olf_b0 = np.where(regulates_voltage, b0, 0.)
            # the range is the SVC's own; without monitoring the automaton is not read
            b0 = np.zeros(nb_svc)
        b_min = (df_svc["b_min"].values.astype(float) + b0) * (svc_reg_vn ** 2) / sn_mva_used
        b_max = (df_svc["b_max"].values.astype(float) + b0) * (svc_reg_vn ** 2) / sn_mva_used
    else:
        target_vm_pu = np.zeros(0)
        b_min = np.zeros(0)
        b_max = np.zeros(0)

    model.init_svcs([int(m) for m in svc_mode],
                    target_vm_pu,
                    svc_target_q_inject,
                    svc_slope_pu,
                    b_min,
                    b_max,
                    svc_reg_bus.astype(np.int32),
                    svc_bus.astype(np.int32))
    for svc_id, disco in enumerate(svc_disco):
        if disco:
            model.deactivate_svc(svc_id)
    model.set_svc_names(df_svc.index)
    own_vn = voltage_levels.loc[df_svc["voltage_level_id"].values, "nominal_v"].to_numpy(float) if nb_svc else np.zeros(0)
    if olf_b0 is not None and (olf_b0 != 0.).any():
        # S -> pu at the SVC's own bus
        model.set_svc_b0(olf_b0 * own_vn ** 2 / sn_mva_used)

    # the SVCs an outer loop froze out of voltage control (what `bake_outer_loops` returns):
    # nothing in the powerflow reads the flags, they only open them to a physical check
    flagged = _aux_svc_can_be_pv_flags(can_be_pv, df_svc.index)
    if flagged is None:
        flagged = np.zeros(nb_svc, dtype=bool)
    try:
        automaton = net.get_extensions("standbyAutomaton")
    except Exception:  # noqa: BLE001 - extension unsupported / absent on old pypowsybl
        automaton = None
    if automaton is not None and automaton.shape[0]:
        in_standby = automaton.index[automaton["standby"].astype(bool)]
    else:
        in_standby = pd.Index([], dtype=object)
    # an idle standby SVC -- one the bake left idle, or one OpenLoadFlow loads as a voltage
    # monitor: the switch on of its automaton. The thresholds are compared with the voltage
    # of the bus it regulates, as OLF does, but in pu of the nominal voltage of the SVC's own
    # voltage level (LfStaticVarCompensatorImpl: the same one for a local SVC).
    standby = (flagged & df_svc.index.isin(in_standby)) | olf_monitor
    if standby.any():
        ids = df_svc.index[standby]
        low_vm_pu = np.full(nb_svc, np.nan)
        high_vm_pu = np.full(nb_svc, np.nan)
        low_vm_pu[standby] = automaton.loc[ids, "low_voltage_threshold"].to_numpy(float) / own_vn[standby]
        high_vm_pu[standby] = automaton.loc[ids, "high_voltage_threshold"].to_numpy(float) / own_vn[standby]
        # a monitor's set-points, for the VoltageMonitoring loop (the bake's idle SVCs have none)
        low_target = np.full(nb_svc, np.nan)
        high_target = np.full(nb_svc, np.nan)
        if olf_monitor.any():
            mon_ids = df_svc.index[olf_monitor]
            low_target[olf_monitor] = automaton.loc[mon_ids, "low_voltage_setpoint"].to_numpy(float) / own_vn[olf_monitor]
            high_target[olf_monitor] = automaton.loc[mon_ids, "high_voltage_setpoint"].to_numpy(float) / own_vn[olf_monitor]
        model.set_svc_standby(standby, low_vm_pu, high_vm_pu, low_target, high_target)
    # any other flagged one: frozen at a reactive limit, the release of the generators' can_be_pv
    release = flagged & ~standby
    if release.any():
        model.set_svc_can_be_pv(release)

    return df_svc
