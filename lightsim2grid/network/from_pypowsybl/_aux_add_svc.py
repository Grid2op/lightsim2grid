# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import copy

import numpy as np

from ._aux_common import _aux_get_bus, _aux_regulated_bus_view_ids


def _aux_add_svc(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl, sn_mva_used):
    """Add every Static Var Compensator (SVC) of ``net`` to ``model``: VOLTAGE
    (local/remote, optional slope), REACTIVE_POWER (fixed Q) or OFF, all solved
    through the bordered VoltageControl NR extension. A grid with no SVC declares
    no controller and stays byte-identical to before this feature. Returns
    ``df_svc`` (its per-substation ids are not needed downstream, unlike every
    other element type here)."""
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
        b_min = df_svc["b_min"].values.astype(float) * (svc_reg_vn ** 2) / sn_mva_used
        b_max = df_svc["b_max"].values.astype(float) * (svc_reg_vn ** 2) / sn_mva_used
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

    return df_svc
