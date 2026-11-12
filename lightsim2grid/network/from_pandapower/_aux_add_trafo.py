# Copyright (c) 2020-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pandapower as pp
from packaging import version

from ._pp_bus_to_ls_bus import pp_bus_to_ls
from ._my_const import _MIN_PP_VERSION_ADV_GRID_MODEL, ALLOWED_PP_ORIG_FILE


def _aux_add_trafo(
    converter,
    model,
    pp_net,
    pp_to_ls,
    pp_orig_file : ALLOWED_PP_ORIG_FILE = "pandapower_v2"
    ):
    """
    Add the transformers of the pp_net into the lightsim2grid "model"

    Parameters
    ----------
    converter
    model
    pp_net

    Returns
    -------

    """
    shift_ = pp_net.trafo["shift_degree"].to_numpy().copy()
    if np.any(~np.isfinite(shift_)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"shift_degree\"], they have been replaced by 0")
    shift_[~np.isfinite(shift_)] = 0.

    # the tap changers, applied one after the other as pandapower's `_calc_tap_from_dataframe`
    # does: the first ("tap_*" columns), and the second ("tap2_*") when there is one
    tap1 = _read_tap_changer(pp_net, "tap")
    tap2 = _read_tap_changer(pp_net, "tap2") if "tap2_pos" in pp_net.trafo else None
    tap_pos = tap1.pos
    # a transformer without a tap side has no tap: it is stored as not on the hv side, as always
    is_tap_hv_side = tap1.on_hv
    is_tabular = tap1.is_tabular
    tap_table = _read_tap_table(pp_net, is_tabular)

    # not attached to the 2-winding checks: those always find their column
    if "trafo3w" in pp_net and "tap_changer_type" in pp_net.trafo3w:
        if np.any(_str_eq(pp_net.trafo3w["tap_changer_type"], "Ideal")):
            raise RuntimeError("Ideal phase shifters are not modeled. Please remove all 3-winding trafos "
                               "with \"tap_changer_type\" set to \"Ideal\".")

    trafo_model_is_t = True
    if "_options" in pp_net and "trafo_model" in pp_net._options:
        trafo_model_is_t = pp_net._options["trafo_model"] == "t"
    
    # compute physical parameters
    use_pp3 = version.parse(pp.__version__) >= _MIN_PP_VERSION_ADV_GRID_MODEL and pp_orig_file == "pandapower_v3"
    if not use_pp3 and not trafo_model_is_t:
        raise RuntimeError("Cannot convert a transformer with model 'pi' to LightSim2grid (using pandapower < 3)")

    # the share of the leakage impedance on the hv side, in the T model (0.5: half and half)
    leak_r = _ratio_column(pp_net, "leakage_resistance_ratio_hv")
    leak_x = _ratio_column(pp_net, "leakage_reactance_ratio_hv")
    leak_asym = trafo_model_is_t & ((np.abs(leak_r - 0.5) > 1e-7) | (np.abs(leak_x - 0.5) > 1e-7))
    if np.any(leak_asym) and not use_pp3:
        warnings.warn("leakage_resistance_ratio_hv / leakage_reactance_ratio_hv != 0.5 are only supported with "
                      "the pandapower 3 converter. They will be replaced by 0.5")
        leak_asym[:] = False
    vn_hv_bus = pp_net.bus.loc[pp_net.trafo["hv_bus"]]["vn_kv"].to_numpy()
    vn_lv_bus = pp_net.bus.loc[pp_net.trafo["lv_bus"]]["vn_kv"].to_numpy()
    # the voltages the tap moves: the rated ones, which may differ from the buses' (off-nominal ratio)
    vn_hv_ref = pp_net.trafo["vn_hv_kv"].to_numpy(dtype=float)
    vn_lv_ref = pp_net.trafo["vn_lv_kv"].to_numpy(dtype=float)
    # `parallel` identical transformers in one row: impedance / parallel, admittance * parallel
    parallel = np.ones(pp_net.trafo.shape[0])
    if "parallel" in pp_net.trafo:
        parallel = pp_net.trafo["parallel"].to_numpy(dtype=float)
    per_trafo = (pp_net.trafo["vk_percent"].to_numpy(),
                 pp_net.trafo["vkr_percent"].to_numpy(),
                 pp_net.trafo["sn_mva"].to_numpy(),
                 pp_net.trafo["pfe_kw"].to_numpy(),
                 pp_net.trafo["i0_percent"].to_numpy())

    def tap_model(trafo_ids, positions):
        """The voltages (kV) of both windings and the phase shift (degree) of the transformers
        `trafo_ids` with their tap at `positions`, as pandapower's `_calc_tap_from_dataframe`."""
        vnh = vn_hv_ref[trafo_ids].astype(float)
        vnl = vn_lv_ref[trafo_ids].astype(float)
        shift = shift_[trafo_ids].astype(float)
        _apply_tap_changer(tap1, trafo_ids, positions, vnh, vnl, shift, tap_table)
        if tap2 is not None:
            _apply_tap_changer(tap2, trafo_ids, tap2.pos[trafo_ids], vnh, vnl, shift, None)
        return vnh, vnl, shift

    def trafo_params(trafo_ids, positions):
        """r, x, the shunt admittance of side 1 and side 2 (pu), the ratio and the phase shift (degree)
        of the transformers `trafo_ids` with their
        taps at `positions` -- the one place they are computed, so that a tap position the
        step tables below describe is exactly the transformer this converter would build at
        that position.

        The ratio is the one `LSGrid.init_trafo` expects on the tap side; the impedances are pandapower's
        `_calc_r_x_y_from_dataframe`: the impedance at the low voltage bus base, scaled by
        `(vn_trafo_lv / vn_lv) ** 2 / parallel` (which commutes with the T to pi conversion)."""
        trafo_ids = np.asarray(trafo_ids, dtype=int)
        vnh, vnl, shift = tap_model(trafo_ids, positions)
        vk, vkr, sn, pfe, i0 = [np.asarray(el, dtype=float)[trafo_ids] for el in per_trafo]
        tabular = is_tabular[trafo_ids]
        if np.any(tabular):
            # pandapower reads them in the table whatever the tap side (`_get_vk_values_from_table`)
            vk_t, vkr_t = tap_table.at(trafo_ids, positions, ("vk_percent", "vkr_percent"), (vk, vkr))
            vk = np.where(tabular, vk_t, vk)
            vkr = np.where(tabular, vkr_t, vkr)
        vn_hv, vn_lv = vn_hv_bus[trafo_ids], vn_lv_bus[trafo_ids]
        no_tap = np.zeros(trafo_ids.shape[0])
        hv_side = [True] * trafo_ids.shape[0]
        if use_pp3:
            # use pandapower version 3 converter in this case.
            # We use it because:
            # - the grid comes from pandapower3 (eg pn.case118() with pandapower 3 installed)
            # - AND the pandapower version is >= 3
            r, x, b = converter.get_trafo_param_pp3(no_tap, no_tap, no_tap, hv_side, vn_hv, vn_lv,
                                                    vk, vkr, sn, pfe, i0, trafo_model_is_t)
        else:
            # default legacy mode: use the pandapower 2 converter
            r, x, b = converter.get_trafo_param_pp2(no_tap, no_tap, no_tap, hv_side, vn_hv, vn_lv,
                                                    vk, vkr, sn, pfe, i0)
        r, x = np.array(r), np.array(x)
        h_1 = 0.5 * np.asarray(b)
        h_2 = h_1.copy()
        asym = leak_asym[trafo_ids]
        if np.any(asym):
            # the T model with its own split of the leakage impedance: from the pi values, as
            # pandapower's `_wye_delta`, which gives a different shunt admittance on each side
            r_pi, x_pi, y_pi = converter.get_trafo_param_pp3(no_tap[asym], no_tap[asym], no_tap[asym],
                                                             hv_side[:int(asym.sum())], vn_hv[asym], vn_lv[asym],
                                                             vk[asym], vkr[asym], sn[asym], pfe[asym], i0[asym],
                                                             False)
            r[asym], x[asym], h_1[asym], h_2[asym] = _wye_delta(np.asarray(r_pi), np.asarray(x_pi),
                                                                np.asarray(y_pi),
                                                                leak_r[trafo_ids][asym], leak_x[trafo_ids][asym])
        z_factor = (vnl / vn_lv) ** 2 / parallel[trafo_ids]
        ratio = (vnh / vnl) / (vn_hv / vn_lv)  # pandapower's, on the high voltage side
        ratio = np.where(is_tap_hv_side[trafo_ids], ratio, 1. / ratio)
        return r * z_factor, x * z_factor, h_1 / z_factor, h_2 / z_factor, ratio, shift

    all_trafos = np.arange(pp_net.trafo.shape[0])
    trafo_r, trafo_x, trafo_h_1, trafo_h_2, trafo_ratio, trafo_shift = trafo_params(all_trafos, tap_pos)
    trafo_b = trafo_h_1 + trafo_h_2

    # initialize the grid
    do_ignore_tap_side_for_phase_shift = True
    model.init_trafo(trafo_r,
                     trafo_x,
                     trafo_b,
                     trafo_ratio,
                     trafo_shift,
                     is_tap_hv_side,
                     pp_bus_to_ls(pp_net.trafo["hv_bus"].values, pp_to_ls),
                     pp_bus_to_ls(pp_net.trafo["lv_bus"].values, pp_to_ls),
                     do_ignore_tap_side_for_phase_shift)
    if np.any(leak_asym):
        model.update_trafos_parameters_full(trafo_r, trafo_x, trafo_h_1, trafo_h_2)

    for tr_id, is_connected in enumerate(pp_net.trafo["in_service"].values):
        if not is_connected:
            # trafo is deactivated
            model.deactivate_trafo(tr_id)

    # the tap changers (what the TransformerVoltageControl outer loop moves), then what they regulate
    with_steps = _aux_trafo_tap_changers(model, pp_net, tap1, tap_table, trafo_params, trafo_r, trafo_x, trafo_b)
    _aux_trafo_tap_regulation(model, pp_net, pp_to_ls, with_steps)

    
def _str_eq(series, value):
    """`series == value` as a plain numpy boolean array, a missing value counting as different.
    Needed since pandapower 4, whose text columns ("tap_side", "tap_changer_type", ...) are of the
    pandas "string" dtype: there a missing value is pd.NA, and comparing to it gives pd.NA rather
    than False ("boolean value of NA is ambiguous" in any later numpy test). copy=True because with
    an object column, pandas 3 can hand back a read-only array."""
    return series.eq(value).to_numpy(dtype=bool, na_value=False, copy=True)
  
  
def _read_tap_changer(pp_net, prefix):
    """One tap changer of the transformers, from the `prefix`_* columns ("tap" or "tap2"), with the
    missing values filled and the kind of tap pandapower applies (`_calc_tap_from_dataframe`):
    a tap needs a side; then, with pandapower 3 data, a "Ratio" or "Symmetrical" tap changer is
    `complex`, an "Ideal" one `ideal` (with pandapower 2 data: tap_phase_shifter), one with a
    tap_dependency_table `tabular`. The others are at ratio 1 whatever their position."""
    trafo = pp_net.trafo
    nb = trafo.shape[0]

    def column(name, fill, what):
        col = f"{prefix}_{name}"
        vals = trafo[col].to_numpy(dtype=float, na_value=np.nan, copy=True) if col in trafo else np.full(nb, np.nan)
        if col in trafo and np.any(~np.isfinite(vals)):
            warnings.warn(f"There were some Nan in the pp_net.trafo[\"{col}\"], they have been replaced by {what}")
        bad = ~np.isfinite(vals)
        vals[bad] = fill[bad] if isinstance(fill, np.ndarray) else fill
        return vals

    tc = SimpleNamespace()
    tc.neutral = column("neutral", 0., "0")
    tc.step_pct = column("step_percent", 0., "0")
    tc.pos = column("pos", tc.neutral, f"the neutral position (pp_net.trafo[\"{prefix}_neutral\"])")
    tc.step_degree = column("step_degree", 0., "0")  # in degree
    if f"{prefix}_side" in trafo:
        side = trafo[f"{prefix}_side"]
        tc.on_hv = _str_eq(side, "hv")
        applies = tc.on_hv | _str_eq(side, "lv")
    else:
        tc.on_hv = np.zeros(nb, dtype=bool)
        applies = np.zeros(nb, dtype=bool)
    has_step = (tc.step_pct != 0.) & (tc.pos != tc.neutral)
    if np.any(~applies & has_step):
        warnings.warn(f"There were some Nan in the pp_net.trafo[\"{prefix}_side\"]: the tap of these "
                      "transformers is ignored, as pandapower does.")

    tc.is_tabular = np.zeros(nb, dtype=bool)
    if f"{prefix}_dependency_table" in trafo:
        tc.is_tabular = trafo[f"{prefix}_dependency_table"].fillna(False).to_numpy(dtype=bool)
        if prefix != "tap" and np.any(tc.is_tabular):
            raise RuntimeError(f"Tabular tap changers (pp_net.trafo[\"{prefix}_dependency_table\"]) are only "
                               "handled for the first tap changer.")
    if f"{prefix}_changer_type" in trafo:
        changer_type = trafo[f"{prefix}_changer_type"]
        is_ideal = _str_eq(changer_type, "Ideal")
        is_complex = _str_eq(changer_type, "Ratio") | _str_eq(changer_type, "Symmetrical")
        if np.any(~is_complex & ~is_ideal & ~tc.is_tabular & applies & has_step):
            warnings.warn(f"Some transformers have a tap but no pp_net.trafo[\"{prefix}_changer_type\"]: their "
                          "tap is ignored, as pandapower does.")
    else:
        is_ideal = np.zeros(nb, dtype=bool)
        if f"{prefix}_phase_shifter" in trafo:
            is_ideal = trafo[f"{prefix}_phase_shifter"].fillna(False).to_numpy(dtype=bool)
        is_complex = ~is_ideal
    tc.tabular = applies & tc.is_tabular
    tc.ideal = applies & is_ideal & ~tc.is_tabular
    tc.complex = applies & is_complex & ~tc.is_tabular
    if np.any(tc.ideal & (tc.step_pct != 0.) & (tc.step_degree != 0.)):
        # as pandapower
        raise RuntimeError(f"Both {prefix}_step_degree and {prefix}_step_percent set for an ideal phase shifter "
                           f"(pp_net.trafo[\"{prefix}_changer_type\"] set to \"Ideal\").")
    return tc


def _apply_tap_changer(tc, trafo_ids, positions, vnh, vnl, shift, tap_table):
    """Move, in place, the winding voltages `vnh`, `vnl` (kV) and the `shift` (degree) of the
    transformers `trafo_ids` by their tap changer `tc` (`_read_tap_changer`) at `positions`."""
    positions = np.asarray(positions, dtype=float)
    diff = positions - tc.neutral[trafo_ids]
    step = tc.step_pct[trafo_ids]
    angle = tc.step_degree[trafo_ids]
    on_hv = tc.on_hv[trafo_ids]
    applies = tc.complex[trafo_ids]
    ideal = tc.ideal[trafo_ids]
    tabular = tc.tabular[trafo_ids]
    for mask, vn, direction in ((applies & on_hv, vnh, 1.), (applies & ~on_hv, vnl, -1.)):
        vn[mask], tap_shift = _complex_tap(vn[mask], step[mask], diff[mask], angle[mask])
        shift[mask] += direction * tap_shift
    for mask, direction in ((ideal & on_hv, 1.), (ideal & ~on_hv, -1.)):
        shift[mask] += direction * _ideal_tap(step[mask], diff[mask], angle[mask])
    if np.any(tabular):
        ratio_t, angle_t = tap_table.at(trafo_ids, positions, ("voltage_ratio", "angle_deg"), (1., 0.))
        for mask, vn, direction in ((tabular & on_hv, vnh, 1.), (tabular & ~on_hv, vnl, -1.)):
            vn[mask] *= ratio_t[mask]
            shift[mask] += direction * angle_t[mask]


class _TapTable:
    """The rows of ``net.trafo_characteristic_table`` of the tabular transformers, by
    (transformer, step)."""
    def __init__(self, id_characteristic, table):
        self.id_characteristic = id_characteristic
        self.rows = {}
        for row in table.to_dict("records"):
            self.rows[(row["id_characteristic"], int(row["step"]))] = row

    def has(self, trafo_id, position):
        return (self.id_characteristic[trafo_id], int(position)) in self.rows

    def at(self, trafo_ids, positions, columns, defaults):
        """The `columns` of the transformers `trafo_ids` at `positions` (one array per column);
        `defaults` (a value or an array like `trafo_ids`) where the table has no such column or
        the transformer no table."""
        positions = np.broadcast_to(np.asarray(positions, dtype=float), np.shape(trafo_ids))
        res = []
        for col, default in zip(columns, defaults):
            out = np.array(np.broadcast_to(np.asarray(default, dtype=float), np.shape(trafo_ids)))
            for k, (tid, pos) in enumerate(zip(trafo_ids, positions)):
                row = self.rows.get((self.id_characteristic[tid], int(pos)))
                if row is not None and col in row and np.isfinite(row[col]):
                    out[k] = row[col]
            res.append(out)
        return res


def _read_tap_table(pp_net, is_tabular):
    """The tap table of the tabular transformers (`_TapTable`), checked as pandapower does: each
    has an id_characteristic_table, and its row at the current tap_pos exists."""
    if not np.any(is_tabular):
        return None
    table = pp_net.get("trafo_characteristic_table", None)
    if table is None or not hasattr(table, "shape") or table.shape[0] == 0:
        raise RuntimeError("Some pp_net.trafo[\"tap_dependency_table\"] are True but there is no "
                           "pp_net.trafo_characteristic_table.")
    id_char = pp_net.trafo["id_characteristic_table"].to_numpy() \
        if "id_characteristic_table" in pp_net.trafo else np.full(is_tabular.shape[0], np.nan)
    if np.any(is_tabular & pd.isna(id_char)):
        raise RuntimeError("Some pp_net.trafo[\"tap_dependency_table\"] are True with no "
                           "pp_net.trafo[\"id_characteristic_table\"].")
    tap_table = _TapTable(id_char, table)
    tap_pos = pp_net.trafo["tap_pos"].to_numpy(dtype=float)
    missing = [pp_net.trafo.index[tid] for tid in np.flatnonzero(is_tabular)
               if not np.isfinite(tap_pos[tid]) or not tap_table.has(tid, tap_pos[tid])]
    if missing:
        raise RuntimeError(f"The tabular transformers {missing} have no row at their tap_pos in "
                           "pp_net.trafo_characteristic_table.")
    return tap_table


def _ratio_column(pp_net, name):
    """a ratio column of pp_net.trafo, 0.5 where missing (pandapower's default)"""
    if name not in pp_net.trafo:
        return np.full(pp_net.trafo.shape[0], 0.5)
    vals = pp_net.trafo[name].to_numpy(dtype=float).copy()
    vals[~np.isfinite(vals)] = 0.5
    return vals


def _wye_delta(r, x, y, r_ratio, x_ratio):
    """pandapower's `_wye_delta`: the pi model of a T one (series r + jx split `r_ratio`,
    `x_ratio` on the hv side, magnetising admittance `y` in the middle). Returns r, x and the shunt
    admittance of each side."""
    r, x = r.astype(float).copy(), x.astype(float).copy()
    h_1 = np.zeros(r.shape[0], dtype=complex)
    h_2 = np.zeros(r.shape[0], dtype=complex)
    tidx = y != 0.
    za_star = r[tidx] * r_ratio[tidx] + 1j * x[tidx] * x_ratio[tidx]
    zb_star = r[tidx] * (1. - r_ratio[tidx]) + 1j * x[tidx] * (1. - x_ratio[tidx])
    zc_star = 1. / y[tidx]
    z_sum = za_star * zb_star + za_star * zc_star + zb_star * zc_star
    zab = z_sum / zc_star
    r[tidx], x[tidx] = zab.real, zab.imag
    h_1[tidx] = zb_star / z_sum  # 1 / zac
    h_2[tidx] = za_star / z_sum  # 1 / zbc
    return r, x, h_1, h_2


def _complex_tap(vn, step_pct, diff, step_degree):
    """pandapower's (non ideal) tap, `_get_trafo_shift`: the tap side voltage `vn` moved by `diff`
    steps of `step_pct` % along `step_degree`. Returns the modulus of the result and its angle
    (degree)."""
    du = vn * 0.01 * step_pct * diff
    angle = np.deg2rad(step_degree)
    v_re = vn + du * np.cos(angle)
    v_im = du * np.sin(angle)
    return np.sqrt(v_re ** 2 + v_im ** 2), np.rad2deg(np.arctan(v_im / v_re))


def _ideal_tap(step_pct, diff, step_degree):
    """pandapower's ideal tap, `_get_trafo_shift(ideal=True)`: the phase shift (degree) of `diff`
    steps of `step_degree`, or, without one, of the chord `step_pct` % of the unit circle."""
    return np.where(step_degree != 0., diff * step_degree,
                    2. * np.rad2deg(np.arcsin(0.01 * diff * step_pct / 2.)))


def _pct_change(new, ref):
    """100 * (new / ref - 1), 0 where ref is 0 (a quantity that is 0 at the current tap is
    0 at every tap for this converter: r, x and the magnetising admittance only scale)."""
    new = np.asarray(new, dtype=float)
    ref = np.asarray(ref, dtype=float)
    out = np.zeros(new.shape)
    nz = ref != 0.
    out[nz] = 100. * (new[nz] / ref[nz] - 1.)
    return out


def _aux_trafo_tap_changers(model, pp_net, tap1, tap_table, trafo_params, trafo_r, trafo_x, trafo_b):
    """Give a tap changer to every transformer whose ``tap_min`` .. ``tap_max`` describe one: one
    step per integer position, with the ratio, the phase shift and the r, x, g, b this converter
    computes at that position. Moving the tap is then the same as converting the grid again with
    ``tap_pos`` at the new position.

    It is a ratio tap changer (``LSGrid.set_trafo_ratio_tap_changer``, what the
    TransformerVoltageControl outer loop moves) when only the ratio moves with the tap, and a phase
    tap changer (``LSGrid.set_trafo_phase_tap_changer``, with a ratio per step too) when the phase
    shift moves as well: an asymmetrical or symmetrical tap (``tap_step_degree`` != 0), an ideal
    phase shifter, a tabular tap changer whose table has an angle.

    Left without one (and so with the fixed ratio and shift of its current tap):

    - a transformer whose tap pandapower ignores, or with no step (``tap_step_percent`` and, for an
      ideal one, ``tap_step_degree`` 0 or NaN);
    - no finite ``tap_min`` / ``tap_max``, or a tap position that is not an integer within them;
    - a tabular tap changer whose table misses a position within ``tap_min`` .. ``tap_max``;
    - the second tap changer (``tap2_*``), which stays at its position.

    Returns the boolean mask of the transformers given a ratio tap changer.
    """
    nb_trafo = pp_net.trafo.shape[0]
    with_steps = np.zeros(nb_trafo, dtype=bool)
    if nb_trafo == 0 or "tap_min" not in pp_net.trafo or "tap_max" not in pp_net.trafo:
        return with_steps
    tap_min = pp_net.trafo["tap_min"].to_numpy(dtype=float)
    tap_max = pp_net.trafo["tap_max"].to_numpy(dtype=float)
    tap_pos = tap1.pos
    with np.errstate(invalid="ignore"):
        ok = (np.isfinite(tap_min) & np.isfinite(tap_max) &
              (tap_min == np.round(tap_min)) & (tap_max == np.round(tap_max)) & (tap_min <= tap_max) &
              (tap_pos == np.round(tap_pos)) & (tap_pos >= tap_min) & (tap_pos <= tap_max) &
              (tap1.neutral == np.round(tap1.neutral)))
    moves_ratio = tap1.complex & (tap1.step_pct != 0.)
    moves_phase = (moves_ratio & (tap1.step_degree != 0.)) | \
        (tap1.ideal & ((tap1.step_pct != 0.) | (tap1.step_degree != 0.)))
    for tid in np.flatnonzero(ok & tap1.tabular):
        positions = range(int(tap_min[tid]), int(tap_max[tid]) + 1)
        if not all(tap_table.has(tid, pos) for pos in positions):
            continue
        moves_ratio[tid] = True
        angles, = tap_table.at(np.full(len(positions), tid), np.array(positions), ("angle_deg",), (0.,))
        moves_phase[tid] = bool(np.any(angles != 0.))
    ok &= moves_ratio | moves_phase
    trafo_ids = np.flatnonzero(ok)
    if trafo_ids.shape[0] == 0:
        return with_steps

    # every (transformer, position) at once, through the same converter as the grid itself
    rep_ids, rep_pos = [], []
    for tid in trafo_ids:
        positions = np.arange(int(tap_min[tid]), int(tap_max[tid]) + 1)
        rep_ids.append(np.full(positions.shape[0], tid))
        rep_pos.append(positions)
    rep_ids = np.concatenate(rep_ids)
    rep_pos = np.concatenate(rep_pos)
    r_all, x_all, h1_all, h2_all, rho_all, shift_all = trafo_params(rep_ids, rep_pos)
    b_all = h1_all + h2_all  # its change: the same on both sides, both scale with the tap

    start = 0
    for tid in trafo_ids:
        n_pos = int(tap_max[tid]) - int(tap_min[tid]) + 1
        sl = slice(start, start + n_pos)
        start += n_pos
        # the corrections apply to the r, x, h the transformer was built with: at its current tap
        corrections = [list(_pct_change(r_all[sl], trafo_r[tid])),
                       list(_pct_change(x_all[sl], trafo_x[tid])),
                       list(_pct_change(b_all[sl].real, np.real(trafo_b[tid]))),
                       list(_pct_change(b_all[sl].imag, np.imag(trafo_b[tid])))]
        if moves_phase[tid]:
            # alpha: the whole shift at each step (shift_degree included), as the transformer's own
            model.set_trafo_phase_tap_changer(int(tid), int(tap_min[tid]), int(tap_pos[tid]),
                                              list(rho_all[sl]), list(shift_all[sl]), *corrections)
        else:
            model.set_trafo_ratio_tap_changer(int(tid), int(tap_min[tid]), int(tap_pos[tid]),
                                              list(rho_all[sl]), *corrections)
            with_steps[tid] = True
    return with_steps


def _tap_controllers(pp_net):
    """The pandapower tap controllers of ``pp_net.controller`` that are in service, as
    ``(controller, is_continuous)``."""
    if "controller" not in pp_net or pp_net.controller.shape[0] == 0:
        return []
    try:
        from pandapower.control import DiscreteTapControl, ContinuousTapControl
    except ImportError:  # pragma: no cover - pandapower without its control package
        return []
    res = []
    for obj, in_service in zip(pp_net.controller["object"], pp_net.controller["in_service"]):
        if not in_service:
            continue
        if isinstance(obj, ContinuousTapControl):
            res.append((obj, True))
        elif isinstance(obj, DiscreteTapControl):
            res.append((obj, False))
    return res


def _aux_trafo_tap_regulation(model, pp_net, pp_to_ls, with_steps):
    """What the ratio tap changers regulate (``LSGrid.set_trafo_ratio_tap_regulation``), read
    from pandapower's ``DiscreteTapControl`` and ``ContinuousTapControl`` controllers: the
    voltage of the transformer's own ``side`` bus, at

    - ``vm_set_pu`` with a deadband of twice ``tol`` for a continuous one;
    - the middle of ``[vm_lower_pu, vm_upper_pu]`` with that band as deadband for a discrete
      one (``vm_set_pu`` when it was built ``from_tap_step_percent``, the same band).

    As pandapower does, a transformer whose regulated bus has an ext_grid does not regulate.
    Controllers of three-winding transformers are skipped (they are not converted), and so are
    those of a transformer without a tap changer (see ``_aux_trafo_tap_changers``).
    """
    controllers = _tap_controllers(pp_net)
    if not controllers:
        return
    ext_grid_bus = set()
    if "ext_grid" in pp_net and pp_net.ext_grid.shape[0]:
        eg = pp_net.ext_grid
        ext_grid_bus = set(eg.loc[eg["in_service"].astype(bool), "bus"].tolist())
    skipped = []
    for obj, continuous in controllers:
        if getattr(obj, "element", "trafo") != "trafo":
            continue
        el_index = np.atleast_1d(obj.element_index)
        n = el_index.shape[0]
        side = np.broadcast_to(np.atleast_1d(obj.side), n)
        if continuous:
            target = np.broadcast_to(np.atleast_1d(np.asarray(obj.vm_set_pu, dtype=float)), n)
            deadband = np.broadcast_to(2. * np.atleast_1d(np.asarray(obj.tol, dtype=float)), n)
        else:
            low = np.broadcast_to(np.atleast_1d(np.asarray(obj.vm_lower_pu, dtype=float)), n)
            high = np.broadcast_to(np.atleast_1d(np.asarray(obj.vm_upper_pu, dtype=float)), n)
            vm_set = getattr(obj, "vm_set_pu", None)
            if vm_set is not None:
                target = np.broadcast_to(np.atleast_1d(np.asarray(vm_set, dtype=float)), n)
            else:
                target = 0.5 * (low + high)
            deadband = high - low
        for k, label in enumerate(el_index):
            if label not in pp_net.trafo.index:
                continue
            tid = int(pp_net.trafo.index.get_loc(label))
            if not with_steps[tid]:
                skipped.append(label)
                continue
            pp_bus = pp_net.trafo[f"{side[k]}_bus"].iloc[tid]
            ls_bus = int(pp_bus_to_ls(np.array([pp_bus]), pp_to_ls)[0])
            regulating = bool(np.isfinite(target[k]) and pp_bus not in ext_grid_bus)
            model.set_trafo_ratio_tap_regulation(tid, regulating,
                                                 float(target[k]) if np.isfinite(target[k]) else 0.,
                                                 float(deadband[k]) if np.isfinite(deadband[k]) else 0.,
                                                 ls_bus)
    if skipped:
        warnings.warn(f"The tap controllers of the transformers {skipped} are ignored: these transformers have "
                      "no tap changer lightsim2grid can model (see `_aux_trafo_tap_changers`).")
