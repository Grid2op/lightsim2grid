# Copyright (c) 2020-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import numpy as np
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
    if "parallel" in pp_net.trafo and np.any(pp_net.trafo["parallel"].to_numpy() != 1):
        raise RuntimeError("Cannot handle 'parallel' trafo columns. Please duplicate the rows if that is the case. "
                           "Some pp_net.trafo[\"parallel\"] != 1 it is not handled by lightsim yet.")

    # fix the missing values
    tap_neutral = pp_net.trafo["tap_neutral"].to_numpy().copy()
    if np.any(~np.isfinite(tap_neutral)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_neutral\"], they have been replaced by 0")
    tap_neutral[~np.isfinite(tap_neutral)] = 0.

    tap_step_pct = pp_net.trafo["tap_step_percent"].to_numpy().copy()
    if np.any(~np.isfinite(tap_step_pct)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_step_percent\"], they have been replaced by 0")
    tap_step_pct[~np.isfinite(tap_step_pct)] = 0.

    tap_pos = pp_net.trafo["tap_pos"].to_numpy().copy()
    if np.any(~np.isfinite(tap_pos)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_pos\"], they have been replaced by the "
                      "neutral position (pp_net.trafo[\"tap_neutral\"])")
    tap_pos[~np.isfinite(tap_pos)] = tap_neutral[~np.isfinite(tap_pos)]
    # pandapower counts the steps from the neutral position, and so does everything below
    tap_rel = tap_pos - tap_neutral

    shift_ = pp_net.trafo["shift_degree"].to_numpy().copy()
    if np.any(~np.isfinite(shift_)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"shift_degree\"], they have been replaced by 0")
    shift_[~np.isfinite(shift_)] = 0.

    # test the missing values on the column itself: `== "hv"` is already False for them
    tap_side_missing = pp_net.trafo["tap_side"].isna().to_numpy()
    if np.any(tap_side_missing):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_side\"], they have been replaced by \"hv\"")
    is_tap_hv_side = (pp_net.trafo["tap_side"].to_numpy() == "hv") | tap_side_missing
    
    if "tap_phase_shifter" in pp_net.trafo:
        if np.any(pp_net.trafo["tap_phase_shifter"].to_numpy()):
            raise RuntimeError("Ideal phase shifters are not modeled. Please remove all trafos with "
                               "pp_net.trafo[\"tap_phase_shifter\"] set to True.")
    elif "tap_changer_type" in pp_net.trafo:
        if np.any(pp_net.trafo["tap_changer_type"].to_numpy() == "Ideal"):
            raise RuntimeError("Ideal phase shifters are not modeled. Please remove all 2-winding trafos "
                               "with \"tap_changer_type\" set to \"Ideal\".")
    elif "tap_changer_type" in pp_net.trafo3w:
        if np.any(pp_net.trafo3w["tap_changer_type"].to_numpy() == "Ideal"):
            raise RuntimeError("Ideal phase shifters are not modeled. Please remove all 3-winding trafos "
                               "with \"tap_changer_type\" set to \"Ideal\".")
            
    tap_angles_ = pp_net.trafo["tap_step_degree"].to_numpy().copy()
    if np.any(~np.isfinite(tap_angles_)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_step_degree\"], they have been replaced by 0")
    tap_angles_[~np.isfinite(tap_angles_)] = 0.
    tap_angles_ = np.deg2rad(tap_angles_)

    if "leakage_resistance_ratio_hv" in pp_net.trafo and (np.abs(pp_net.trafo["leakage_resistance_ratio_hv"].to_numpy() - 0.5) > 1e-7).any():
        warnings.warn("leakage_resistance_ratio_hv != 0.5 is not supported by this converter at the moment. It will be replaced by 0.5")
        
    trafo_model_is_t = True
    if "_options" in pp_net and "trafo_model" in pp_net._options:
        trafo_model_is_t = pp_net._options["trafo_model"] == "t"
    
    # compute physical parameters
    use_pp3 = version.parse(pp.__version__) >= _MIN_PP_VERSION_ADV_GRID_MODEL and pp_orig_file == "pandapower_v3"
    if not use_pp3 and not trafo_model_is_t:
        raise RuntimeError("Cannot convert a transformer with model 'pi' to LightSim2grid (using pandapower < 3)")
    per_trafo = (tap_step_pct,
                 tap_angles_,  # in radian !
                 is_tap_hv_side,
                 pp_net.bus.loc[pp_net.trafo["hv_bus"]]["vn_kv"].to_numpy(),
                 pp_net.bus.loc[pp_net.trafo["lv_bus"]]["vn_kv"].to_numpy(),
                 pp_net.trafo["vk_percent"].to_numpy(),
                 pp_net.trafo["vkr_percent"].to_numpy(),
                 pp_net.trafo["sn_mva"].to_numpy(),
                 pp_net.trafo["pfe_kw"].to_numpy(),
                 pp_net.trafo["i0_percent"].to_numpy())

    def trafo_params(trafo_ids, rel_positions):
        """r, x, b (pu) of the transformers `trafo_ids` with their taps `rel_positions`
        steps away from neutral -- the one place they are computed, so that a tap
        position the step tables below describe is exactly the transformer this
        converter would build at that position."""
        step, angle, hv_side, vn_hv, vn_lv, vk, vkr, sn, pfe, i0 = [np.asarray(el)[trafo_ids] for el in per_trafo]
        rel_positions = np.asarray(rel_positions, dtype=float)
        if use_pp3:
            # use pandapower version 3 converter in this case.
            # We use it because:
            # - the grid comes from pandapower3 (eg pn.case118() with pandapower 3 installed)
            # - AND the pandapower version is >= 3
            return converter.get_trafo_param_pp3(step, rel_positions, angle, list(hv_side), vn_hv, vn_lv,
                                                 vk, vkr, sn, pfe, i0, trafo_model_is_t)
        # default legacy mode: use the pandapower 2 converter
        return converter.get_trafo_param_pp2(step, rel_positions, angle, list(hv_side), vn_hv, vn_lv,
                                             vk, vkr, sn, pfe, i0)

    all_trafos = np.arange(pp_net.trafo.shape[0])
    trafo_r, trafo_x, trafo_b = trafo_params(all_trafos, tap_rel)

    # initialize the grid
    do_ignore_tap_side_for_phase_shift = True
    model.init_trafo_pandapower(trafo_r,
                                trafo_x,
                                trafo_b,
                                tap_step_pct,
                                tap_rel,
                                shift_,
                                is_tap_hv_side,
                                pp_bus_to_ls(pp_net.trafo["hv_bus"].values, pp_to_ls),
                                pp_bus_to_ls(pp_net.trafo["lv_bus"].values, pp_to_ls),
                                do_ignore_tap_side_for_phase_shift)

    for tr_id, is_connected in enumerate(pp_net.trafo["in_service"].values):
        if not is_connected:
            # trafo is deactivated
            model.deactivate_trafo(tr_id)

    # the tap changers (what the TransformerVoltageControl outer loop moves), then what they regulate
    with_steps = _aux_trafo_tap_changers(model, pp_net, tap_step_pct, tap_angles_, tap_neutral, tap_pos,
                                         trafo_params, trafo_r, trafo_x, trafo_b)
    _aux_trafo_tap_regulation(model, pp_net, pp_to_ls, with_steps)


def _pct_change(new, ref):
    """100 * (new / ref - 1), 0 where ref is 0 (a quantity that is 0 at the current tap is
    0 at every tap for this converter: r, x and the magnetising admittance only scale)."""
    new = np.asarray(new, dtype=float)
    ref = np.asarray(ref, dtype=float)
    out = np.zeros(new.shape)
    nz = ref != 0.
    out[nz] = 100. * (new[nz] / ref[nz] - 1.)
    return out


def _aux_trafo_tap_changers(model, pp_net, tap_step_pct, tap_angles_rad, tap_neutral, tap_pos,
                            trafo_params, trafo_r, trafo_x, trafo_b):
    """Give a ratio tap changer (``LSGrid.set_trafo_ratio_tap_changer``) to every transformer
    whose ``tap_min`` .. ``tap_max`` describe one: one step per integer position, with the
    ratio and the r, x, g, b this converter computes at that position. Moving the tap is then
    the same as converting the grid again with ``tap_pos`` at the new position.

    Left without one (and so with the fixed ratio of its current tap, as before):

    - a transformer with no step (``tap_step_percent`` 0 or NaN) or no finite ``tap_min`` / ``tap_max``;
    - a tap position that is not an integer within ``[tap_min, tap_max]``;
    - a non-zero ``tap_step_degree``: pandapower moves the phase shift with the tap there, which
      this converter does not model;
    - a tabular tap changer (``tap_dependency_table``), whose table is not read.

    Returns the boolean mask of the transformers given one.
    """
    nb_trafo = pp_net.trafo.shape[0]
    with_steps = np.zeros(nb_trafo, dtype=bool)
    if nb_trafo == 0 or "tap_min" not in pp_net.trafo or "tap_max" not in pp_net.trafo:
        return with_steps
    tap_min = pp_net.trafo["tap_min"].to_numpy(dtype=float)
    tap_max = pp_net.trafo["tap_max"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        ok = (np.isfinite(tap_min) & np.isfinite(tap_max) & (tap_step_pct != 0.) &
              (tap_min == np.round(tap_min)) & (tap_max == np.round(tap_max)) & (tap_min <= tap_max) &
              (tap_pos == np.round(tap_pos)) & (tap_pos >= tap_min) & (tap_pos <= tap_max) &
              (tap_neutral == np.round(tap_neutral)))
    ok &= tap_angles_rad == 0.
    if "tap_dependency_table" in pp_net.trafo:
        ok &= ~pp_net.trafo["tap_dependency_table"].fillna(False).to_numpy(dtype=bool)
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
    rel = rep_pos - tap_neutral[rep_ids]
    r_all, x_all, b_all = trafo_params(rep_ids, rel)
    b_all = np.asarray(b_all)
    # what init_trafo_pandapower does with the ratio: 1 + step % . (position - neutral)
    rho_all = 1. + 0.01 * tap_step_pct[rep_ids] * rel

    start = 0
    for tid in trafo_ids:
        n_pos = int(tap_max[tid]) - int(tap_min[tid]) + 1
        sl = slice(start, start + n_pos)
        start += n_pos
        # the corrections apply to the r, x, h the transformer was built with: at its current tap
        model.set_trafo_ratio_tap_changer(int(tid), int(tap_min[tid]), int(tap_pos[tid]),
                                          list(rho_all[sl]),
                                          list(_pct_change(r_all[sl], trafo_r[tid])),
                                          list(_pct_change(x_all[sl], trafo_x[tid])),
                                          list(_pct_change(b_all[sl].real, np.real(trafo_b[tid]))),
                                          list(_pct_change(b_all[sl].imag, np.imag(trafo_b[tid]))))
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
