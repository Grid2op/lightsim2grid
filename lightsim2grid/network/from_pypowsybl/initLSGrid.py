# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Use the pypowsybl converter to properly initialize a LSGrid c++ object.

``init`` orchestrates one phase per element type, each factored into its own
``_aux_add_*.py`` module (mirroring the ``from_pandapower`` / ``from_powermodels``
converters' layout): buses/substations first (every other phase depends on its
``voltage_levels``/``bus_df``/``first_bus_per_vl`` output to resolve its own
elements' buses), then generators, loads, lines, transformers, shunts, SVCs,
HVDC lines, storage units, and finally slack assignment.
"""

from typing import Dict, Iterable, Optional, Union

import numpy as np
import pandas as pd
import pypowsybl as pypo

from ...lightsim2grid_cpp import LSGrid  # type: ignore

from ._aux_common import _aux_ensure_net_pu, _aux_get_operational_limits_current
from ._aux_add_buses import _aux_add_buses
from ._aux_add_generators import _aux_add_generators
from ._aux_add_loads import _aux_add_loads
from ._aux_add_lines import _aux_add_lines
from ._aux_add_trafos import _aux_add_trafos
from ._aux_add_shunts import _aux_add_shunts
from ._aux_add_svc import _aux_add_svc
from ._aux_add_hvdc import _aux_add_hvdc
from ._aux_add_storage import _aux_add_storage
from ._aux_add_slack import _aux_add_slack


def init(net : pypo.network.Network,
         gen_slack_id: Union[int, str, Iterable[str], Dict[str, float]] = None,
         slack_bus_id: int = None,
         sn_mva : float = 100.,  # only used if not present in the grid
         sort_index : bool =True,
         f_hz : float = 50.,  # unused
         net_pu : Optional[pypo.network.Network] = None,
         only_main_component : bool =True,
         return_sub_id: bool=False,
         n_busbar_per_sub: Optional[int]=None,  # new in 0.9.1
         buses_for_sub:Optional[bool]=None,  # new in 0.9.1
         init_vm_pu:float=1.06,
         keep_half_open_lines: bool=False,
         convert_dangling_lines: bool=False,
         fuse_zero_impedance_branches: bool=False,
         zero_impedance_threshold_pu: float=1e-8,
         ) -> LSGrid:
    """
    This function is available under the `init_from_pypowsybl` in lightsim2grid


    .. code-block:: python

        from lightsim2grid.network import init_from_pypowsybl

    .. warning::
        It is not available if the `pypowsybl` python package is not installed.

    :param net: The pypowsybl network
    :type net: pypo.network.Network

    :param gen_slack_id: The id of the generator that should be used as the slack
                         (either it's given by id (int) or by name (str))
    :type gen_slack_id: Union[int, str]

    :param slack_bus_id: If you don't provide a generator ID as a slack bus, you can
            provide a bus id (int). We do not recommend setting the slack this way.
    :type slack_bus_id: int

    :param sn_mva: The nominal apparent power used when converting the grid to
                   per unit. It is only used if the pypowsybl grid
                   has no `_nominal_apparent_power` attribute.
                   **Advanced usage**.
    :type sn_mva: float

    :param sort_index: Whether you want to sort the indexes of all the
                       pypowsybl tables (*eg* get_loads() or *get_buses()*) or not.
                       Sorting the grid tables is preferable if you want to be
                       "future proof" and don't want to depend on pandas version
                       (same order is guaranteed). Not sorting the grid will give
                       easier comparison of results with pypowsybl.
    :type sn_mva: bool

    :param f_hz: Not used currently (frequency of the grid)
    :type net_pu: float

    :param net_pu: If you have already converted the grid in "per unit" then
                   you can pass it as the `net_pu` argument. Otherwise this
                   function will do it.
                   **Advanced usage**.
    :type net_pu: Optional[pypo.network.Network]

    :param only_main_component: If this is True, then only the main component (*ie*
                                the one containing the slack bus) will be used. All
                                equipments not part of this component will be
                                deactivated (switched-off). **NB** currently
                                lightsim2grid will diverge if the grid is not connected,
                                this option might then "hide" some equipements from
                                the grid (silently) but you have higher chances of
                                convergence.
    :type only_main_component: bool

    :param return_sub_id: **Advanced usage**. If you want to retrieve the id of the
                          equipments as "tables". Used only for `LightSimBackend`
    :type return_sub_id: bool

    :param n_busbar_per_sub: Currently, lightsim2grid works well with a constant
                             number of independant buses that can be made at each
                             substations. It can be infered from the grid or
                             set with this attribute. We recommend to leave it
                             to `None` (which corresponds to the "infer it from
                             the grid" behaviour) in most cases.
    :type n_busbar_per_sub: Opional[int]

    :param buses_for_sub: Whether the lightsim2grid substation will correspond to buses
                          of the pypowsybl grid (if buses_for_sub is `True`).
                          Alternatively, if buses_for_sub is `False`, the
                          lightsim2grid susbtation will correspond to
                          pypowsybl voltage level (read from net.get_voltage_levels()).
                          buses_for_sub==`True` is a "legacy" behaviour.
    :type buses_for_sub: bool

    :param init_vm_pu: The voltage magnitude with which the init vector of AC powerflow
                       will be set.
    :type init_vm_pu: float

    :param keep_half_open_lines: If True, a powerline or transformer connected on only
                       one terminal (``connected1 != connected2``, *eg* a dangling
                       boundary stub in a real grid) is modeled as "half-open": the
                       energized side is kept in the admittance matrix and the open end
                       is Kron-reduced out, instead of deactivating the whole branch.
                       This sets ``synch_status_both_side=False`` on the returned model,
                       so a later one-sided topology change is no longer mirrored to the
                       other side. Branches disconnected on *both* sides are still fully
                       deactivated. Default False (whole-branch deactivation, as before).
    :type keep_half_open_lines: bool

    :param convert_dangling_lines: If True, every IIDM ``DanglingLine`` (*eg*
                       produced by ``network.reduce_by_ids_and_depths(...,
                       with_boundary_lines=True)`` when zooming into a
                       sub-area) is converted to its equivalent branch +
                       constant-power load: a fictitious 1-bus "substation" at
                       the boundary end, a line carrying the dangling line's own
                       r/x/g/b (shunt entirely on the local side, matching
                       transformers' single ``h``), and a load consuming its
                       p0/q0. Without this (the default), dangling lines are
                       silently ignored -- fine for real full grids, which
                       never have any, but it drops a real boundary injection
                       (can be hundreds of MW) whenever they do appear. Off by
                       default to keep existing behaviour unchanged; the
                       reduce/validate debug scripts turn it on.
    :type convert_dangling_lines: bool

    :param fuse_zero_impedance_branches: If True, a line or 2-winding transformer
                       whose per-unit impedance is (near-)zero -- ``|r_pu| <
                       zero_impedance_threshold_pu`` and ``|x_pu| <
                       zero_impedance_threshold_pu`` -- has its two terminal buses
                       fused into a single electrical node instead of contributing
                       a ``1/Z`` admittance (which is ``Inf`` for an exact zero, and
                       breaks the sparse LU factorization outright). This mirrors
                       OpenLoadFlow's ``lowImpedanceBranchMode`` /
                       ``lowImpedanceThreshold``. A zero-impedance transformer is
                       only fused if it is also at (near-)neutral tap (``rho`` close
                       to 1 and ``alpha`` close to 0) *and* its two sides are at the
                       same nominal voltage: pypowsybl's per-unit ``rho`` is the
                       deviation from the transformer's *own* rated ratio (the
                       tap-changer effect), not its absolute turns ratio, so
                       ``rho~=1`` alone does not mean "no transformation" for a
                       genuine step-down/up transformer -- otherwise it is a real
                       ideal ratio/phase-shifting element, not a same-node short, and
                       is left untouched. A zero-impedance *line* spanning two
                       *different* nominal voltages raises a ``RuntimeError``
                       (inconsistent grid data) -- unlike for transformers, this is
                       never legitimate for a line. The fusing branch itself is kept in
                       the model (for topology-vector bookkeeping) but deactivated,
                       same as any other disconnected branch; substation/topology
                       identity (``glop_sub_id``) of elements on the two original
                       buses is *not* changed, only the internal solver bus id --
                       grid2op topology actions still see the original, distinct
                       substations. Off by default to keep existing behaviour
                       unchanged.

                       .. warning::
                            This is a static, import-time decision. If a fused
                            transformer's tap is later moved away from neutral
                            during a simulation, the two buses stay fused in
                            lightsim2grid even though the transformer is no longer
                            electrically a plain wire. Best suited for static /
                            diagnostic use, or grid2op environments known not to
                            actuate the affected transformer's tap.
    :type fuse_zero_impedance_branches: bool

    :param zero_impedance_threshold_pu: Per-unit impedance magnitude threshold used
                       by ``fuse_zero_impedance_branches`` (ignored otherwise).
                       Matches OpenLoadFlow's ``lowImpedanceThreshold`` default.
    :type zero_impedance_threshold_pu: float

    :return: The properly initialized network.
    :rtype: :class:`LSGrid`
    """
    model = LSGrid()
    if hasattr(net, "_nominal_apparent_power"):
        sn_mva_used = getattr(net, "_nominal_apparent_power")
    else:
        sn_mva_used = float(sn_mva)
    model.set_sn_mva(sn_mva_used)
    model.set_init_vm_pu(float(init_vm_pu))
    if keep_half_open_lines:
        # allow branches connected on a single terminal (the open end is Kron-reduced
        # in the C++ model); a one-sided disconnection is no longer mirrored to the
        # other side.
        model.set_synch_status_both_side(False)

    if gen_slack_id is not None and slack_bus_id is not None:
        raise RuntimeError("Impossible to intialize a grid with both gen_slack_id and slack_bus_id")

    # buses / substations: every other phase below depends on this phase's
    # voltage_levels / bus_df / first_bus_per_vl to resolve its own elements' buses.
    voltage_levels, bus_df, first_bus_per_vl, df_dl, net_pu, fused_line_ids, fused_trafo_ids = _aux_add_buses(
        model, net, net_pu, sort_index, buses_for_sub, n_busbar_per_sub,
        convert_dangling_lines, fuse_zero_impedance_branches, zero_impedance_threshold_pu,
    )

    # generators
    df_gen, gen_sub = _aux_add_generators(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl)

    # loads
    df_load, load_sub = _aux_add_loads(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl, df_dl)

    # thermal (current) limits for lines / trafos, shared by both phases below
    ol_current = _aux_get_operational_limits_current(net)

    # per unit view, shared by both phases below (built once: the bus-fusion phase
    # above may already have built it if fuse_zero_impedance_branches=True)
    net_pu = _aux_ensure_net_pu(net, net_pu)

    # lines
    df_line, lor_sub, lex_sub = _aux_add_lines(
        model, net, net_pu, sort_index, voltage_levels, bus_df, first_bus_per_vl,
        df_dl, ol_current, keep_half_open_lines, fuse_zero_impedance_branches, fused_line_ids,
    )

    # trafos
    df_trafo, tor_sub, tex_sub = _aux_add_trafos(
        model, net, net_pu, sort_index, voltage_levels, bus_df, first_bus_per_vl,
        ol_current, keep_half_open_lines, fuse_zero_impedance_branches, fused_trafo_ids,
    )

    # shunts
    df_shunt, sh_sub = _aux_add_shunts(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl)

    # SVCs
    df_svc = _aux_add_svc(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl, sn_mva_used)

    # HVDC lines
    df_dc, hvdc_sub_from_id, hvdc_sub_to_id = _aux_add_hvdc(
        model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl,
    )

    # storage units
    df_batt, batt_sub = _aux_add_storage(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl)

    # slack bus(es)
    gen_slack_ids_int = _aux_add_slack(model, net, df_gen, gen_slack_id, slack_bus_id)

    # TODO checks
    # no 3windings trafo and other exotic stuff

    # and now deactivate all elements and nodes not in the main component
    if only_main_component:
        model.consider_only_main_component()
    else:
        # automatically disconnect non connected buses
        # (this is automatically done by consider_only_main_component)
        model.init_bus_status()

    gen_sub = pd.DataFrame(index=df_gen.index, data={"sub_id": gen_sub})
    gen_sub["desired_slack"] = False
    gen_sub.loc[gen_sub.index[gen_slack_ids_int], "desired_slack"] = True
    load_sub = pd.DataFrame(index=df_load.index, data={"sub_id": load_sub})
    lor_sub = pd.DataFrame(index=df_line.index, data={"sub_id": lor_sub})
    lex_sub = pd.DataFrame(index=df_line.index, data={"sub_id": lex_sub})
    tor_sub = pd.DataFrame(index=df_trafo.index, data={"sub_id": tor_sub})
    tex_sub = pd.DataFrame(index=df_trafo.index, data={"sub_id": tex_sub})
    batt_sub = pd.DataFrame(index=df_batt.index, data={"sub_id": batt_sub})
    sh_sub = pd.DataFrame(index=df_shunt.index, data={"sub_id": sh_sub})
    hvdc_sub_from_id = pd.DataFrame(index=df_dc.index, data={"sub_id": hvdc_sub_from_id})
    hvdc_sub_to_id = pd.DataFrame(index=df_dc.index, data={"sub_id": hvdc_sub_to_id})

    # set the substation ID to which each object belong
    model.set_gen_to_subid(gen_sub["sub_id"].values)
    model.set_load_to_subid(load_sub["sub_id"].values)
    model.set_storage_to_subid(batt_sub["sub_id"].values)
    model.set_shunt_to_subid(sh_sub["sub_id"].values)
    model.set_line_to_sub1_id(lor_sub["sub_id"].values)
    model.set_line_to_sub2_id(lex_sub["sub_id"].values)
    model.set_trafo_to_sub1_id(tor_sub["sub_id"].values)
    model.set_trafo_to_sub2_id(tex_sub["sub_id"].values)
    if not return_sub_id:
        return model
    else:
        return model, (gen_sub, load_sub, (lor_sub, tor_sub), (lex_sub, tex_sub), batt_sub, sh_sub, hvdc_sub_from_id, hvdc_sub_to_id)
