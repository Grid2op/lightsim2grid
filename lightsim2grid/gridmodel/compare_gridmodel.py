# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

from functools import partial
import numpy as np


from lightsim2grid_cpp import GridModel


ATTR_SUBSTATIONS_INPUT = [
    "name",
    "id",
    "nb_max_busbars",
    "vn_kv"
]


ATTR_1SIDE_INPUT = [
    "id",
    "name",
    "sub_id",
    "pos_topo_vect",
    "connected",
    "bus_id",
    "target_p_mw",
]


ATTR_GENS_INPUT = [
    "is_slack",
    'slack_weight',
    "voltage_regulator_on",
    "target_vm_pu",
    "target_p_mw",
    "target_q_mvar",
    "min_q_mvar",
    "max_q_mvar"
]


ATTR_SGENS_INPUT = [
    "min_p_mw",
    "max_p_mw",
    "target_p_mw",
    "target_q_mvar",
    "min_q_mvar",
    "max_q_mvar"
]


ATTR_LOADS_INPUT = [
    "target_p_mw",
    "target_q_mvar",
]


ATTR_STORAGES_INPUT = ATTR_LOADS_INPUT


ATTR_SHUNTS_INPUT = [
    "target_p_mw",
    "target_q_mvar",
]


ATTR_2SIDES_INPUT = [
    "name",
    "sub1_id",
    "sub2_id",
    "pos1_topo_vect",
    "pos2_topo_vect",
    "connected_global",
    "connected1",
    "connected2",
    "bus1_id",
    "bus2_id",
]


ATTR_LINES_INPUT = (ATTR_2SIDES_INPUT + 
    [
         "r_pu",
         "x_pu",
         "h1_pu",
         "h2_pu",
         "yac_11",
         "yac_12",
         "yac_21",
         "yac_22",
         "ydc_11",
         "ydc_12",
         "ydc_21",
         "ydc_22",
    ]
)


ATTR_TAFO_INPUT = (ATTR_LINES_INPUT + 
    [
        "is_tap_hv_side",
        "ratio",
        "shift_rad",
    ]
)


ATTR_DCLINE_INPUT = (ATTR_2SIDES_INPUT + 
    [
        "target_p1_mw",
        "p2_mw",
        "target_vm1_pu",
        "target_vm2_pu",
        "loss_pct",
        "loss_mw",
    ]
)


def __are_float_different(attr1, attr2, tol) -> bool:
    return np.abs(attr1 - attr2) > tol


def __are_attr_different(attr1, attr2) -> bool:
    return attr1 != attr2


def _aux_compare_one_el(tmp, li_attrs, el1, el2, tol):
    for attr_nm in li_attrs:
        attr1 = getattr(el1, attr_nm)
        attr2 = getattr(el2, attr_nm)
        if isinstance(attr1, (float, np.float16, np.float32, np.float64, complex, np.complex64, np.complex128)):
            fun = partial(__are_float_different, tol=tol)
        else:
            fun = __are_attr_different
        if fun(attr1, attr2):
            tmp[attr_nm] = attr1, attr2
    
    
def _aux_compare(
    gridmodel1: GridModel,
    gridmodel2: GridModel,
    meth_nm,
    li_attrs,
    tol):
    res = {}
    if len(getattr(gridmodel1, meth_nm)()) != len(getattr(gridmodel2, meth_nm)()):
        res["size"] = len(getattr(gridmodel1, meth_nm)()), len(getattr(gridmodel2, meth_nm)())
        return res
    
    for num, (el1, el2) in enumerate(zip(getattr(gridmodel1, meth_nm)(), getattr(gridmodel2, meth_nm)())):
        tmp = {}
        _aux_compare_one_el(tmp, li_attrs, el1, el2, tol)
        if len(tmp) > 0: 
            res[f"{num}"] = tmp
    return res
    
        
def _compare_substations(gridmodel1: GridModel, gridmodel2: GridModel, tol=1e-8):
    res_sub = _aux_compare(
        gridmodel1,
        gridmodel2,
        "get_substations",
        ATTR_SUBSTATIONS_INPUT,
        tol,
    )
    return res_sub
    
    
def _compare_lines(gridmodel1: GridModel, gridmodel2: GridModel, tol=1e-8):
    res_lines =  _aux_compare(
        gridmodel1,
        gridmodel2,
        "get_lines",
        ATTR_LINES_INPUT,
        tol,
    )
    return res_lines
    
    
def _compare_trafos(gridmodel1: GridModel, gridmodel2: GridModel, tol=1e-8):
    res_trafos =  _aux_compare(
        gridmodel1,
        gridmodel2,
        "get_trafos",
        ATTR_TAFO_INPUT,
        tol
    )
    return res_trafos


def _compare_dclines(gridmodel1: GridModel, gridmodel2: GridModel, tol=1e-8):
    meth_nm = "get_dclines"
    res_dclines =  _aux_compare(
        gridmodel1,
        gridmodel2,
        meth_nm,
        ATTR_DCLINE_INPUT,
        tol,
    )
    if "size" in res_dclines:
        return res_dclines

    for num, (el1, el2) in enumerate(zip(getattr(gridmodel1, meth_nm)(), getattr(gridmodel2, meth_nm)())):
        gen1 = el1.gen_side_1
        gen2 = el2.gen_side_1
        tmp = {}
        _aux_compare_one_el(tmp, ATTR_GENS_INPUT, gen1, gen2, tol)
        if len(tmp) > 0:
            res_dclines[f"{num}_gen_side_1"] = tmp
            
        gen1 = el1.gen_side_2
        gen2 = el2.gen_side_2
        tmp = {}
        _aux_compare_one_el(tmp, ATTR_GENS_INPUT, gen1, gen2, tol)
        if len(tmp) > 0:
            res_dclines[f"{num}_gen_side_2"] = tmp
        
    return res_dclines


def _compare_generators(gridmodel1: GridModel, gridmodel2: GridModel, tol=1e-8):
    res_gens =  _aux_compare(
        gridmodel1,
        gridmodel2,
        "get_generators",
        ATTR_GENS_INPUT,
        tol
    )
    return res_gens


def _compare_static_generators(gridmodel1: GridModel, gridmodel2: GridModel, tol=1e-8):
    res_sgens =  _aux_compare(
        gridmodel1,
        gridmodel2,
        "get_static_generators",
        ATTR_SGENS_INPUT,
        tol,
    )
    return res_sgens


def _compare_loads(gridmodel1: GridModel, gridmodel2: GridModel, tol=1e-8):
    res_loads =  _aux_compare(
        gridmodel1,
        gridmodel2,
        "get_loads",
        ATTR_LOADS_INPUT,
        tol,
    )
    return res_loads


def _compare_storages(gridmodel1: GridModel, gridmodel2: GridModel, tol=1e-8):
    res_sto =  _aux_compare(
        gridmodel1,
        gridmodel2,
        "get_storages",
        ATTR_STORAGES_INPUT,
        tol,
    )
    return res_sto


def _compare_shunts(gridmodel1: GridModel, gridmodel2: GridModel, tol=1e-8):
    res_sto =  _aux_compare(
        gridmodel1,
        gridmodel2,
        "get_shunts",
        ATTR_SHUNTS_INPUT,
        tol,
    )
    return res_sto


def compare_gridmodel_input(gridmodel1: GridModel, gridmodel2: GridModel, tol=1e-8):
    """
    This function tests that the two gridmodels as argument have the same underlying grid.
    
    This means the same structure (same elements, connected to the same substations, with the same 
    physical properties) and also the same "setpoints" (*eg* same target active and reactive power 
    for every loads, same voltage setpoints for each generators etc.)
    
    .. warning::
        It does not compare the output of a powerflows: no check is done for the "results" of powerflow,
        such as "res_p_mw" for generators.
        
    """
    res = {}
    subs = _compare_substations(gridmodel1, gridmodel2, tol)
    if len(subs) > 0:
        res["substations"] = subs
    lines = _compare_lines(gridmodel1, gridmodel2, tol)
    if len(lines) > 0:
        res["lines"] = lines
    trafos = _compare_trafos(gridmodel1, gridmodel2, tol)
    if len(trafos) > 0:
        res["trafos"] = trafos
    dclines = _compare_dclines(gridmodel1, gridmodel2, tol)
    if len(dclines) > 0:
        res["dclines"] = dclines
    gens = _compare_generators(gridmodel1, gridmodel2, tol)
    if len(gens) > 0:
        res["generators"] = gens
    sgens = _compare_static_generators(gridmodel1, gridmodel2, tol)
    if len(gens) > 0:
        res["static_generators"] = sgens
    loads = _compare_loads(gridmodel1, gridmodel2, tol)
    if len(loads) > 0:
        res["loads"] = loads
    stos = _compare_storages(gridmodel1, gridmodel2, tol)
    if len(stos) > 0:
        res["storages"] = stos
    shunts = _compare_shunts(gridmodel1, gridmodel2, tol)
    if len(shunts) > 0:
        res["shunts"] = shunts
    return res
