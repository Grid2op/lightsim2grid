"""
Use the pandapower converter to properly initialized a GridModel c++ object.
"""

import numpy as np
from pyklu_cpp import GridModel, PandaPowerConverter


def init(pp_net):
    """
    Convert a pandapower network as input into a GridModel.

    This does not throw any error at the moment when the conversion is not possible.

    Cases for which conversion is not possible include, but are not limited to:

    - the pandapower grid has 3 winding transformers
    - the pandapower grid has xwards
    - the pandapower grid any parrallel "elements" (at least one of the column "parrallel" is not 1)
    - some `g_us_per_km` for some lines are not zero
    - some `p_mw` for some shunts are not zero
    - some `tap_step_degre` are non zero for some trafo
    - no "ext_grid" is reported on the initial grid

    if you really need any of the above, please submit a github issue and we will work on their support.

    This conversion has been extensively studied for the case118() of pandapower.networks and should work
    really well for this grid. Actually, this grid is used for testing the GridModel class.

    Parameters
    ----------
    pp_net: :class:`pandapower.grid`
        The initial pandapower network you want to convert

    Returns
    -------
    model: :class:`GridModel`
        The initialize gridmodel

    """
    # initialize and use converters
    converter = PandaPowerConverter()
    converter.set_sn_mva(pp_net.sn_mva)  # TODO raise an error if not set !
    converter.set_f_hz(pp_net.f_hz)
    line_r, line_x, line_h = \
        converter.get_line_param(
            pp_net.line["r_ohm_per_km"].values * pp_net.line["length_km"].values,
            pp_net.line["x_ohm_per_km"].values * pp_net.line["length_km"].values,
            pp_net.line["c_nf_per_km"].values * pp_net.line["length_km"].values,
            pp_net.line["g_us_per_km"].values * pp_net.line["length_km"].values,
            pp_net.bus.loc[pp_net.line["from_bus"]]["vn_kv"],
            pp_net.bus.loc[pp_net.line["to_bus"]]["vn_kv"]
        )
    trafo_r, trafo_x, trafo_b = \
        converter.get_trafo_param(pp_net.trafo["vn_hv_kv"].values,
                                  pp_net.trafo["vn_lv_kv"].values,
                                  pp_net.trafo["vk_percent"].values,
                                  pp_net.trafo["vkr_percent"].values,
                                  pp_net.trafo["sn_mva"].values,
                                  pp_net.trafo["pfe_kw"].values,
                                  pp_net.trafo["i0_percent"].values,
                                  pp_net.bus.loc[pp_net.trafo["lv_bus"]]["vn_kv"]
                                       )

    # set up the data model accordingly
    model = GridModel()
    tmp_bus_ind = np.argsort(pp_net.bus.index)
    model.init_bus(pp_net.bus.iloc[tmp_bus_ind]["vn_kv"].values,
                   pp_net.line.shape[0],
                   pp_net.trafo.shape[0])

    model.init_powerlines(line_r, line_x, line_h,
                          pp_net.line["from_bus"].values,
                          pp_net.line["to_bus"].values
                               )

    # init the shunts
    model.init_shunt(pp_net.shunt["p_mw"].values,
                     pp_net.shunt["q_mvar"].values,
                     pp_net.shunt["bus"].values
                          )

    tap_step_pct = pp_net.trafo["tap_step_percent"].values
    tap_step_pct[~np.isfinite(tap_step_pct)] = 0.

    tap_pos = pp_net.trafo["tap_pos"].values
    tap_pos[~np.isfinite(tap_pos)] = 0.

    is_tap_hv_side = pp_net.trafo["tap_side"].values == "hv"
    is_tap_hv_side[~np.isfinite(tap_pos)] = True
    model.init_trafo(trafo_r,
                     trafo_x,
                     trafo_b,
                     tap_step_pct,
                     tap_pos,
                     is_tap_hv_side,
                     pp_net.trafo["hv_bus"].values,
                     pp_net.trafo["lv_bus"].values)

    model.init_loads(pp_net.load["p_mw"].values,
                     pp_net.load["q_mvar"].values,
                     pp_net.load["bus"].values
                          )
    model.init_generators(pp_net.gen["p_mw"].values,
                          pp_net.gen["vm_pu"].values,
                          pp_net.gen["bus"].values
                               )

    # TODO handle that better maybe
    model.add_slackbus(pp_net.ext_grid["bus"].values[0])
    return model