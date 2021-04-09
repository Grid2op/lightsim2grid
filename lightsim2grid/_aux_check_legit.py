# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
import warnings
import numpy as np


def _aux_check_legit(pp_net):
    """
    Check that all the dataframes of the pandapower network can be handled in lightsim2grid.

    Parameters
    ----------
    pp_net

    Returns
    -------

    """

    if "trafo3w" in pp_net and pp_net.trafo3w.shape[0]:
        raise RuntimeError("Unsupported element found (Three Winding Transformer - \"pp_net.trafo3w\") "
                           "in pandapower network")
    if "switch" in pp_net and pp_net.switch.shape[0]:
        warnings.warn("There are switches on the pandapower network, they will not be used in lightsim2grid.")
    if "motor" in pp_net and pp_net.motor.shape[0]:
        raise RuntimeError("Unsupported element found (motor - \"pp_net.motor\") in pandapower network")
    if "asymmetric_load" in pp_net and pp_net.asymmetric_load.shape[0]:
        raise RuntimeError("Unsupported element found (Asymmetric Load - \"pp_net.asymmetric_load\") "
                           "in pandapower network")
    if "impedance" in pp_net and pp_net.impedance.shape[0]:
        raise RuntimeError("Unsupported element found (Impedance - \"pp_net.impedance\") "
                           "in pandapower network")
    if "ward" in pp_net and pp_net.ward.shape[0]:
        raise RuntimeError("Unsupported element found (Ward - \"pp_net.ward\") "
                           "in pandapower network")
    if "xward" in pp_net and pp_net.xward.shape[0]:
        raise RuntimeError("Unsupported element found (Extended Ward - \"pp_net.xward\") "
                           "in pandapower network")
    if "dcline" in pp_net and pp_net.dcline.shape[0]:
        raise RuntimeError("Unsupported element found (DC Line - \"pp_net.dcline\") "
                           "in pandapower network")

    # if pp_net.sn_mva != 1.:
    #     warnings.warn("Pandapower network with sn_mva != 1 are not fully supported yet, especially in the "
    #                   "conversion for some trafo (and probably somewhere else too)")
    if "_options" in pp_net and \
       "trafo_model" in pp_net["_options"] and \
       pp_net["_options"]["trafo_model"] != "t":
        raise RuntimeError("Cannot handle pp_net with trafo not modeled in \"t\". Please change "
                           "your network and load it again.")
