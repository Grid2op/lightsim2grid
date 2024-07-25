# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

def pp_bus_to_ls(pp_bus_id, pp_to_ls_converter):
    if pp_to_ls_converter is None:
        res = pp_bus_id
    else:
        res = np.array([pp_to_ls_converter[pp_id] for pp_id in pp_bus_id])
    return res
