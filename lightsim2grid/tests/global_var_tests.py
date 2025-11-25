# Copyright (c) 2020-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

from packaging import version
import pandapower as pp
import pypowsybl as pypow

MAX_PP_DATAREADER_NOT_BROKEN = version.parse("2.16")  
MIN_PYPO_DC_NOT_WORKING = version.parse("1.13.0")

# minimum version of pandapower required to run the tests
# pandapower version with more advanced grid modelling
CURRENT_PP_VERSION = version.parse(pp.__version__)
CURRENT_PYPOW_VERSION = version.parse(pypow.__version__)
