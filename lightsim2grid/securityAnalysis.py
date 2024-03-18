# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP, ContingencyAnalysis

# Deprecated now, will be removed
SecurityAnalysisCPP = ContingencyAnalysisCPP
SecurityAnalysis = ContingencyAnalysis
