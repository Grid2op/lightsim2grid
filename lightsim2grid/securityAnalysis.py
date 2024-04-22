# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
# Deprecated now, will be removed
warnings.warn("You are using old names. Please upgrade to SecurityAnalysisCPP > ContingencyAnalysisCPP"
              " and SecurityAnalysis > ContingencyAnalysis instead.",
              category=DeprecationWarning)
__all__ = ["SecurityAnalysisCPP"]


from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP

try:
    from lightsim2grid.contingencyAnalysis import ContingencyAnalysis
    GRID2OP_INSTALLED = True
except ImportError as exc_:
    GRID2OP_INSTALLED = False


SecurityAnalysisCPP = ContingencyAnalysisCPP


if GRID2OP_INSTALLED:
    SecurityAnalysis = ContingencyAnalysis
    __all__.append("SecurityAnalysis")
