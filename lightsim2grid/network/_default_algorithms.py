# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

from lightsim2grid.algorithm import AlgorithmType


def use_fastest_default_algorithms(model):
    """Select ``NR_KLU`` and ``DC_KLU`` on a freshly loaded grid when this build has them.

    A grid starts on ``NR_SparseLU`` / ``DC_SparseLU``, which are always compiled in. KLU
    is not always (it is optional at build time), but where it is, it solves the same
    systems faster, and to the same answer up to rounding. The loaders call this last, so
    that a grid handed to a user, or copied by a batch class (which inherits its AC
    algorithm), is on the faster one without anybody having to ask.
    ``model.change_algorithm`` still selects anything else afterwards.
    """
    available = model.available_default_algorithms()
    if AlgorithmType.NR_KLU in available:
        model.change_algorithm(AlgorithmType.NR_KLU)
    if AlgorithmType.DC_KLU in available:
        model.change_algorithm(AlgorithmType.DC_KLU)
