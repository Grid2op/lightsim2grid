# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

# `SolverType` itself is *not* deprecated to use (isinstance checks against it, e.g. in
# LightSimBackend's back-compat bridging, are legitimate and should not print a warning).
# It lives here, rather than in lightsim2grid.solver, precisely so that importing it does
# not trigger lightsim2grid.solver's module-level "you are using the deprecated module"
# DeprecationWarning. lightsim2grid.solver imports it from here and re-exports it, so
# `from lightsim2grid.solver import SolverType` still works (and still warns) exactly as
# before -- this module is just where the class is defined now.
#
# Deliberately NOT under lightsim2grid._utils (or any other subpackage): importing a
# submodule always executes its package's __init__.py first, and lightsim2grid._utils's
# __init__.py does `from grid2op.Backend import PandaPowerBackend` -- which is only safe
# lazily, well after grid2op/pandapower have finished importing (see the comment at the
# top of that file: "grid2op -> pandapower -> lightsim2grid -> grid2op"). Importing this
# module eagerly from LightSimBackend's own top level, when it lived inside
# lightsim2grid._utils, forced that package's __init__.py to run reentrantly *while*
# grid2op.Backend.pandaPowerBackend was itself still mid-import (pandapower's own import
# chain probes for lightsim2grid), raising "cannot import name 'PandaPowerBackend' from
# partially initialized module 'grid2op.Backend'" -- silently swallowed by that file's
# `except ImportError: pass`, so `_DoNotUseAnywherePandaPowerBackend` was never defined
# and every grid2op.make(..., backend=LightSimBackend()) failed. A plain top-level module
# (this one) has no such package __init__.py to reenter.
#
# This module is internal: do not import it from anywhere except lightsim2grid.solver and
# LightSimBackend's own back-compat bridging.

from enum import Enum

from lightsim2grid.algorithm import AlgorithmType


class SolverType(Enum):
    """For backward compatibility, please use
    `from lightsim2grid.algorithm import AlgorithmType` instead now.
    """
    GaussSeidel = AlgorithmType.GaussSeidel
    GaussSeidelSynch = AlgorithmType.GaussSeidelSynch

    SparseLU = AlgorithmType.NR_SparseLU
    SparseLUSingleSlack = AlgorithmType.NRSing_SparseLU
    DC = AlgorithmType.DC_SparseLU
    FDPF_XB_SparseLU = AlgorithmType.FDPF_XB_SparseLU
    FDPF_BX_SparseLU = AlgorithmType.FDPF_BX_SparseLU

    KLU = AlgorithmType.NR_KLU
    KLUSingleSlack = AlgorithmType.NRSing_KLU
    KLUDC = AlgorithmType.DC_KLU
    FDPF_XB_KLU = AlgorithmType.FDPF_XB_KLU
    FDPF_BX_KLU = AlgorithmType.FDPF_BX_KLU

    NICSLU = AlgorithmType.NR_NICSLU
    NICSLUSingleSlack = AlgorithmType.NR_NICSLU
    NICSLUDC = AlgorithmType.NR_NICSLU
    FDPF_XB_NICSLU = AlgorithmType.NR_NICSLU
    FDPF_BX_NICSLU = AlgorithmType.NR_NICSLU

    CKTSO = AlgorithmType.NR_CKTSO
    CKTSOSingleSlack = AlgorithmType.NR_CKTSO
    CKTSODC = AlgorithmType.NR_CKTSO
    FDPF_XB_CKTSO = AlgorithmType.NR_CKTSO
    FDPF_BX_CKTSO = AlgorithmType.NR_CKTSO
