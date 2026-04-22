// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"

#ifndef KLU_SOLVER_AVAILABLE
#define this_KLU_SOLVER_AVAILABLE 0
#else
#define this_KLU_SOLVER_AVAILABLE 1
#endif
#ifndef NICSLU_SOLVER_AVAILABLE
#define this_NICSLU_SOLVER_AVAILABLE 0
#else
#define this_NICSLU_SOLVER_AVAILABLE 1
#endif
#ifndef CKTSO_SOLVER_AVAILABLE
#define this_CKTSO_SOLVER_AVAILABLE 0
#else
#define this_CKTSO_SOLVER_AVAILABLE 1
#endif
#ifndef __COMPILE_MARCHNATIVE
#define this__COMPILE_MARCHNATIVE 0
#else
#define this__COMPILE_MARCHNATIVE 1
#endif
#ifndef __O3_OPTIM
#define this__O3_OPTIM 0
#else
#define this__O3_OPTIM 1
#endif
#ifndef VERSION
#define this_VERSION "unknown"
#else
#define this_VERSION VERSION
#endif
#ifdef NICSLU_PATH
#define this_NICSLU_PATH NICSLU_PATH
#endif
#ifdef CKTSO_PATH
#define this_CKTSO_PATH CKTSO_PATH
#endif

PYBIND11_MODULE(lightsim2grid_cpp, m)
{
    // constant and compilation information
    m.attr("klu_solver_available") = py::bool_(this_KLU_SOLVER_AVAILABLE);
    m.attr("nicslu_solver_available") = py::bool_(this_NICSLU_SOLVER_AVAILABLE);
    m.attr("cktso_solver_available") = py::bool_(this_CKTSO_SOLVER_AVAILABLE);
    m.attr("compiled_march_native") = py::bool_(this__COMPILE_MARCHNATIVE);
    m.attr("compiled_o3_optim") = py::bool_(this__O3_OPTIM);
    m.attr("version") = py::str(this_VERSION);
    #ifdef NICSLU_PATH
    m.attr("nicslu_lib") = py::str(this_NICSLU_PATH);
    #endif
    #ifdef CKTSO_PATH
    m.attr("cktso_lib") = py::str(this_CKTSO_PATH);
    #endif

    bind_enums(m);
    bind_solvers(m);
    bind_containers(m);
    bind_misc(m);
    bind_gridmodel(m);
    bind_batch(m);
}
