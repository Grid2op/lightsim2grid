// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "Ls2gAbiTag.hpp"

#include <stdexcept>

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
    // Guard against the exact bug that motivated this check: lightsim2grid_core
    // and lightsim2grid_cpp are two separate shared libraries, and if they were
    // compiled with different -march (or otherwise disagree on Eigen's SIMD/
    // alignment settings), an Eigen object allocated in one and freed in the
    // other silently corrupts the heap (see docs/solver_plugin.rst). Comparing
    // the two ABI tags here -- bindings' own view vs. core's compiled-in view,
    // queried via a real cross-.so function call -- catches the mismatch at
    // import time instead. PYBIND11_MODULE wraps this function body in a
    // try/catch, so throwing here surfaces as a clean Python ImportError.
    {
        const ls2g::Ls2gAbiTag bindings_tag = ls2g::ls2g_current_abi_tag();
        const ls2g::Ls2gAbiTag core_tag = ls2g::core_abi_tag();
        if (bindings_tag != core_tag) {
            throw std::runtime_error(
                "lightsim2grid: lightsim2grid_cpp (the Python bindings) and lightsim2grid_core were "
                "compiled with different Eigen SIMD/alignment settings -- loading them together would "
                "silently corrupt the heap the first time an Eigen object crosses the module boundary. "
                "lightsim2grid_core was built with [" + core_tag.describe() + "], lightsim2grid_cpp was "
                "built with [" + bindings_tag.describe() + "]. This should not happen with the provided "
                "build system; if you customized the build, make sure both are compiled with identical "
                "-march flags -- see docs/solver_plugin.rst, section \"Matching build flags\".");
        }
    }

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
