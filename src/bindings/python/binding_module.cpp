// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "Ls2gAbiTag.hpp"
#include "AlgorithmRegistry.hpp"

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>

#if defined(_WIN32)
#  include <windows.h>
#else
#  include <dlfcn.h>
#endif

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

// Behind lightsim2grid.load_algorithm_plugin(). Load the shared library, look
// up its exported ls2g_register_plugin entry point, and call it -- all from a
// real C++ frame that pybind11 wraps in a try/catch, so any registration
// failure (ABI mismatch, duplicate solver name, missing entry point, missing
// file) surfaces as a normal Python exception. This is the whole safety
// improvement over dlopen'ing the plugin directly from ctypes: there, an
// exception thrown by a static-constructor registration would unwind into the
// C loader frames and std::terminate() the interpreter. A rejected plugin is
// unloaded again so it leaves no trace; a successfully-registered one is
// intentionally kept mapped for the life of the process, since the registry now
// holds factories whose code lives in it.
// Size of the stack diagnostic buffer handed to a plugin's entry point. The
// entry point is told this size and its contract is to write at most this many
// bytes, NUL included; plugins built with LS2G_PLUGIN_ENTRY honour it by
// construction (the bounded write in ls2g::detail::ls2g_run_plugin_entry).
static const std::size_t max_errbuf_size = 512;

static void load_algorithm_plugin_impl(const std::string& path)
{
    // Platform-specific load + symbol lookup. Both branches leave `entry` set
    // to the resolved entry point (or null) and `unload` able to undo the load,
    // so the entry call and all buffer handling below happen exactly once.
    ls2g::ls2g_plugin_entry_fn entry = nullptr;
    std::function<void()> unload;
#if defined(_WIN32)
    HMODULE handle = ::LoadLibraryA(path.c_str());
    if (handle == nullptr) {
        throw std::runtime_error("lightsim2grid.load_algorithm_plugin: could not load '" + path +
                                 "' (LoadLibrary failed, error code " + std::to_string(::GetLastError()) + ").");
    }
    unload = [handle]{ ::FreeLibrary(handle); };
    entry = reinterpret_cast<ls2g::ls2g_plugin_entry_fn>(
        ::GetProcAddress(handle, "ls2g_register_plugin"));
#else
    ::dlerror();  // clear any stale error state
    void* handle = ::dlopen(path.c_str(), RTLD_NOW | RTLD_GLOBAL);
    if (handle == nullptr) {
        const char* e = ::dlerror();
        throw std::runtime_error("lightsim2grid.load_algorithm_plugin: could not load '" + path +
                                 "': " + (e != nullptr ? e : "unknown dlopen error"));
    }
    unload = [handle]{ ::dlclose(handle); };
    ::dlerror();
    entry = reinterpret_cast<ls2g::ls2g_plugin_entry_fn>(
        ::dlsym(handle, "ls2g_register_plugin"));
    if (::dlerror() != nullptr) entry = nullptr;
#endif
    if (entry == nullptr) {
        unload();
        throw std::runtime_error("lightsim2grid.load_algorithm_plugin: '" + path +
                                 "' is not a lightsim2grid plugin (no 'ls2g_register_plugin' entry point).");
    }

    // errbuf is written only by `entry`, which is told its full size and must
    // not exceed it. Force a NUL into the last slot afterwards regardless: even
    // a misbehaving plugin that ignored the size contract then cannot make the
    // C-string read below run past the buffer -- the host side cannot over-read
    // by construction, independently of what the plugin wrote.
    char errbuf[max_errbuf_size];
    errbuf[0] = '\0';
    const int rc = entry(&ls2g::AlgorithmRegistry::instance(), errbuf, max_errbuf_size);
    errbuf[max_errbuf_size - 1] = '\0';

    if (rc != 0) {
        unload();  // reject: unload so a failed plugin leaves nothing behind
        throw std::runtime_error("lightsim2grid.load_algorithm_plugin: '" + path +
                                 "' failed to register: " + errbuf);
    }
    // success: keep the library mapped for the life of the process, since the
    // registry now holds factories whose code lives in it.
}

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

    m.def("_load_algorithm_plugin", &load_algorithm_plugin_impl, py::arg("path"),
          "Internal helper behind lightsim2grid.load_algorithm_plugin(): load the given "
          "shared library, run its 'ls2g_register_plugin' entry point, and raise a Python "
          "exception (rather than crashing the interpreter) if registration fails.");

    bind_enums(m);
    bind_solvers(m);
    bind_containers(m);
    bind_misc(m);
    bind_gridmodel(m);
    bind_batch(m);
}
