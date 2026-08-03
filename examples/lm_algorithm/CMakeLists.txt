cmake_minimum_required(VERSION 3.15)
project(lm_algorithm CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -----------------------------------------------------------------------
# Strategy 1: find via installed CMake config
#   pip install lightsim2grid  →  lightsim2grid.get_cmake_dir()
#   cmake --install             →  _install/lib/cmake/lightsim2grid_core/
# Pass -DLIGHTSIM2GRID_CMAKE_DIR=<path> to hint the location, or rely on
# CMAKE_PREFIX_PATH.
# -----------------------------------------------------------------------
find_package(lightsim2grid_core CONFIG QUIET
    HINTS "${LIGHTSIM2GRID_CMAKE_DIR}")

# -----------------------------------------------------------------------
# Strategy 2: fall back to source tree (in-repo development builds).
# -----------------------------------------------------------------------
if(NOT lightsim2grid_core_FOUND)
    if(NOT DEFINED LIGHTSIM2GRID_SRC)
        set(LIGHTSIM2GRID_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../../src/core")
    endif()
    if(NOT DEFINED Eigen3_INCLUDE)
        set(Eigen3_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../../eigen")
    endif()
    if(NOT EXISTS "${LIGHTSIM2GRID_SRC}/AlgorithmRegistry.hpp")
        message(FATAL_ERROR
            "lightsim2grid_core not found.\n"
            "Either install lightsim2grid and pass:\n"
            "  -DLIGHTSIM2GRID_CMAKE_DIR=<cmake-dir>\n"
            "(obtain the path via: python -c \"import lightsim2grid; print(lightsim2grid.get_cmake_dir())\")\n"
            "or pass -DLIGHTSIM2GRID_SRC=<path/to/src/core> for a source-tree build.")
    endif()
    add_library(lightsim2grid_core_iface INTERFACE)
    target_include_directories(lightsim2grid_core_iface INTERFACE
        "${LIGHTSIM2GRID_SRC}"
        "${Eigen3_INCLUDE}"
    )
    add_library(lightsim2grid::core ALIAS lightsim2grid_core_iface)
    message(STATUS "lightsim2grid_core: using source tree at ${LIGHTSIM2GRID_SRC}")
endif()

# The plugin is a MODULE (loadable at runtime via dlopen/LoadLibrary).
add_library(lm_algorithm MODULE LMNRAlgo_plugin.cpp)
target_link_libraries(lm_algorithm PRIVATE lightsim2grid::core)

# Match lightsim2grid_core's -march=native / -O3 -- see docs/solver_plugin.rst
# and examples/cmake/MatchLightsim2gridBuildFlags.cmake for why this matters
# (skipping it is a real, silent heap-corruption risk, not just a perf loss).
include("${CMAKE_CURRENT_SOURCE_DIR}/../cmake/MatchLightsim2gridBuildFlags.cmake")
ls2g_match_core_build_flags(lm_algorithm)

# -----------------------------------------------------------------------
# KLUSolver.hpp (pulled in transitively via Solvers.hpp) needs SuiteSparse's
# own headers (cs.h, klu.h, amd.h, colamd.h, btf.h, SuiteSparse_config.h) to
# parse the KLULinearSolver class declaration. Unlike
# examples/dist_slack_algorithm (which reuses an existing `extern template`
# NRAlgo<..., KLULinearSolver, ...> instantiation and links straight against
# it), this plugin instantiates its own fresh
# LMNRAlgo<LinearSolverPolicy<KLULinearSolver>, ...> / NRSystem<...,
# DiagonalEntry> from headers -- but LinearSolverPolicy<KLULinearSolver> is a
# thin counting/timing wrapper: its methods just forward into
# KLULinearSolver's own analyze/factorize/refactorize/solve, which remain
# ordinary, already-compiled-into-core symbols. So we still only need KLU's
# *headers* here, never its implementation. Neither an installed
# lightsim2grid_core CMake package nor a plain source-tree checkout exposes
# those on its public include path (a private implementation detail of the
# core .so), so point at the SuiteSparse submodule directly, header-only, as
# a repo-relative fallback.
# -----------------------------------------------------------------------
if(NOT DEFINED SUITESPARSE_SRC)
    set(SUITESPARSE_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../../SuiteSparse")
endif()
if(EXISTS "${SUITESPARSE_SRC}/KLU/Include/klu.h")
    target_include_directories(lm_algorithm PRIVATE
        "${SUITESPARSE_SRC}/KLU/Include"
        "${SUITESPARSE_SRC}/AMD/Include"
        "${SUITESPARSE_SRC}/COLAMD/Include"
        "${SUITESPARSE_SRC}/BTF/Include"
        "${SUITESPARSE_SRC}/SuiteSparse_config"
        "${SUITESPARSE_SRC}/CXSparse/Include"
    )
endif()

if(WIN32)
    # find_package provides the import lib via the IMPORTED target — nothing extra needed.
elseif(UNIX AND NOT APPLE)
    if(NOT lightsim2grid_core_FOUND)
        # Source-tree mode: no shared library to link against.
        # Symbols are resolved from the already-loaded lightsim2grid_cpp.so at runtime.
        target_link_options(lm_algorithm PRIVATE -Wl,--allow-shlib-undefined)
    endif()
    set_target_properties(lm_algorithm PROPERTIES PREFIX "lib" SUFFIX ".so")
elseif(APPLE)
    if(NOT lightsim2grid_core_FOUND)
        target_link_options(lm_algorithm PRIVATE -undefined dynamic_lookup)
    endif()
    set_target_properties(lm_algorithm PROPERTIES PREFIX "lib" SUFFIX ".so")
endif()
