cmake_minimum_required(VERSION 3.15)
project(linear_solver_plugin C CXX)
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

# -----------------------------------------------------------------------
# The linear solver library. Here, the stand-in of the mylu/ folder, built with the plugin
# so that the example is self-contained. With a real library, find it instead, eg:
#   find_library(MYLU_LIB NAMES mylu HINTS "${MYLU_ROOT}/lib" REQUIRED)
#   find_path(MYLU_INCLUDE mylu.h HINTS "${MYLU_ROOT}/include" REQUIRED)
# and link ${MYLU_LIB} / include ${MYLU_INCLUDE} below. A static library (or one built with
# position-independent code) ends up inside the plugin; a shared one must be found by the
# dynamic loader when the plugin is loaded.
# -----------------------------------------------------------------------
add_library(mylu STATIC mylu/mylu.c)
set_target_properties(mylu PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(mylu PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/mylu")
if(UNIX)
    target_link_libraries(mylu PRIVATE m)
endif()

# The plugin is a MODULE (loadable at runtime via dlopen/LoadLibrary).
add_library(linear_solver_plugin MODULE MyLU_plugin.cpp)
target_link_libraries(linear_solver_plugin PRIVATE lightsim2grid::core mylu)

# Match lightsim2grid_core's -march=native / -O3 -- see docs/solver_plugin.rst
# and examples/cmake/MatchLightsim2gridBuildFlags.cmake for why this matters
# (skipping it is a real, silent heap-corruption risk, not just a perf loss).
include("${CMAKE_CURRENT_SOURCE_DIR}/../cmake/MatchLightsim2gridBuildFlags.cmake")
ls2g_match_core_build_flags(linear_solver_plugin)

if(WIN32)
    # find_package provides the import lib via the IMPORTED target — nothing extra needed.
elseif(UNIX AND NOT APPLE)
    if(NOT lightsim2grid_core_FOUND)
        # Source-tree mode: no shared library to link against.
        # Symbols are resolved from the already-loaded lightsim2grid_cpp.so at runtime.
        target_link_options(linear_solver_plugin PRIVATE -Wl,--allow-shlib-undefined)
    endif()
    set_target_properties(linear_solver_plugin PROPERTIES PREFIX "lib" SUFFIX ".so")
elseif(APPLE)
    if(NOT lightsim2grid_core_FOUND)
        target_link_options(linear_solver_plugin PRIVATE -undefined dynamic_lookup)
    endif()
    set_target_properties(linear_solver_plugin PROPERTIES PREFIX "lib" SUFFIX ".so")
endif()
