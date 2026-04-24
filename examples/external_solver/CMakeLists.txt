cmake_minimum_required(VERSION 3.15)
project(dummy_solver CXX)
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
    if(NOT EXISTS "${LIGHTSIM2GRID_SRC}/SolverRegistry.hpp")
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
add_library(dummy_solver MODULE DummyExternalSolver.cpp)
target_link_libraries(dummy_solver PRIVATE lightsim2grid::core)

if(WIN32)
    # find_package provides the import lib via the IMPORTED target — nothing extra needed.
elseif(UNIX AND NOT APPLE)
    if(NOT lightsim2grid_core_FOUND)
        # Source-tree mode: no shared library to link against.
        # Symbols are resolved from the already-loaded lightsim2grid_cpp.so at runtime.
        target_link_options(dummy_solver PRIVATE -Wl,--allow-shlib-undefined)
    endif()
    set_target_properties(dummy_solver PROPERTIES PREFIX "lib" SUFFIX ".so")
elseif(APPLE)
    if(NOT lightsim2grid_core_FOUND)
        target_link_options(dummy_solver PRIVATE -undefined dynamic_lookup)
    endif()
    set_target_properties(dummy_solver PROPERTIES PREFIX "lib" SUFFIX ".so")
endif()
