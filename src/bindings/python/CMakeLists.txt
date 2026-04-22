# =============================================================================
# lightsim2grid Python bindings
# =============================================================================
# This file builds the lightsim2grid_cpp pybind11 extension module.
#
# It follows the same pattern as the KLU detection in the root CMakeLists:
#   1. If lightsim2grid_core is already a known CMake target (e.g. when included
#      from the root project), use it directly.
#   2. If a pre-built lightsim2grid_core shared library can be found on disk,
#      import it as an IMPORTED target.
#   3. Otherwise, build lightsim2grid_core from src/core/ via add_subdirectory.
# =============================================================================

if(NOT TARGET lightsim2grid_core)
    find_library(LS2G_CORE_PREBUILT lightsim2grid_core
        HINTS
            "${CMAKE_CURRENT_BINARY_DIR}/../../core"    # sibling build tree
            "${CMAKE_BINARY_DIR}"                       # top-level build tree
            "${CMAKE_CURRENT_SOURCE_DIR}/../../../lightsim2grid"  # installed wheel
    )
    if(LS2G_CORE_PREBUILT)
        message(STATUS "lightsim2grid_core: using pre-built library (${LS2G_CORE_PREBUILT})")
        add_library(lightsim2grid_core SHARED IMPORTED GLOBAL)
        set_target_properties(lightsim2grid_core PROPERTIES
            IMPORTED_LOCATION "${LS2G_CORE_PREBUILT}"
            INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/../../core"
        )
    else()
        message(STATUS "lightsim2grid_core: not found pre-built, building from src/core")
        add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../core"
                         "${CMAKE_CURRENT_BINARY_DIR}/core")
    endif()
endif()

# Ensure Python and pybind11 are available (no-op when already found by parent)
if(NOT TARGET pybind11::pybind11)
    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
    find_package(pybind11 CONFIG REQUIRED)
endif()

pybind11_add_module(lightsim2grid_cpp
    "${CMAKE_CURRENT_SOURCE_DIR}/binding_module.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/binding_enums.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/binding_solvers.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/binding_containers.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/binding_misc.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/binding_gridmodel.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/binding_batch.cpp"
)

target_include_directories(lightsim2grid_cpp PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../core"
)

target_link_libraries(lightsim2grid_cpp PRIVATE lightsim2grid_core)

if(APPLE)
    set(_LS2G_RPATH "@loader_path")
elseif(NOT WIN32)
    set(_LS2G_RPATH "$ORIGIN")
endif()
if(NOT WIN32)
    set_target_properties(lightsim2grid_cpp PROPERTIES
        INSTALL_RPATH "${_LS2G_RPATH}"
        BUILD_RPATH   "${_LS2G_RPATH}"
    )
endif()
