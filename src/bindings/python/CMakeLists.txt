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
    "${CMAKE_CURRENT_SOURCE_DIR}/binding_lsgrid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/binding_batch.cpp"
)

target_include_directories(lightsim2grid_cpp PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../core"
)

target_link_libraries(lightsim2grid_cpp PRIVATE lightsim2grid_core)

# The env vars are read again here (rather than relying on the
# USE_MARCH_NATIVE/USE_O3_OPTIM/USE_SANITIZERS/USE_DEBUG_ASSERTS options
# defined in src/core/CMakeLists.txt) so the flags also apply when
# lightsim2grid_core is imported pre-built and src/core is never processed.
if(NOT MSVC)
    # -march=native changes which SIMD instruction sets Eigen sees enabled
    # (__AVX__, __AVX2__, ...), which changes EIGEN_MAX_ALIGN_BYTES and thus
    # how Eigen aligns/allocates/frees its dynamic-size matrices. lightsim2grid_cpp
    # instantiates the same Eigen::Matrix<...> types as lightsim2grid_core (e.g.
    # via pybind11's Eigen return-value caster) whenever a function's return
    # value crosses the module boundary, so if the two are compiled with
    # different -march flags, an Eigen object allocated in one and freed in the
    # other computes a different alignment offset -- heap corruption
    # ("free(): invalid size" / "double free or corruption") at the first such
    # crossing. Must stay identical to src/core/CMakeLists.txt's USE_MARCH_NATIVE.
    # target_compile_definitions (not just _options) below: binding_module.cpp's
    # compiled_march_native/compiled_o3_optim (reported via lightsim2grid.
    # compilation_options) check the __COMPILE_MARCHNATIVE/__O3_OPTIM *macros*,
    # separately from whether the compiler flags are actually active.
    set(ENV_COMPILE_MARCHNATIVE "$ENV{__COMPILE_MARCHNATIVE}")
    if(ENV_COMPILE_MARCHNATIVE STREQUAL "1" OR ENV_COMPILE_MARCHNATIVE STREQUAL "True")
        target_compile_options(lightsim2grid_cpp PRIVATE -march=native)
        target_compile_definitions(lightsim2grid_cpp PRIVATE __COMPILE_MARCHNATIVE)
    endif()

    set(ENV_O3_OPTIM "$ENV{__O3_OPTIM}")
    if(ENV_O3_OPTIM STREQUAL "1" OR ENV_O3_OPTIM STREQUAL "True")
        target_compile_options(lightsim2grid_cpp PRIVATE -O3)
        target_compile_definitions(lightsim2grid_cpp PRIVATE __O3_OPTIM)
    endif()

    set(ENV_SANITIZE "$ENV{__SANITIZE}")
    if(ENV_SANITIZE STREQUAL "1" OR ENV_SANITIZE STREQUAL "True")
        target_compile_options(lightsim2grid_cpp PRIVATE
            -fsanitize=address,undefined -fno-omit-frame-pointer -g -O1)
        target_link_options(lightsim2grid_cpp PRIVATE -fsanitize=address,undefined)
    endif()
    set(ENV_DEBUG_ASSERTS "$ENV{__DEBUG_ASSERTS}")
    if(ENV_DEBUG_ASSERTS STREQUAL "1" OR ENV_DEBUG_ASSERTS STREQUAL "True")
        target_compile_options(lightsim2grid_cpp PRIVATE -UNDEBUG -D_GLIBCXX_ASSERTIONS)
    endif()
endif()

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
