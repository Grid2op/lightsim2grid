# ---- Standalone guard -------------------------------------------------------
# Mirrors src/core/CMakeLists.txt: when included via add_subdirectory() from
# the root CMakeLists (BUILD_TESTING=ON), the project() call is skipped. When
# invoked as a standalone "cmake -S src/tests", the full project is set up so
# the C++ unit tests build without python, pybind11 or the pip machinery --
# this is how the cpp_unit_tests CI workflow builds them.
if(NOT PROJECT_NAME)
    cmake_minimum_required(VERSION 3.16...4.99)
    # the VERSION is required: src/core (added below) generates its CMake
    # package version file from PROJECT_VERSION
    project(lightsim2grid_tests VERSION 0.0.0 LANGUAGES CXX)

    # LS2G_CXX_STANDARD: force a specific C++ standard instead of auto-detecting
    # the newest one the compiler supports. Used by the cpp-standards CI
    # workflow to pin this standalone test build to C++14/C++26, e.g.:
    #   cmake -S src/tests -B build -DLS2G_CXX_STANDARD=14
    set(LS2G_CXX_STANDARD "" CACHE STRING "Force a specific C++ standard (14/17/20/23/26) instead of auto-detecting the newest supported")
    if(LS2G_CXX_STANDARD)
        if(NOT "cxx_std_${LS2G_CXX_STANDARD}" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
            message(FATAL_ERROR "LS2G_CXX_STANDARD=${LS2G_CXX_STANDARD} requested but the compiler does not support C++${LS2G_CXX_STANDARD}.")
        endif()
        set(CMAKE_CXX_STANDARD ${LS2G_CXX_STANDARD})
    # cxx_std_26 only appears in CMAKE_CXX_COMPILE_FEATURES when both CMake
    # (>= 3.30) and the compiler know C++26, so the check is safe with the
    # 3.16 minimum above: elsewhere it is simply never selected.
    elseif("cxx_std_26" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
        set(CMAKE_CXX_STANDARD 26)
    elseif("cxx_std_23" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
        set(CMAKE_CXX_STANDARD 23)
    elseif("cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
        set(CMAKE_CXX_STANDARD 20)
    elseif("cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
        set(CMAKE_CXX_STANDARD 17)
    elseif("cxx_std_14" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
        set(CMAKE_CXX_STANDARD 14)
    else()
        message(FATAL_ERROR "lightsim2grid requires at least C++14")
    endif()
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Eigen is two levels up (repo root / eigen/)
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../eigen")

    enable_testing()
endif()

# Catch2 (repo root / Catch2/ submodule). The guard keeps a future second
# consumer (or a parent that already provides Catch2) from adding it twice.
if(NOT TARGET Catch2::Catch2WithMain)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../Catch2" "${CMAKE_CURRENT_BINARY_DIR}/catch2")
endif()

# The full core library: needed since the tests exercise LSGrid (powerflow,
# solvers). Same guard as src/bindings/python/CMakeLists.txt: when added from
# the root CMakeLists the bindings have already created the target; standalone
# it is built from src/core (KLU detection included -- the core degrades
# gracefully to Eigen SparseLU when KLU is not found, and the LSGrid tests
# only rely on the always-available SparseLU algorithms).
if(NOT TARGET lightsim2grid_core)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../core" "${CMAKE_CURRENT_BINARY_DIR}/core")
endif()

add_executable(lightsim2grid_unit_tests
    test_binary_archive.cpp
    test_binary_archive_stateres.cpp
    test_binary_archive_corruption.cpp
    test_lsgrid.cpp
    test_case_exotic_elements.cpp
    test_ls2g_abi_tag.cpp
    test_plugin_registration.cpp
    test_check_grid.cpp
    test_state_poisoning.cpp
    test_powerflow_algorithm.cpp
)
target_link_libraries(lightsim2grid_unit_tests PRIVATE lightsim2grid_core Catch2::Catch2WithMain)

# CI hardening builds (see .github/workflows/sanitizers.yml). USE_SANITIZERS /
# USE_DEBUG_ASSERTS are the cached options src/core/CMakeLists.txt derives from
# the __SANITIZE / __DEBUG_ASSERTS environment variables; they are set there
# (added above) before this point, and simply undefined -- hence false -- when
# this file is consumed by a parent project that never created them.
#
# The core library carries the same flags already, but the test binary needs
# its OWN copy of them: ASan/UBSan only instrument the translation units they
# compile, so an out-of-bounds access made by test code (or by an Eigen
# expression instantiated in a test TU) is invisible unless the test objects
# are instrumented too -- and the executable must be LINKED with
# -fsanitize=... so the ASan runtime is the first DSO of the process (no
# LD_PRELOAD needed, unlike the python-extension jobs).
if(NOT MSVC)
    if(USE_SANITIZERS)
        target_compile_options(lightsim2grid_unit_tests PRIVATE
            -fsanitize=address,undefined -fno-omit-frame-pointer -g -O1)
        target_link_options(lightsim2grid_unit_tests PRIVATE -fsanitize=address,undefined)
    endif()
    if(USE_DEBUG_ASSERTS)
        # -UNDEBUG comes after the build type's -DNDEBUG on the command line,
        # so Eigen's bounds assertions come back to life in the test objects
        target_compile_options(lightsim2grid_unit_tests PRIVATE -UNDEBUG -D_GLIBCXX_ASSERTIONS)
    endif()
endif()

# one ctest entry per TEST_CASE
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../Catch2/extras")
include(Catch)
catch_discover_tests(lightsim2grid_unit_tests)
