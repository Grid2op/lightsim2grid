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

    # cxx_std_26 only appears in CMAKE_CXX_COMPILE_FEATURES when both CMake
    # (>= 3.30) and the compiler know C++26, so the check is safe with the
    # 3.16 minimum above: elsewhere it is simply never selected.
    if("cxx_std_26" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
        set(CMAKE_CXX_STANDARD 26)
    elseif("cxx_std_23" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
        set(CMAKE_CXX_STANDARD 23)
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
)
target_link_libraries(lightsim2grid_unit_tests PRIVATE lightsim2grid_core Catch2::Catch2WithMain)

# one ctest entry per TEST_CASE
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../Catch2/extras")
include(Catch)
catch_discover_tests(lightsim2grid_unit_tests)
