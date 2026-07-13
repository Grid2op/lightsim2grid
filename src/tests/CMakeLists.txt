# ---- Standalone guard -------------------------------------------------------
# Mirrors src/core/CMakeLists.txt: when included via add_subdirectory() from
# the root CMakeLists (BUILD_TESTING=ON), the project() call is skipped. When
# invoked as a standalone "cmake -S src/tests", the full project is set up so
# the C++ unit tests build without python, pybind11 or the pip machinery --
# this is how the cpp_unit_tests CI workflow builds them.
if(NOT PROJECT_NAME)
    cmake_minimum_required(VERSION 3.16...4.99)
    project(lightsim2grid_tests LANGUAGES CXX)

    if("cxx_std_23" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
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

# The archive layer under test is deliberately compiled into the test binary
# instead of linking lightsim2grid_core: BinaryArchive.cpp is self-contained
# (see its header), so the tests need neither KLU/SuiteSparse detection nor a
# shared-library build, and compile in seconds. When solver-level unit tests
# arrive, link lightsim2grid_core here instead of growing this source list.
add_executable(lightsim2grid_unit_tests
    "${CMAKE_CURRENT_SOURCE_DIR}/../core/BinaryArchive.cpp"
    test_binary_archive.cpp
    test_binary_archive_stateres.cpp
    test_binary_archive_corruption.cpp
)
target_include_directories(lightsim2grid_unit_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../core")
# LS2G_API must resolve to dllexport (not dllimport) on MSVC, since the core
# translation unit is built into this binary rather than imported from a DLL
target_compile_definitions(lightsim2grid_unit_tests PRIVATE LS2G_BUILDING_CORE)
target_link_libraries(lightsim2grid_unit_tests PRIVATE Catch2::Catch2WithMain)

# one ctest entry per TEST_CASE
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../Catch2/extras")
include(Catch)
catch_discover_tests(lightsim2grid_unit_tests)
