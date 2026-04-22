cmake_minimum_required(VERSION 3.15...4.99)

project(lightsim2grid
        VERSION "${SKBUILD_PROJECT_VERSION}"
        LANGUAGES CXX)

if("cxx_std_23" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(CMAKE_CXX_STANDARD 23)
elseif("cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(CMAKE_CXX_STANDARD 20)
elseif("cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(CMAKE_CXX_STANDARD 17)
elseif("cxx_std_14" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(CMAKE_CXX_STANDARD 14)
else()
    message(FATAL_ERROR "lightsim2grid requires at least C++14, but the compiler does not support it.")
endif()
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Dependencies ---
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# Eigen (header-only)
include_directories("eigen")

# --- Options that affect only the Python layer ---

# Documentation build
if(DEFINED ENV{_READ_THE_DOCS})
    option(BUILD_DOCS_ONLY "Flag for documentation generation" ON)
else()
    option(BUILD_DOCS_ONLY "Flag for documentation generation" OFF)
endif()

# Legacy install: also place lightsim2grid_cpp at site-packages root so that
# "from lightsim2grid_cpp import XXX" still works (deprecated).
option(INSTALL_LEGACY_LIGHTSIM2GRID_CPP
    "Also install lightsim2grid_cpp at site-packages root for legacy compatibility (deprecated)"
    ON)

# --- Build ---
# All solver detection (KLU/NICSLU/CKTSO), build flags (O3, march-native), and
# the lightsim2grid_core target are defined in src/core/CMakeLists.txt.
# src/bindings/python/CMakeLists.txt includes core if needed, then builds the
# pybind11 extension.
add_subdirectory("src/bindings/python")

# --- Python-only compile definitions (applied after targets exist) ---

# Read The Docs — only affects binding stubs
if(BUILD_DOCS_ONLY)
    target_compile_definitions(lightsim2grid_cpp PRIVATE _READ_THE_DOCS)
endif()

# Version macros — used in pickle_helpers.hpp (binding side only)
target_compile_definitions(lightsim2grid_cpp PRIVATE
    VERSION="${SKBUILD_PROJECT_VERSION_FULL}"
    VERSION_MAJOR="${CMAKE_PROJECT_VERSION_MAJOR}"
    VERSION_MEDIUM="${CMAKE_PROJECT_VERSION_MINOR}"
    VERSION_MINOR="${CMAKE_PROJECT_VERSION_PATCH}"
)

# --- Install rules ---

# Core shared library — inside the Python wheel (SKBUILD) or lib/ (standalone)
install(TARGETS lightsim2grid_core
    LIBRARY DESTINATION $<IF:$<BOOL:${SKBUILD}>,${SKBUILD_PROJECT_NAME},lib>
    ARCHIVE DESTINATION $<IF:$<BOOL:${SKBUILD}>,${SKBUILD_PROJECT_NAME},lib>
    RUNTIME DESTINATION $<IF:$<BOOL:${SKBUILD}>,${SKBUILD_PROJECT_NAME},bin>
)

# Python extension module
install(TARGETS lightsim2grid_cpp DESTINATION ${SKBUILD_PROJECT_NAME})

# Windows: install the import library (.lib) so plugin authors can link against it
if(WIN32)
    install(TARGETS lightsim2grid_cpp ARCHIVE DESTINATION ${SKBUILD_PROJECT_NAME})
endif()

if(INSTALL_LEGACY_LIGHTSIM2GRID_CPP)
    install(TARGETS lightsim2grid_cpp DESTINATION .)
endif()
