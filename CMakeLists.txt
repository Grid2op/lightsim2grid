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

# --- Core install-dir overrides for Python wheel layout ---
# src/core/CMakeLists.txt defaults to a standard FHS layout (lib/, include/, ...).
# When building a wheel via scikit-build-core, everything must land inside the
# lightsim2grid/ package directory so it's importable after pip install.
if(SKBUILD)
    set(LS2G_INSTALL_LIBDIR     "${SKBUILD_PROJECT_NAME}"
        CACHE STRING "" FORCE)
    set(LS2G_INSTALL_INCLUDEDIR "${SKBUILD_PROJECT_NAME}/include"
        CACHE STRING "" FORCE)
    set(LS2G_INSTALL_CMAKEDIR   "${SKBUILD_PROJECT_NAME}/share/cmake/lightsim2grid_core"
        CACHE STRING "" FORCE)
endif()

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
# lightsim2grid_core + its headers + CMake config are installed by src/core/CMakeLists.txt.

# Python extension module
install(TARGETS lightsim2grid_cpp DESTINATION ${SKBUILD_PROJECT_NAME})

if(INSTALL_LEGACY_LIGHTSIM2GRID_CPP)
    install(TARGETS lightsim2grid_cpp DESTINATION .)
    # The root-level copy has rpath $ORIGIN/@loader_path (= site-packages root),
    # but liblightsim2grid_core lives in the lightsim2grid/ subdirectory.
    # Add a second rpath entry so auditwheel/delocate and the dynamic linker find it.
    if(APPLE)
        install(CODE "
            file(GLOB _legacy_sos
                \"\${CMAKE_INSTALL_PREFIX}/lightsim2grid_cpp*.so\")
            foreach(_so IN LISTS _legacy_sos)
                execute_process(
                    COMMAND install_name_tool
                        -add_rpath \"@loader_path/${SKBUILD_PROJECT_NAME}\"
                        \"\${_so}\"
                    RESULT_VARIABLE _ret)
                if(NOT _ret EQUAL 0)
                    message(WARNING \"install_name_tool failed for \${_so}\")
                endif()
            endforeach()
        ")
    elseif(UNIX)
        install(CODE "
            file(GLOB _legacy_sos
                \"\${CMAKE_INSTALL_PREFIX}/lightsim2grid_cpp*.so\")
            foreach(_so IN LISTS _legacy_sos)
                execute_process(
                    COMMAND patchelf --add-rpath \"\$ORIGIN/${SKBUILD_PROJECT_NAME}\"
                        \"\${_so}\"
                    RESULT_VARIABLE _ret)
                if(NOT _ret EQUAL 0)
                    message(WARNING \"patchelf failed for \${_so}\")
                endif()
            endforeach()
        ")
    endif()
endif()
