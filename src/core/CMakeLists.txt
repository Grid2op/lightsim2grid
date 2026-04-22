# ---- Standalone guard -------------------------------------------------------
# When included via add_subdirectory() from a parent project (e.g. the Python
# bindings CMakeLists or the root CMakeLists), the project() call is skipped.
# When invoked as a standalone "cmake src/core", the full project is set up so
# that pure-C++ users can build and install just the core library.
if(NOT PROJECT_NAME)
    cmake_minimum_required(VERSION 3.15...4.99)
    project(lightsim2grid_core VERSION 0.0.0 LANGUAGES CXX)

    if("cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
        set(CMAKE_CXX_STANDARD 17)
    elseif("cxx_std_14" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
        set(CMAKE_CXX_STANDARD 14)
    else()
        message(FATAL_ERROR "lightsim2grid_core requires at least C++14")
    endif()
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Eigen is two levels up (repo root / eigen/)
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../eigen")
endif()


# =============================================================================
# KLU (SuiteSparse) detection — 4 strategies, same logic as historical root
# =============================================================================
set(KLU_SOLVER_FOUND FALSE)

# Strategy 1: modern CMake config (e.g. libsuitesparse-dev on Debian/Ubuntu)
find_package(KLU QUIET)
if(KLU_FOUND)
    set(KLU_SOLVER_FOUND TRUE)
    message(STATUS "KLU: found via system CMake config")
    if(TARGET SuiteSparse::KLU_static)
        set(KLU_LINK_TARGET SuiteSparse::KLU_static)
    elseif(TARGET SuiteSparse::KLU)
        set(KLU_LINK_TARGET SuiteSparse::KLU)
    endif()
    find_package(CXSparse QUIET)
    if(TARGET SuiteSparse::CXSparse_static)
        set(CXSPARSE_LINK_TARGET SuiteSparse::CXSparse_static)
    elseif(TARGET SuiteSparse::CXSparse)
        set(CXSPARSE_LINK_TARGET SuiteSparse::CXSparse)
    else()
        find_library(CXSPARSE_FALLBACK_LIB NAMES cxsparse)
        find_path(CXSPARSE_FALLBACK_INCLUDE NAMES cs.h PATH_SUFFIXES suitesparse)
        if(CXSPARSE_FALLBACK_LIB AND CXSPARSE_FALLBACK_INCLUDE)
            add_library(CXSparse_fallback UNKNOWN IMPORTED)
            set_target_properties(CXSparse_fallback PROPERTIES
                IMPORTED_LOCATION "${CXSPARSE_FALLBACK_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES "${CXSPARSE_FALLBACK_INCLUDE}")
            set(CXSPARSE_LINK_TARGET CXSparse_fallback)
        endif()
    endif()
endif()

# Strategy 2: legacy find_library
if(NOT KLU_SOLVER_FOUND)
    find_library(KLU_SYSTEM_LIB NAMES klu)
    find_path(KLU_SYSTEM_INCLUDE NAMES klu.h PATH_SUFFIXES suitesparse)
    if(KLU_SYSTEM_LIB AND KLU_SYSTEM_INCLUDE)
        set(KLU_SOLVER_FOUND TRUE)
        message(STATUS "KLU: found via system libraries (${KLU_SYSTEM_LIB})")
        add_library(KLU_system UNKNOWN IMPORTED)
        set_target_properties(KLU_system PROPERTIES
            IMPORTED_LOCATION "${KLU_SYSTEM_LIB}"
            INTERFACE_INCLUDE_DIRECTORIES "${KLU_SYSTEM_INCLUDE}")
        set(KLU_LINK_TARGET KLU_system)
        find_library(CXSPARSE_SYSTEM_LIB NAMES cxsparse)
        find_path(CXSPARSE_SYSTEM_INCLUDE NAMES cs.h PATH_SUFFIXES suitesparse)
        if(CXSPARSE_SYSTEM_LIB AND CXSPARSE_SYSTEM_INCLUDE)
            add_library(CXSparse_system UNKNOWN IMPORTED)
            set_target_properties(CXSparse_system PROPERTIES
                IMPORTED_LOCATION "${CXSPARSE_SYSTEM_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES "${CXSPARSE_SYSTEM_INCLUDE}")
            set(CXSPARSE_LINK_TARGET CXSparse_system)
        endif()
    endif()
endif()

# Strategy 3: pre-built static libraries in the SuiteSparse submodule
if(NOT KLU_SOLVER_FOUND)
    set(_SS "${CMAKE_CURRENT_SOURCE_DIR}/../../SuiteSparse")
    set(_KLU_PREBUILT_LIB  "${_SS}/KLU/Lib/libklu.a")
    set(_KLU_PREBUILT_INC  "${_SS}/KLU/Include")
    if(EXISTS "${_KLU_PREBUILT_LIB}" AND EXISTS "${_KLU_PREBUILT_INC}/klu.h")
        message(STATUS "KLU: found pre-built static library in SuiteSparse submodule")
        set(KLU_SOLVER_FOUND TRUE)
        add_library(KLU_submodule STATIC IMPORTED)
        set_target_properties(KLU_submodule PROPERTIES
            IMPORTED_LOCATION "${_KLU_PREBUILT_LIB}"
            INTERFACE_INCLUDE_DIRECTORIES "${_KLU_PREBUILT_INC}")
        set(KLU_LINK_TARGET KLU_submodule)
        foreach(_dep AMD BTF COLAMD)
            string(TOLOWER "${_dep}" _dep_lc)
            set(_lib "${_SS}/${_dep}/Lib/lib${_dep_lc}.a")
            set(_inc "${_SS}/${_dep}/Include")
            if(EXISTS "${_lib}")
                add_library(${_dep}_submodule STATIC IMPORTED)
                set_target_properties(${_dep}_submodule PROPERTIES
                    IMPORTED_LOCATION "${_lib}"
                    INTERFACE_INCLUDE_DIRECTORIES "${_inc}")
                target_link_libraries(KLU_submodule INTERFACE ${_dep}_submodule)
            endif()
        endforeach()
        set(_ss_config_lib "${_SS}/SuiteSparse_config/libsuitesparseconfig.a")
        set(_ss_config_inc "${_SS}/SuiteSparse_config")
        if(EXISTS "${_ss_config_lib}")
            add_library(SuiteSparseConfig_submodule STATIC IMPORTED)
            set_target_properties(SuiteSparseConfig_submodule PROPERTIES
                IMPORTED_LOCATION "${_ss_config_lib}"
                INTERFACE_INCLUDE_DIRECTORIES "${_ss_config_inc}")
            target_link_libraries(KLU_submodule INTERFACE SuiteSparseConfig_submodule)
        endif()
        set(_cxsparse_lib "${_SS}/CXSparse/Lib/libcxsparse.a")
        set(_cxsparse_inc "${_SS}/CXSparse/Include")
        if(EXISTS "${_cxsparse_lib}")
            add_library(CXSparse_submodule STATIC IMPORTED)
            set_target_properties(CXSparse_submodule PROPERTIES
                IMPORTED_LOCATION "${_cxsparse_lib}"
                INTERFACE_INCLUDE_DIRECTORIES "${_cxsparse_inc}")
            set(CXSPARSE_LINK_TARGET CXSparse_submodule)
        endif()
    endif()
    unset(_SS)
    unset(_KLU_PREBUILT_LIB)
    unset(_KLU_PREBUILT_INC)
endif()

# Strategy 4: build KLU from the SuiteSparse git submodule via CMake
if(NOT KLU_SOLVER_FOUND)
    set(_SS_CMAKE "${CMAKE_CURRENT_SOURCE_DIR}/../../SuiteSparse/KLU/CMakeLists.txt")
    if(EXISTS "${_SS_CMAKE}")
        if(CMAKE_VERSION VERSION_LESS "3.22")
            message(WARNING
                "Building KLU from submodule requires CMake >= 3.22 "
                "(current: ${CMAKE_VERSION}). KLU solver will not be available.")
        else()
            message(STATUS "KLU: building static libraries from SuiteSparse submodule")
            set(BUILD_SHARED_LIBS        OFF CACHE BOOL "" FORCE)
            set(BUILD_STATIC_LIBS        ON  CACHE BOOL "" FORCE)
            set(SUITESPARSE_USE_OPENMP   OFF CACHE BOOL "" FORCE)
            set(SUITESPARSE_REQUIRE_BLAS OFF CACHE BOOL "" FORCE)
            set(SUITESPARSE_USE_CUDA     OFF CACHE BOOL "" FORCE)
            if(NOT BLA_VENDOR)
                set(BLA_VENDOR "Generic" CACHE STRING "" FORCE)
            endif()
            set(KLU_USE_CHOLMOD          OFF CACHE BOOL "" FORCE)
            set(SUITESPARSE_ROOT_CMAKELISTS ON)
            set(CMAKE_POSITION_INDEPENDENT_CODE ON)
            set(_SS "${CMAKE_CURRENT_SOURCE_DIR}/../../SuiteSparse")
            add_subdirectory(${_SS}/SuiteSparse_config ${CMAKE_CURRENT_BINARY_DIR}/SuiteSparse_config EXCLUDE_FROM_ALL)
            add_library(SuiteSparse::SuiteSparseConfig        ALIAS SuiteSparseConfig_static)
            add_library(SuiteSparse::SuiteSparseConfig_static ALIAS SuiteSparseConfig_static)
            add_subdirectory(${_SS}/AMD    ${CMAKE_CURRENT_BINARY_DIR}/AMD    EXCLUDE_FROM_ALL)
            add_library(SuiteSparse::AMD        ALIAS AMD_static)
            add_library(SuiteSparse::AMD_static ALIAS AMD_static)
            add_subdirectory(${_SS}/BTF    ${CMAKE_CURRENT_BINARY_DIR}/BTF    EXCLUDE_FROM_ALL)
            add_library(SuiteSparse::BTF        ALIAS BTF_static)
            add_library(SuiteSparse::BTF_static ALIAS BTF_static)
            add_subdirectory(${_SS}/COLAMD ${CMAKE_CURRENT_BINARY_DIR}/COLAMD EXCLUDE_FROM_ALL)
            add_library(SuiteSparse::COLAMD        ALIAS COLAMD_static)
            add_library(SuiteSparse::COLAMD_static ALIAS COLAMD_static)
            add_subdirectory(${_SS}/CXSparse ${CMAKE_CURRENT_BINARY_DIR}/CXSparse EXCLUDE_FROM_ALL)
            add_library(SuiteSparse::CXSparse        ALIAS CXSparse_static)
            add_library(SuiteSparse::CXSparse_static ALIAS CXSparse_static)
            add_subdirectory(${_SS}/KLU    ${CMAKE_CURRENT_BINARY_DIR}/KLU    EXCLUDE_FROM_ALL)
            set(KLU_SOLVER_FOUND     TRUE)
            set(KLU_LINK_TARGET      KLU_static)
            set(CXSPARSE_LINK_TARGET CXSparse_static)
            unset(_SS)
        endif()
    else()
        message(WARNING
            "KLU not found. Run: git submodule update --init SuiteSparse\n"
            "The KLU linear solver will not be available.")
    endif()
    unset(_SS_CMAKE)
endif()


# =============================================================================
# NICSLU
# =============================================================================
set(PATH_NICSLU "$ENV{PATH_NICSLU}" CACHE PATH "Path to NICSLU root directory")
option(ENABLE_NICSLU "Enable NICSLU solver" OFF)
if(PATH_NICSLU)
    if(WIN32)
        set(LIB_MAYBE_NICSLU ${PATH_NICSLU}/windows/lib_win7_x64_fma/int32)
    else()
        set(LIB_MAYBE_NICSLU ${PATH_NICSLU}/linux/lib_centos6_x64_gcc482_fma/int32)
    endif()
    find_library(LIB_NICSLU NAMES nicslu HINTS ${LIB_MAYBE_NICSLU})
    if(LIB_NICSLU)
        set(ENABLE_NICSLU ON)
        message(STATUS "NICSLU enabled via path: ${LIB_NICSLU}")
    else()
        message(WARNING "NICSLU disabled: no library found at ${LIB_MAYBE_NICSLU}")
    endif()
endif()


# =============================================================================
# CKTSO
# =============================================================================
set(PATH_CKTSO "$ENV{PATH_CKTSO}" CACHE PATH "Path to CKTSO root directory")
option(ENABLE_CKTSO "Enable CKTSO solver" OFF)
if(PATH_CKTSO)
    if(WIN32)
        set(LIB_MAYBE_CKTSO ${PATH_CKTSO}/windows_x64)
    else()
        set(LIB_MAYBE_CKTSO ${PATH_CKTSO}/ubuntu2004_x64_gcc940)
    endif()
    find_library(LIB_CKTSO NAMES cktso HINTS ${LIB_MAYBE_CKTSO})
    if(LIB_CKTSO)
        set(ENABLE_CKTSO ON)
        message(STATUS "CKTSO enabled via path: ${LIB_CKTSO}")
    else()
        message(WARNING "CKTSO disabled: no library found at ${LIB_MAYBE_CKTSO}")
    endif()
endif()


# =============================================================================
# Build flags
# =============================================================================
set(__COUT_TIMES "$ENV{__COUT_TIMES}" CACHE STRING "Print detailed execution timings")
if(__COUT_TIMES STREQUAL "1" OR __COUT_TIMES STREQUAL "True" OR __COUT_TIMES STREQUAL "true")
    option(ENABLE_COUT_TIMES "Enable execution time logging" ON)
else()
    option(ENABLE_COUT_TIMES "Enable execution time logging" OFF)
endif()

set(ENV_COMPILE_MARCHNATIVE "$ENV{__COMPILE_MARCHNATIVE}")
if(ENV_COMPILE_MARCHNATIVE STREQUAL "1" OR ENV_COMPILE_MARCHNATIVE STREQUAL "True")
    option(USE_MARCH_NATIVE "Enable -march=native optimization" ON)
else()
    option(USE_MARCH_NATIVE "Enable -march=native optimization" OFF)
endif()

set(ENV_O3_OPTIM "$ENV{__O3_OPTIM}")
if(ENV_O3_OPTIM STREQUAL "0" OR ENV_O3_OPTIM STREQUAL "False")
    option(USE_O3_OPTIM "Enable O3 optimization" OFF)
else()
    option(USE_O3_OPTIM "Enable O3 optimization" ON)
endif()


# =============================================================================
# Core library sources
# =============================================================================
set(SOURCES "")
add_subdirectory(linear_solvers)
add_subdirectory(batch_algorithm)
add_subdirectory(element_container)
add_subdirectory(powerflow_algorithm)

add_library(lightsim2grid_core SHARED
    "${CMAKE_CURRENT_SOURCE_DIR}/help_fun_msg.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BaseConstants.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GridModel.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ChooseSolver.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SolverRegistry.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SolverTypeNames.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BuiltinSolversRegistration.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SolverInstantiations.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/DataConverter.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SubstationContainer.cpp"
    ${SOURCES}
)

set_target_properties(lightsim2grid_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    EXPORT_NAME core
)

# Install-dir variables — can be overridden by a parent CMakeLists (e.g. for wheel builds)
set(LS2G_INSTALL_LIBDIR     "lib"
    CACHE STRING "Install destination for the shared library (relative to prefix)")
set(LS2G_INSTALL_INCLUDEDIR "include/lightsim2grid"
    CACHE STRING "Install destination for public headers (relative to prefix)")
set(LS2G_INSTALL_CMAKEDIR   "lib/cmake/lightsim2grid_core"
    CACHE STRING "Install destination for CMake config files (relative to prefix)")

target_include_directories(lightsim2grid_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${LS2G_INSTALL_INCLUDEDIR}>
)

target_compile_definitions(lightsim2grid_core PRIVATE LS2G_BUILDING_CORE)


# =============================================================================
# Link optional solvers
# KLU/CXSparse are STATIC libs whose symbols are fully embedded in the
# lightsim2grid_core SHARED library.  We use $<BUILD_INTERFACE:> so that:
#   - build-tree consumers (lightsim2grid_cpp bindings) get the include paths
#     and compile definitions they need to compile KLU-using headers;
#   - the EXPORTED CMake target does NOT carry KLU_submodule / CXSparse_submodule
#     in INTERFACE_LINK_LIBRARIES (those internal targets don't exist in an
#     installed or pre-built context, which would cause find_package consumers
#     to fail during their configure step).
# =============================================================================
if(KLU_SOLVER_FOUND)
    target_compile_definitions(lightsim2grid_core PUBLIC KLU_SOLVER_AVAILABLE)
    target_link_libraries(lightsim2grid_core PUBLIC $<BUILD_INTERFACE:${KLU_LINK_TARGET}>)
    if(DEFINED CXSPARSE_LINK_TARGET)
        target_link_libraries(lightsim2grid_core PUBLIC $<BUILD_INTERFACE:${CXSPARSE_LINK_TARGET}>)
    endif()
endif()

if(ENABLE_NICSLU)
    target_link_directories(lightsim2grid_core PRIVATE ${LIB_MAYBE_NICSLU})
    target_link_libraries(lightsim2grid_core PRIVATE nicslu)
    target_include_directories(lightsim2grid_core PRIVATE ${PATH_NICSLU}/include)
    target_compile_definitions(lightsim2grid_core PUBLIC  NICSLU_SOLVER_AVAILABLE)
    target_compile_definitions(lightsim2grid_core PRIVATE NICSLU_PATH="${LIB_MAYBE_NICSLU}")
endif()

if(LIB_CKTSO)
    target_link_directories(lightsim2grid_core PRIVATE ${LIB_MAYBE_CKTSO})
    target_link_libraries(lightsim2grid_core PRIVATE cktso)
    target_include_directories(lightsim2grid_core PRIVATE ${PATH_CKTSO}/include)
    target_compile_definitions(lightsim2grid_core PUBLIC  CKTSO_SOLVER_AVAILABLE)
    target_compile_definitions(lightsim2grid_core PRIVATE CKTSO_PATH="${LIB_MAYBE_CKTSO}")
endif()


# =============================================================================
# Compile-time flags
# =============================================================================
if(ENABLE_COUT_TIMES)
    target_compile_definitions(lightsim2grid_core PRIVATE __COUT_TIMES)
endif()

if(USE_O3_OPTIM)
    target_compile_definitions(lightsim2grid_core PRIVATE __O3_OPTIM)
endif()

if(NOT MSVC)
    if(USE_O3_OPTIM)
        target_compile_options(lightsim2grid_core PRIVATE -O3)
    endif()
    if(USE_MARCH_NATIVE)
        target_compile_options(lightsim2grid_core PRIVATE -march=native)
    endif()
else()
    if(USE_O3_OPTIM)
        target_compile_options(lightsim2grid_core PRIVATE /O2)
    endif()
endif()


# =============================================================================
# RPATH — both .so files find each other regardless of install location
# =============================================================================
if(APPLE)
    set(_LS2G_RPATH "@loader_path")
elseif(NOT WIN32)
    set(_LS2G_RPATH "$ORIGIN")
endif()
if(NOT WIN32)
    set_target_properties(lightsim2grid_core PROPERTIES
        INSTALL_RPATH "${_LS2G_RPATH}"
        BUILD_RPATH   "${_LS2G_RPATH}"
    )
endif()

if(WIN32)
    set_target_properties(lightsim2grid_core PROPERTIES ENABLE_EXPORTS ON)
endif()


# =============================================================================
# Install rules + CMake package config
# =============================================================================
include(CMakePackageConfigHelpers)

# Library + export set
install(TARGETS lightsim2grid_core
    EXPORT  lightsim2grid_coreTargets
    LIBRARY DESTINATION "${LS2G_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${LS2G_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "bin"
)

install(EXPORT lightsim2grid_coreTargets
    FILE      lightsim2grid_coreTargets.cmake
    NAMESPACE lightsim2grid::
    DESTINATION "${LS2G_INSTALL_CMAKEDIR}"
)

# Public headers (.hpp and .tpp template implementations, excluding cmake/)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    DESTINATION "${LS2G_INSTALL_INCLUDEDIR}"
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.tpp"
    PATTERN "cmake" EXCLUDE
)

# Bundle the Eigen header-only library alongside the ls2g headers.
# Eigen is part of the public interface (public headers #include "Eigen/Core"),
# so consumers need it at the same directory level as the ls2g headers.
# We install Eigen/ next to the ls2g headers so that the already-set
# INTERFACE_INCLUDE_DIRECTORIES (LS2G_INSTALL_INCLUDEDIR) covers both.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../eigen/Eigen")
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../eigen/Eigen"
        DESTINATION "${LS2G_INSTALL_INCLUDEDIR}"
    )
endif()

# Generate and install Config + ConfigVersion files
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/lightsim2grid_coreConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/lightsim2grid_coreConfig.cmake"
    INSTALL_DESTINATION "${LS2G_INSTALL_CMAKEDIR}"
    PATH_VARS LS2G_INSTALL_INCLUDEDIR
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/lightsim2grid_coreConfigVersion.cmake"
    VERSION       "${CMAKE_PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/lightsim2grid_coreConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/lightsim2grid_coreConfigVersion.cmake"
    DESTINATION "${LS2G_INSTALL_CMAKEDIR}"
)
