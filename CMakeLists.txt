cmake_minimum_required(VERSION 3.15...4.0)

project(lightsim2grid
        VERSION "${SKBUILD_PROJECT_VERSION}"
        LANGUAGES CXX)

# --- Dépendances ---
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)


# --- Handle dependencies ---
# Eigen
# add_subdirectory("eigen")
include_directories("eigen")

# --- Détection de KLU (SuiteSparse) ---
# Strategy 1: modern CMake config (e.g. from libsuitesparse-dev on Debian/Ubuntu)
set(KLU_SOLVER_FOUND FALSE)
find_package(KLU QUIET)
if(KLU_FOUND)
    set(KLU_SOLVER_FOUND TRUE)
    message(STATUS "KLU: found via system CMake config")
    if(TARGET SuiteSparse::KLU_static)
        set(KLU_LINK_TARGET SuiteSparse::KLU_static)
    elseif(TARGET SuiteSparse::KLU)
        set(KLU_LINK_TARGET SuiteSparse::KLU)
    endif()
    # CXSparse provides cs.h used by KLUSolver.hpp
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

# Strategy 2: legacy find_library (older SuiteSparse without CMake config files)
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

# Strategy 3: build KLU from the SuiteSparse git submodule
if(NOT KLU_SOLVER_FOUND)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/KLU/CMakeLists.txt")
        if(CMAKE_VERSION VERSION_LESS "3.22")
            message(FATAL_ERROR
                "Building KLU from the SuiteSparse submodule requires CMake >= 3.22. "
                "Either upgrade CMake or install SuiteSparse system-wide "
                "(e.g. apt install libsuitesparse-dev).")
        endif()
        message(STATUS "KLU: building static libraries from SuiteSparse submodule")
        # Force static-only build for all SuiteSparse components
        set(BUILD_SHARED_LIBS        OFF CACHE BOOL "" FORCE)
        set(BUILD_STATIC_LIBS        ON  CACHE BOOL "" FORCE)
        # Tell each component it is being built from a parent CMakeLists so it
        # skips its internal find_package() calls for sibling dependencies.
        set(SUITESPARSE_ROOT_CMAKELISTS ON)
        # Disable the optional CHOLMOD dependency inside KLU.
        set(KLU_USE_CHOLMOD          OFF CACHE BOOL "" FORCE)
        # Disable OpenMP: SuiteSparse_config's find_package(OpenMP) creates a
        # directory-local target on CMake < 3.24, which is invisible to sibling
        # subdirectories (AMD, BTF, …) and causes generator-expression errors.
        # KLU does not need OpenMP, so this has no functional impact.
        set(SUITESPARSE_USE_OPENMP   OFF CACHE BOOL "" FORCE)
        # Disable mandatory BLAS: SuiteSparseBLAS.cmake ends with
        # find_package(BLAS REQUIRED) which fails on systems without a BLAS.
        # KLU/AMD/BTF/COLAMD/CXSparse do not use BLAS; it is only needed by
        # CHOLMOD/UMFPACK/SPQR which we are not building.
        set(SUITESPARSE_REQUIRE_BLAS OFF CACHE BOOL "" FORCE)

        # Build components in dependency order and expose SuiteSparse:: aliases
        # so that later components can reference them via target_link_libraries.
        add_subdirectory(SuiteSparse/SuiteSparse_config EXCLUDE_FROM_ALL)
        add_library(SuiteSparse::SuiteSparseConfig        ALIAS SuiteSparseConfig_static)
        add_library(SuiteSparse::SuiteSparseConfig_static ALIAS SuiteSparseConfig_static)

        add_subdirectory(SuiteSparse/AMD EXCLUDE_FROM_ALL)
        add_library(SuiteSparse::AMD        ALIAS AMD_static)
        add_library(SuiteSparse::AMD_static ALIAS AMD_static)

        add_subdirectory(SuiteSparse/BTF EXCLUDE_FROM_ALL)
        add_library(SuiteSparse::BTF        ALIAS BTF_static)
        add_library(SuiteSparse::BTF_static ALIAS BTF_static)

        add_subdirectory(SuiteSparse/COLAMD EXCLUDE_FROM_ALL)
        add_library(SuiteSparse::COLAMD        ALIAS COLAMD_static)
        add_library(SuiteSparse::COLAMD_static ALIAS COLAMD_static)

        add_subdirectory(SuiteSparse/CXSparse EXCLUDE_FROM_ALL)
        add_library(SuiteSparse::CXSparse        ALIAS CXSparse_static)
        add_library(SuiteSparse::CXSparse_static ALIAS CXSparse_static)

        add_subdirectory(SuiteSparse/KLU EXCLUDE_FROM_ALL)

        set(KLU_SOLVER_FOUND    TRUE)
        set(KLU_LINK_TARGET     KLU_static)
        set(CXSPARSE_LINK_TARGET CXSparse_static)
    else()
        message(WARNING
            "KLU not found system-wide and the SuiteSparse submodule is not initialised. "
            "Run: git submodule update --init SuiteSparse\n"
            "The KLU solver will not be available.")
    endif()
endif()

# --- Gestion des Options & Variables d'Environnement ---

# 1. NICSLU
set(PATH_NICSLU "$ENV{PATH_NICSLU}" CACHE PATH "Path to NICSLU root directory")
option(ENABLE_NICSLU "Enable NICSLU solver" OFF)

if(PATH_NICSLU)
    if(WIN32)
        set(LIB_MAYBE_NICSLU ${PATH_NICSLU}/windows/lib_win7_x64_fma/int32)
    else()
        set(LIB_MAYBE_NICSLU ${PATH_NICSLU}/linux/lib_centos6_x64_gcc482_fma/int32)
    endif()

    find_library(LIB_NICSLU 
                 NAMES nicslu
                 HINTS ${LIB_MAYBE_NICSLU})
    if(LIB_NICSLU)
        set(ENABLE_NICSLU ON)
        message(STATUS "NICSLU enabled via path: ${LIB_NICSLU}")
    else()
        message(WARNING "NICSLU disable no library found at ${LIB_MAYBE_NICSLU}")
    endif()
endif()

# 2. CKTSO
set(PATH_CKTSO "$ENV{PATH_CKTSO}" CACHE PATH "Path to CKTSO root directory")
option(ENABLE_CKTSO "Enable CKTSO solver" OFF)

if(PATH_CKTSO)
    if(WIN32)
        set(LIB_MAYBE_CKTSO ${PATH_CKTSO}/windows_x64)
    else()
        set(LIB_MAYBE_CKTSO ${PATH_CKTSO}/ubuntu2004_x64_gcc940)
    endif()

    find_library(LIB_CKTSO 
                 NAMES cktso
                 HINTS ${LIB_MAYBE_CKTSO})
    if(LIB_CKTSO)
        # add_library(libcktso SHARED IMPORTED)
        set(ENABLE_CKTSO ON)
        message(STATUS "CKTSO enabled via path: ${LIB_CKTSO}")
    else()
        message(WARNING "CKTSO disable no library found at ${LIB_MAYBE_CKTSO}")
    endif()
endif()

# 3. Flags divers
# Debug / Profiling time
set(__COUT_TIMES "$ENV{__COUT_TIMES}" CACHE PATH "Print detail execution timings in std cout")
if(__COUT_TIMES STREQUAL "1" OR __COUT_TIMES STREQUAL "True" OR __COUT_TIMES STREQUAL "true") 
    option(ENABLE_COUT_TIMES "Enable execution time logging" ON)
else()
    option(ENABLE_COUT_TIMES "Enable execution time logging" OFF)
endif()

# Documentation build
if(DEFINED ENV{_READ_THE_DOCS})
    option(BUILD_DOCS_ONLY "Flag for documentation generation" ON)
else()
    option(BUILD_DOCS_ONLY "Flag for documentation generation" OFF)
endif()

# Optimisation March Native
if(DEFINED ENV{__COMPILE_MARCHNATIVE})
    option(USE_MARCH_NATIVE "Enable -march=native optimization" ON)
else()
    option(USE_MARCH_NATIVE "Enable -march=native optimization" OFF)
endif()

# Optimisation O3 (activated by default)
set(ENV_O3_OPTIM "$ENV{__O3_OPTIM}")
if(ENV_O3_OPTIM STREQUAL "0" OR ENV_O3_OPTIM STREQUAL "False")
    option(USE_O3_OPTIM "Enable O3 optimization" OFF)
else()
    option(USE_O3_OPTIM "Enable O3 optimization" ON)
endif()

# --- Compile lightsim2grid sources ---
add_subdirectory("src")

# --- Build Options ---

# 1. Configuration NICSLU
if(ENABLE_NICSLU)
    target_link_directories(lightsim2grid_cpp PRIVATE ${LIB_MAYBE_NICSLU})
    target_link_libraries(lightsim2grid_cpp PRIVATE nicslu)
    target_include_directories(lightsim2grid_cpp PRIVATE ${PATH_NICSLU}/include)
    target_compile_definitions(lightsim2grid_cpp PRIVATE NICSLU_PATH="${LIB_MAYBE_NICSLU}")
    target_compile_definitions(lightsim2grid_cpp PRIVATE NICSLU_SOLVER_AVAILABLE)
endif()

# 2. Configuration CKTSO
if(LIB_CKTSO)
    target_link_directories(lightsim2grid_cpp PRIVATE ${LIB_MAYBE_CKTSO})
    target_link_libraries(lightsim2grid_cpp PRIVATE cktso)
    target_include_directories(lightsim2grid_cpp PRIVATE ${PATH_CKTSO}/include)
    target_compile_definitions(lightsim2grid_cpp PRIVATE CKTSO_PATH="${LIB_MAYBE_CKTSO}")
    target_compile_definitions(lightsim2grid_cpp PRIVATE CKTSO_SOLVER_AVAILABLE)
endif()

# 3. Configuration KLU (SuiteSparse)
if(KLU_SOLVER_FOUND)
    target_compile_definitions(lightsim2grid_cpp PRIVATE KLU_SOLVER_AVAILABLE)
    target_link_libraries(lightsim2grid_cpp PRIVATE ${KLU_LINK_TARGET})
    if(DEFINED CXSPARSE_LINK_TARGET)
        target_link_libraries(lightsim2grid_cpp PRIVATE ${CXSPARSE_LINK_TARGET})
    endif()
endif()

# 5. Read The Docs
if(BUILD_DOCS_ONLY)
    target_compile_definitions(lightsim2grid_cpp PRIVATE _READ_THE_DOCS)
endif()

# 6. Cout Times
if(ENABLE_COUT_TIMES)
    target_compile_definitions(lightsim2grid_cpp PRIVATE __COUT_TIMES)
endif()

# version
target_compile_definitions(lightsim2grid_cpp PRIVATE 
    VERSION="${SKBUILD_PROJECT_VERSION_FULL}"
)
target_compile_definitions(lightsim2grid_cpp PRIVATE 
    VERSION_MAJOR="${CMAKE_PROJECT_VERSION_MAJOR}"
)
target_compile_definitions(lightsim2grid_cpp PRIVATE 
    VERSION_MEDIUM="${CMAKE_PROJECT_VERSION_MINOR}"
)
target_compile_definitions(lightsim2grid_cpp PRIVATE 
    VERSION_MINOR="${CMAKE_PROJECT_VERSION_PATCH}"
)

# --- Compilation options

# 7. Optimisations
if(USE_O3_OPTIM)
    target_compile_definitions(lightsim2grid_cpp PRIVATE 
        __O3_OPTIM
    )
endif()

if(NOT MSVC)
    if(USE_O3_OPTIM)
        target_compile_options(lightsim2grid_cpp PRIVATE -O3)
    endif()
    
    if(USE_MARCH_NATIVE)
        target_compile_options(lightsim2grid_cpp PRIVATE -march=native)
    endif()
else()
    if(USE_O3_OPTIM)
        target_compile_options(lightsim2grid_cpp PRIVATE /O2)
    endif()
endif()

# tells cmake to install the cpp extension within the python project
# so that we can do things (from python) like
# from lightsim2grid.lightsim2grid_cpp import XXX
# instead of from lightsim2grid_cpp import XXX
install(TARGETS lightsim2grid_cpp DESTINATION ${SKBUILD_PROJECT_NAME})
