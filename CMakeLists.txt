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

# KLU TODO
# --- Détection de KLU (SuiteSparse) ---
# set(SUITESPARSE_ENABLE_PROJECTS "klu;btf;amd;colamd;cxsparse")
# add_subdirectory("SuiteSparse")
# install(TARGETS klu btf amd colamd cxsparse DESTINATION ${SKBUILD_PROJECT_NAME})

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

# 3. Read The Docs
if(BUILD_DOCS_ONLY)
    target_compile_definitions(lightsim2grid_cpp PRIVATE _READ_THE_DOCS)
endif()

# 4. Cout Times
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

# 5. Optimisations
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
