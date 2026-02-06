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


# --- Gestion des Options & Variables d'Environnement ---

# 1. NICSLU
set(PATH_NICSLU "$ENV{PATH_NICSLU}" CACHE PATH "Path to NICSLU root directory")
option(ENABLE_NICSLU "Enable NICSLU solver" OFF)

if(PATH_NICSLU)
    set(ENABLE_NICSLU ON)
    message(STATUS "NICSLU enabled via path: ${PATH_NICSLU}")
endif()

# 2. CKTSO
set(PATH_CKTSO "$ENV{PATH_CKTSO}" CACHE PATH "Path to CKTSO root directory")
option(ENABLE_CKTSO "Enable CKTSO solver" OFF)

if(PATH_CKTSO)
    set(ENABLE_CKTSO ON)
    message(STATUS "CKTSO enabled via path: ${PATH_CKTSO}")
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
    target_compile_definitions(lightsim2grid_cpp PRIVATE NICSLU_SOLVER_AVAILABLE)
    
    # Ajout des headers et libs externes
    target_include_directories(lightsim2grid_cpp PRIVATE ${PATH_NICSLU}/include)
    
    # Configuration du Linker
    # Similaire à ce que faisait setup.py avec library_dirs
    target_link_directories(lightsim2grid_cpp PRIVATE ${PATH_NICSLU}/lib)
    # Assurez-vous de linker la lib nicslu (ajustez le nom 'nicslu' selon le nom réel du fichier .so/.a)
    target_link_libraries(lightsim2grid_cpp PRIVATE nicslu) 
    
    # Ajout d'un fichier source supplémentaire si nécessaire
    # target_sources(lightsim2grid_cpp PRIVATE src/NicsluInterface.cpp)
endif()

# 2. Configuration CKTSO
if(ENABLE_CKTSO)
    target_compile_definitions(lightsim2grid_cpp PRIVATE CKTSO_SOLVER_AVAILABLE)
    target_include_directories(lightsim2grid_cpp PRIVATE ${PATH_CKTSO}/include)
    target_link_directories(lightsim2grid_cpp PRIVATE ${PATH_CKTSO}/lib)
    target_link_libraries(lightsim2grid_cpp PRIVATE cktso) # Nom de lib à vérifier
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
