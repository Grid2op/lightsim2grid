
set(BINDING_SOURCES
    "bindings/binding_module.cpp"
    "bindings/binding_enums.cpp"
    "bindings/binding_solvers.cpp"
    "bindings/binding_containers.cpp"
    "bindings/binding_misc.cpp"
    "bindings/binding_gridmodel.cpp"
    "bindings/binding_batch.cpp"
)

set(SOURCES
    "help_fun_msg.cpp"
    "BaseConstants.cpp"
    "GridModel.cpp"
    "ChooseSolver.cpp"
    "SolverRegistry.cpp"
    "SolverTypeNames.cpp"
    "BuiltinSolversRegistration.cpp"
    "SolverInstantiations.cpp"
    "Utils.cpp"
    "DataConverter.cpp"
    "SubstationContainer.cpp"
)
add_subdirectory(linear_solvers)
add_subdirectory(batch_algorithm)
add_subdirectory(element_container)
add_subdirectory(powerflow_algorithm)


# ---- Pure C++ shared library (no Python dependency) ----
add_library(lightsim2grid_core SHARED ${SOURCES})

set_target_properties(lightsim2grid_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(lightsim2grid_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/lightsim2grid>
)

# LS2G_BUILDING_CORE must be defined when compiling the core library
# so that LS2G_API expands to __declspec(dllexport) on Windows.
target_compile_definitions(lightsim2grid_core PRIVATE LS2G_BUILDING_CORE)


# ---- pybind11 wrapper — thin glue layer linking against core ----
pybind11_add_module(lightsim2grid_cpp ${BINDING_SOURCES})

target_include_directories(lightsim2grid_cpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings
)

target_link_libraries(lightsim2grid_cpp PRIVATE lightsim2grid_core)

# On Linux/macOS, embed $ORIGIN so both .so files find each other regardless of
# where the wheel is installed (site-packages, venv, editable source tree, etc.)
if(NOT WIN32)
    set_target_properties(lightsim2grid_core PROPERTIES
        INSTALL_RPATH "$ORIGIN"
        BUILD_RPATH   "$ORIGIN"
    )
    set_target_properties(lightsim2grid_cpp PROPERTIES
        INSTALL_RPATH "$ORIGIN"
        BUILD_RPATH   "$ORIGIN"
    )
endif()

# Windows: export symbols from core and generate import lib so plugin DLLs
# can link against lightsim2grid_core.
if(WIN32)
    set_target_properties(lightsim2grid_core PROPERTIES ENABLE_EXPORTS ON)
endif()
