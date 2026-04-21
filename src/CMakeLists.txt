
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


include_directories(".")
include_directories("bindings")

# python_add_library(lightsim2grid_cpp MODULE ${SOURCES} WITH_SOABI)

# target_include_directories(lightsim2grid_cpp PRIVATE 
    # ${CMAKE_CURRENT_SOURCE_DIR}/src
# )


pybind11_add_module(lightsim2grid_cpp ${BINDING_SOURCES} ${SOURCES})

# Windows: mark exported symbols with __declspec(dllexport) and generate an
# import library (.lib) so that plugin DLLs can link against this module.
if(WIN32)
    target_compile_definitions(lightsim2grid_cpp PRIVATE LS2G_BUILDING_CORE)
    set_target_properties(lightsim2grid_cpp PROPERTIES ENABLE_EXPORTS ON)
endif()

# target_link_libraries(lightsim2grid_cpp PRIVATE pybind11::headers)



# python_add_library(lightsim2grid_cpp MODULE ${SOURCES} WITH_SOABI)
# target_link_libraries(lightsim2grid_cpp PRIVATE pybind11::headers)
