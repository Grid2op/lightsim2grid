
set(SOURCES
    "main.cpp"
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

# python_add_library(lightsim2grid_cpp MODULE ${SOURCES} WITH_SOABI)

# target_include_directories(lightsim2grid_cpp PRIVATE 
    # ${CMAKE_CURRENT_SOURCE_DIR}/src
# )


pybind11_add_module(lightsim2grid_cpp ${SOURCES})

# Windows: mark exported symbols with __declspec(dllexport) and generate an
# import library (.lib) so that plugin DLLs can link against this module.
if(WIN32)
    target_compile_definitions(lightsim2grid_cpp PRIVATE LS2G_BUILDING_CORE)
    set_target_properties(lightsim2grid_cpp PROPERTIES ENABLE_EXPORTS ON)
endif()

# target_link_libraries(lightsim2grid_cpp PRIVATE pybind11::headers)



# python_add_library(lightsim2grid_cpp MODULE ${SOURCES} WITH_SOABI)
# target_link_libraries(lightsim2grid_cpp PRIVATE pybind11::headers)
