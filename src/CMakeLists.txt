
set(SOURCES
    "main.cpp"
    "help_fun_msg.cpp"
    "BaseConstants.cpp"
    "GridModel.cpp"
    "ChooseSolver.cpp"
    "Solvers.cpp"
    "Utils.cpp"
    "DataConverter.cpp"
    "SubstationContainer.cpp"
)
message(SOURCES 1="${SOURCES}")
add_subdirectory(linear_solvers)
message(SOURCES 2="${SOURCES}")
add_subdirectory(batch_algorithm)
add_subdirectory(element_container)
add_subdirectory(powerflow_algorithm)


include_directories(".")

# python_add_library(lightsim2grid_cpp MODULE ${SOURCES} WITH_SOABI)

# target_include_directories(lightsim2grid_cpp PRIVATE 
    # ${CMAKE_CURRENT_SOURCE_DIR}/src
# )


pybind11_add_module(lightsim2grid_cpp ${SOURCES})

# target_link_libraries(lightsim2grid_cpp PRIVATE pybind11::headers)



# python_add_library(lightsim2grid_cpp MODULE ${SOURCES} WITH_SOABI)
# target_link_libraries(lightsim2grid_cpp PRIVATE pybind11::headers)
