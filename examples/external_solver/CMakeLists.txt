cmake_minimum_required(VERSION 3.15)
project(dummy_solver CXX)

# -----------------------------------------------------------------------
# Locate the installed lightsim2grid headers.
# Pass -DLIGHTSIM2GRID_SRC=<path> on the cmake command line if the source
# tree is not at the default relative location.
# -----------------------------------------------------------------------
if(NOT DEFINED LIGHTSIM2GRID_SRC)
    # Default: two levels up from this CMakeLists.txt, then into src/
    set(LIGHTSIM2GRID_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../../src")
endif()

if(NOT DEFINED Eigen3_INCLUDE)
    set(Eigen3_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../../eigen")
endif()

if(NOT EXISTS "${LIGHTSIM2GRID_SRC}/SolverRegistry.hpp")
    message(FATAL_ERROR
        "lightsim2grid src headers not found at '${LIGHTSIM2GRID_SRC}'.\n"
        "Pass -DLIGHTSIM2GRID_SRC=<path> to cmake.")
endif()

# Use the same C++ standard as the parent build.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The plugin is a MODULE (loadable at runtime via dlopen).
add_library(dummy_solver MODULE DummyExternalSolver.cpp)

# Expose lightsim2grid src/ headers directly so #include <SolverRegistry.hpp> works,
# and expose the repo root so #include <src/SolverRegistry.hpp> also works.
target_include_directories(dummy_solver PRIVATE
    "${LIGHTSIM2GRID_SRC}"               # #include <SolverRegistry.hpp>
    "${CMAKE_CURRENT_SOURCE_DIR}/../.."  # #include <src/SolverRegistry.hpp>
    "${Eigen3_INCLUDE}"
)

if(WIN32)
    # On Windows the plugin DLL must link against lightsim2grid_cpp's import
    # library (.lib), which is installed alongside lightsim2grid_cpp.pyd.
    # We try to locate it automatically via Python; override with
    # -DLIGHTSIM2GRID_CPP_LIB=<path> if auto-detection fails.
    if(NOT DEFINED LIGHTSIM2GRID_CPP_LIB)
        find_package(Python REQUIRED COMPONENTS Interpreter)
        execute_process(
            COMMAND "${Python_EXECUTABLE}" -c
                "import importlib.util, pathlib; \
                 spec = importlib.util.find_spec('lightsim2grid.lightsim2grid_cpp'); \
                 print(pathlib.Path(spec.origin).parent / 'lightsim2grid_cpp.lib')"
            OUTPUT_VARIABLE _ls2g_lib
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(EXISTS "${_ls2g_lib}")
            set(LIGHTSIM2GRID_CPP_LIB "${_ls2g_lib}")
            message(STATUS "lightsim2grid_cpp.lib found at: ${LIGHTSIM2GRID_CPP_LIB}")
        else()
            message(FATAL_ERROR
                "lightsim2grid_cpp.lib not found (tried '${_ls2g_lib}').\n"
                "Pass -DLIGHTSIM2GRID_CPP_LIB=<path> to cmake.\n"
                "The .lib is installed alongside lightsim2grid_cpp.pyd when "
                "lightsim2grid is built from source on Windows.")
        endif()
    endif()
    target_link_libraries(dummy_solver PRIVATE "${LIGHTSIM2GRID_CPP_LIB}")
    # Windows default: no prefix, .dll suffix — no extra properties needed.
elseif(UNIX AND NOT APPLE)
    # Undefined symbols are resolved at load time from the main .so.
    target_link_options(dummy_solver PRIVATE -Wl,--allow-shlib-undefined)
    set_target_properties(dummy_solver PROPERTIES PREFIX "lib" SUFFIX ".so")
elseif(APPLE)
    target_link_options(dummy_solver PRIVATE -undefined dynamic_lookup)
    set_target_properties(dummy_solver PROPERTIES PREFIX "lib" SUFFIX ".so")
endif()
