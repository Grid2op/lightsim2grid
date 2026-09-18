# Match lightsim2grid_core's -march=native / -O3 on a plugin target.
#
# lightsim2grid_core and a plugin module (e.g. this one) are two *separate*
# shared libraries. -march=native changes which SIMD instruction sets Eigen
# sees enabled (__AVX__, __AVX2__, ...), which changes EIGEN_MAX_ALIGN_BYTES
# and thus how Eigen aligns/allocates/frees its dynamic-size matrices. Any
# owning Eigen object (RealVect, CplxVect, ...) allocated under one alignment
# assumption and freed under another -- which can happen the moment such an
# object crosses the BaseAlgo virtual interface -- silently corrupts the
# heap ("double free or corruption", usually far away from the real cause).
# See docs/solver_plugin.rst for the full writeup.
#
# Usage, after find_package(lightsim2grid_core ...) and add_library(<target> ...):
#   include("${CMAKE_CURRENT_SOURCE_DIR}/../cmake/MatchLightsim2gridBuildFlags.cmake")
#   ls2g_match_core_build_flags(<target>)
function(ls2g_match_core_build_flags target)
    if(lightsim2grid_core_FOUND)
        # Strategy 1 (installed package): the CMake config exports the
        # flags lightsim2grid_core was actually built with -- no guessing.
        set(_ls2g_march_native "${lightsim2grid_core_MARCH_NATIVE}")
        set(_ls2g_o3_optim     "${lightsim2grid_core_O3_OPTIM}")
    else()
        # Strategy 2 (source tree): no installed package to query -- fall
        # back to the same env vars lightsim2grid itself reads, so building
        # both the plugin and lightsim2grid from source with
        # `__COMPILE_MARCHNATIVE=1` set stays consistent.
        if("$ENV{__COMPILE_MARCHNATIVE}" STREQUAL "1" OR "$ENV{__COMPILE_MARCHNATIVE}" STREQUAL "True")
            set(_ls2g_march_native ON)
        else()
            set(_ls2g_march_native OFF)
        endif()
        # -O3 is ON by default in lightsim2grid_core (src/core/CMakeLists.txt):
        # only __O3_OPTIM=0 / False turns it off, an unset env var means ON.
        if("$ENV{__O3_OPTIM}" STREQUAL "0" OR "$ENV{__O3_OPTIM}" STREQUAL "False")
            set(_ls2g_o3_optim OFF)
        else()
            set(_ls2g_o3_optim ON)
        endif()
    endif()

    # -march=native has no MSVC equivalent; lightsim2grid_core itself only
    # ever applies it under `if(NOT MSVC)` (see src/core/CMakeLists.txt), so
    # a Windows-built lightsim2grid_core is never USE_MARCH_NATIVE in the
    # first place.
    if(NOT MSVC)
        if(_ls2g_march_native)
            message(STATUS "${target}: lightsim2grid_core was built with -march=native -- matching it (required, see docs/solver_plugin.rst)")
            target_compile_options(${target} PRIVATE -march=native)
        endif()
        if(_ls2g_o3_optim)
            message(STATUS "${target}: lightsim2grid_core was built with -O3 -- matching it")
            target_compile_options(${target} PRIVATE -O3)
        endif()
    endif()
endfunction()
