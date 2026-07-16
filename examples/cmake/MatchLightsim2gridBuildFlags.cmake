# Match a solver plugin target's -march=native / -O3 flags to whatever
# lightsim2grid_core was actually built with.
#
# This is a correctness requirement, not a performance nicety:
# lightsim2grid_core and a plugin module are two *separate* shared
# libraries. -march=native changes which SIMD ISA Eigen sees enabled
# (__AVX__, __AVX2__, ...), which changes EIGEN_MAX_ALIGN_BYTES and thus
# how Eigen aligns/allocates/frees its dynamic-size matrices (RealVect,
# CplxVect, ...). If the two libraries disagree, an Eigen object allocated
# under one alignment assumption and freed under another (possible the
# moment such an object crosses the BaseAlgo interface) silently corrupts
# the heap. See docs/solver_plugin.rst, section
# "Matching build flags (-march=native / -O3)".
function(ls2g_match_core_build_flags target)
    if(lightsim2grid_core_FOUND)
        # Installed package: it exports the flags it was actually built
        # with -- no guessing needed.
        set(_ls2g_march_native "${lightsim2grid_core_MARCH_NATIVE}")
        set(_ls2g_o3_optim     "${lightsim2grid_core_O3_OPTIM}")
    else()
        # Source tree: no installed package to query -- fall back to the
        # same env vars lightsim2grid itself reads.
        if("$ENV{__COMPILE_MARCHNATIVE}" STREQUAL "1" OR "$ENV{__COMPILE_MARCHNATIVE}" STREQUAL "True")
            set(_ls2g_march_native ON)
        else()
            set(_ls2g_march_native OFF)
        endif()
        if("$ENV{__O3_OPTIM}" STREQUAL "1" OR "$ENV{__O3_OPTIM}" STREQUAL "True")
            set(_ls2g_o3_optim ON)
        else()
            set(_ls2g_o3_optim OFF)
        endif()
    endif()

    if(NOT MSVC)
        if(_ls2g_march_native)
            message(STATUS "${target}: lightsim2grid_core was built with -march=native -- matching it")
            target_compile_options(${target} PRIVATE -march=native)
        endif()
        if(_ls2g_o3_optim)
            message(STATUS "${target}: lightsim2grid_core was built with -O3 -- matching it")
            target_compile_options(${target} PRIVATE -O3)
        endif()
    endif()
endfunction()
