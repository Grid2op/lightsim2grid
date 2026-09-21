# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
#
# Standalone build of the profiling drivers (cached powerflow, batch algorithms). Needs no python,
# no pybind11 and no pip: it links the core library straight from src/core.
#
#   cmake -S benchmarks/cache_profiling -B build_profile -DCMAKE_BUILD_TYPE=Release
#   cmake --build build_profile -j
#
# The drivers are meant to be run under callgrind (see run_profile.sh and
# run_profile_batch.sh), so the
# build keeps -O3 (what the wheels ship) and adds -g: debug info changes no
# codegen but is what lets callgrind_annotate attribute instructions to source
# lines instead of to bare addresses.

cmake_minimum_required(VERSION 3.16...4.99)
project(ls2g_cache_profiling VERSION 0.0.0 LANGUAGES CXX)

if("cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(CMAKE_CXX_STANDARD 17)
else()
    set(CMAKE_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../eigen")

if(NOT TARGET lightsim2grid_core)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../src/core"
                     "${CMAKE_CURRENT_BINARY_DIR}/lightsim2grid_core")
endif()

# Two drivers, built the same way: the single-solve audit and the batch one.
# -g everywhere (core included) so callgrind can name source lines. Frame
# pointers keep the call graph honest under -O3.
if(TARGET lightsim2grid_core)
    get_target_property(_imported lightsim2grid_core IMPORTED)
    if(NOT _imported)
        target_compile_options(lightsim2grid_core PRIVATE -g -fno-omit-frame-pointer)
    endif()
endif()

# valgrind/callgrind.h -- the drivers degrade to a plain timing run without it.
include(CheckIncludeFileCXX)
check_include_file_cxx("valgrind/callgrind.h" HAVE_CALLGRIND_H)
if(NOT HAVE_CALLGRIND_H)
    message(WARNING "valgrind/callgrind.h not found: the drivers will run but "
                    "will not be able to restrict callgrind's collection to the "
                    "solves themselves. Install valgrind's headers.")
endif()

foreach(driver profile_cached_pf profile_batch)
    add_executable(${driver} "${CMAKE_CURRENT_SOURCE_DIR}/${driver}.cpp")
    target_link_libraries(${driver} PRIVATE lightsim2grid_core)
    target_compile_options(${driver} PRIVATE -g -fno-omit-frame-pointer)
    if(HAVE_CALLGRIND_H)
        target_compile_definitions(${driver} PRIVATE LS2G_HAS_CALLGRIND=1)
    endif()
endforeach()
