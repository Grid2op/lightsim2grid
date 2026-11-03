# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

# The A/B module of README.md: lightsim2grid_core (as src/tests builds it) and one pybind11
# module binding its result getters with each return value policy.
#
#   cmake -S benchmarks/view_policy_ab -B build_ab -DCMAKE_BUILD_TYPE=Release \
#         -Dpybind11_DIR=$(python -m pybind11 --cmakedir)
#   cmake --build build_ab -j
cmake_minimum_required(VERSION 3.16...4.99)
# the VERSION is required: src/core generates its CMake package version file from it
project(lightsim2grid_view_policy_ab VERSION 0.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../eigen")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../src/core" "${CMAKE_CURRENT_BINARY_DIR}/core")
find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(ab_views ab_views.cpp)
target_link_libraries(ab_views PRIVATE lightsim2grid_core)
