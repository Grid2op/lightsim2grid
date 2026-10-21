// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// NROuterRebuild_* plugin: OpenLoadFlow's outer loops (the NROuter_* driver and loops) on a
// Newton whose Jacobian is rebuilt -- and analyzed again -- whenever the active equations change,
// as OpenLoadFlow does, see RebuildOuterInner.hpp. Registers through LS2G_PLUGIN_ENTRY:
//   "NROuterRebuild_SparseLU"
//   "NROuterRebuild_KLU"       (only if KLU is available)
//
// Build (after pip install lightsim2grid, or from a source-tree checkout):
//   LS2G_CMAKE=$(python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())")
//   cmake -S examples/rebuild_outer_algorithm -B examples/rebuild_outer_algorithm/build -DLIGHTSIM2GRID_CMAKE_DIR="$LS2G_CMAKE"
//   cmake --build examples/rebuild_outer_algorithm/build
//
// Python usage:
//   import lightsim2grid
//   lightsim2grid.load_algorithm_plugin("examples/rebuild_outer_algorithm/build/librebuild_outer_algorithm.so")
//   grid.change_algorithm("NROuterRebuild_KLU")   # then add_outer_loop(...) as with NROuter_KLU

#include <AlgorithmRegistry.hpp>

#include "RebuildOuterInner.hpp"

static void register_plugin(ls2g::PluginRegistrar& reg) {
    reg.add("NROuterRebuild_SparseLU",
        []{ return std::unique_ptr<ls2g::BaseAlgo>(
                new ls2g::OuterLoopAlgo<ls2g::RebuildOuterInner<ls2g::LinearSolverPolicy<ls2g::SparseLULinearSolver> > >()); });

#ifdef KLU_SOLVER_AVAILABLE
    reg.add("NROuterRebuild_KLU",
        []{ return std::unique_ptr<ls2g::BaseAlgo>(
                new ls2g::OuterLoopAlgo<ls2g::RebuildOuterInner<ls2g::LinearSolverPolicy<ls2g::KLULinearSolver> > >()); });
#endif
}
LS2G_PLUGIN_ENTRY(register_plugin)
