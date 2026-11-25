// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// A linear solver brought in through an algorithm plugin: the algorithms of lightsim2grid
// take their linear solver as a template parameter, so instantiating them on
// MyLULinearSolver here -- from the installed headers -- gives algorithms that use it,
// without touching lightsim2grid itself. Registers through LS2G_PLUGIN_ENTRY:
//   "NR_MyLU"        Newton-Raphson, distributed slack
//   "NRSing_MyLU"    Newton-Raphson, single slack
//   "NROuter_MyLU"   Newton-Raphson with the outer loops (see docs/outer_loops.rst)
//   "DC_MyLU"        DC approximation
// The fast-decoupled algorithms (FDPF_*) cannot be registered from a plugin: part of their
// code is only compiled inside lightsim2grid_core.
//
// Build (after pip install lightsim2grid):
//   LS2G_CMAKE=$(python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())")
//   cmake -S examples/linear_solver_plugin -B examples/linear_solver_plugin/build -DLIGHTSIM2GRID_CMAKE_DIR="$LS2G_CMAKE"
//   cmake --build examples/linear_solver_plugin/build
//
// Python usage:
//   import lightsim2grid
//   lightsim2grid.load_algorithm_plugin("examples/linear_solver_plugin/build/liblinear_solver_plugin.so")
//   grid.change_algorithm("NR_MyLU")

#include <memory>

#include <AlgorithmRegistry.hpp>
// the algorithm templates, without Solvers.hpp: that one also includes the headers of
// KLU, NICSLU and CKTSO, which an installed lightsim2grid does not ship
#include <linear_solvers/LinearSolverPolicy.hpp>
#include <powerflow_algorithm/BaseDCAlgo.hpp>
#include <powerflow_algorithm/NRAlgo.hpp>
#include <powerflow_algorithm/NROuterInner.hpp>
#include <powerflow_algorithm/OuterLoopAlgo.hpp>

#include "MyLULinearSolver.hpp"

namespace {
    using MyLU = ls2g::LinearSolverPolicy<ls2g::MyLULinearSolver>;
    using NR_MyLU = ls2g::NRAlgo<MyLU, ls2g::MultiSlackNRSystem>;
    using NRSing_MyLU = ls2g::NRAlgo<MyLU, ls2g::SingleSlackNRSystem>;
    using NROuter_MyLU = ls2g::OuterLoopAlgo<ls2g::NROuterInner<MyLU>>;
    using DC_MyLU = ls2g::BaseDCAlgo<MyLU>;
}

static void register_plugin(ls2g::PluginRegistrar& reg) {
    reg.add("NR_MyLU",      []{ return std::unique_ptr<ls2g::BaseAlgo>(new NR_MyLU()); });
    reg.add("NRSing_MyLU",  []{ return std::unique_ptr<ls2g::BaseAlgo>(new NRSing_MyLU()); });
    reg.add("NROuter_MyLU", []{ return std::unique_ptr<ls2g::BaseAlgo>(new NROuter_MyLU()); });
    reg.add("DC_MyLU",      []{ return std::unique_ptr<ls2g::BaseAlgo>(new DC_MyLU()); });
}
LS2G_PLUGIN_ENTRY(register_plugin)
