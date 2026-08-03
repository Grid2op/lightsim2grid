// NRAlgoDistSlack plugin registration.
//
// Registers two names with the C++ AlgorithmRegistry through the exported
// ls2g_register_plugin entry point (LS2G_PLUGIN_ENTRY -- same mechanism as
// examples/external_algorithm/DummyExternalAlgo.cpp):
//   "NRDistSlack_SparseLU" -> NRAlgoDistSlack<SparseLULinearSolver>
//   "NRDistSlack_KLU"      -> NRAlgoDistSlack<KLULinearSolver>  (only if KLU is available)
//
// Build (after pip install lightsim2grid, or from a source-tree checkout):
//   LS2G_CMAKE=$(python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())")
//   cmake -S examples/dist_slack_algorithm -B examples/dist_slack_algorithm/build -DLIGHTSIM2GRID_CMAKE_DIR="$LS2G_CMAKE"
//   cmake --build examples/dist_slack_algorithm/build
//
// Python usage:
//   import lightsim2grid
//   lightsim2grid.load_algorithm_plugin("examples/dist_slack_algorithm/build/libdist_slack_algorithm.so")
//   grid.change_algorithm("NRDistSlack_KLU")

#include <AlgorithmRegistry.hpp>

#include "NRAlgoDistSlack.hpp"

static void register_plugin(ls2g::PluginRegistrar& reg) {
    reg.add("NRDistSlack_SparseLU",
        []{ return std::unique_ptr<ls2g::BaseAlgo>(new ls2g::NRAlgoDistSlack<ls2g::SparseLULinearSolver>()); });

#ifdef KLU_SOLVER_AVAILABLE
    reg.add("NRDistSlack_KLU",
        []{ return std::unique_ptr<ls2g::BaseAlgo>(new ls2g::NRAlgoDistSlack<ls2g::KLULinearSolver>()); });
#endif
}
LS2G_PLUGIN_ENTRY(register_plugin)
