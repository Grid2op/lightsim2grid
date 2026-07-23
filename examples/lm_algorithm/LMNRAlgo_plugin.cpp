// LMNRAlgo plugin registration.
//
// Registers two names with the C++ AlgorithmRegistry through the exported
// ls2g_register_plugin entry point (LS2G_PLUGIN_ENTRY -- same mechanism as
// examples/dist_slack_algorithm/NRAlgoDistSlack_plugin.cpp and
// examples/external_algorithm/DummyExternalAlgo.cpp):
//   "NR_LM_SparseLU" -> LMNRAlgo<LinearSolverPolicy<SparseLULinearSolver>, SingleSlackNRSystem>
//   "NR_LM_KLU"      -> LMNRAlgo<LinearSolverPolicy<KLULinearSolver>, SingleSlackNRSystem>  (only if KLU is available)
//
// Both names are registered atomically: if either is already taken, neither is
// added and load_algorithm_plugin() raises instead of half-registering.
//
// Single-slack only for this first cut; LMNRAlgo is generic over its
// NRSystem template parameter, so registering MultiSlackNRSystem variants
// later is a one-line addition.
//
// Build (after pip install lightsim2grid, or from a source-tree checkout):
//   LS2G_CMAKE=$(python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())")
//   cmake -S examples/lm_algorithm -B examples/lm_algorithm/build -DLIGHTSIM2GRID_CMAKE_DIR="$LS2G_CMAKE"
//   cmake --build examples/lm_algorithm/build
//
// Python usage:
//   import lightsim2grid
//   lightsim2grid.load_algorithm_plugin("examples/lm_algorithm/build/liblm_algorithm.so")
//   grid.change_algorithm("NR_LM_KLU")

#include <AlgorithmRegistry.hpp>

#include "LMNRAlgo.hpp"

static void register_plugin(ls2g::PluginRegistrar& reg) {
    reg.add("NR_LM_SparseLU",
        []{ return std::unique_ptr<ls2g::BaseAlgo>(
            new ls2g::LMNRAlgo<ls2g::LinearSolverPolicy<ls2g::SparseLULinearSolver>, ls2g::SingleSlackNRSystem>()); });

#ifdef KLU_SOLVER_AVAILABLE
    reg.add("NR_LM_KLU",
        []{ return std::unique_ptr<ls2g::BaseAlgo>(
            new ls2g::LMNRAlgo<ls2g::LinearSolverPolicy<ls2g::KLULinearSolver>, ls2g::SingleSlackNRSystem>()); });
#endif
}
LS2G_PLUGIN_ENTRY(register_plugin)
