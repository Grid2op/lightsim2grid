// NRAlgoDistSlack plugin registration.
//
// Registers two names with the C++ AlgorithmRegistry at dlopen() time (via a
// static AlgorithmRegistrar in an anonymous namespace -- same mechanism as
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

namespace {
    ls2g::AlgorithmRegistrar _dist_slack_sparselu_registrar(
        "NRDistSlack_SparseLU",
        []{ return std::unique_ptr<ls2g::BaseAlgo>(new ls2g::NRAlgoDistSlack<ls2g::SparseLULinearSolver>()); }
    );

#ifdef KLU_SOLVER_AVAILABLE
    ls2g::AlgorithmRegistrar _dist_slack_klu_registrar(
        "NRDistSlack_KLU",
        []{ return std::unique_ptr<ls2g::BaseAlgo>(new ls2g::NRAlgoDistSlack<ls2g::KLULinearSolver>()); }
    );
#endif
}
