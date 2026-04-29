// Minimal lightsim2grid solver plugin example.
//
// To register this solver with the C++ SolverRegistry at dlopen() time, we
// declare a static SolverRegistrar in an anonymous namespace.  No macro is
// needed — the object's constructor fires when the shared library is loaded,
// which registers "DummyExternal" into the singleton registry.
//
// Build (after pip install lightsim2grid):
//   LS2G_CMAKE=$(python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())")
//   cmake -S . -B build -DLIGHTSIM2GRID_CMAKE_DIR="$LS2G_CMAKE"
//   cmake --build build
//
// Python usage (from examples/external_solver/):
//   import lightsim2grid
//   lightsim2grid.load_solver_plugin("build/libdummy_solver.so")
//   grid.change_solver("DummyExternal")

#include <AlgorithmRegistry.hpp>
#include <powerflow_algorithm/BaseAlgo.hpp>

// ---------------------------------------------------------------------------
// A trivial AC solver that always "converges" on the first call by returning
// the initial voltage vector unchanged.  Useful as a smoke-test for the
// plugin mechanism; not suitable for real power-flow calculations.
// ---------------------------------------------------------------------------
class DummyExternalAlgo : public ls2g::BaseAlgo {
public:
    DummyExternalAlgo() : ls2g::BaseAlgo(/*is_ac=*/true) {}

    bool compute_pf(
        const Eigen::SparseMatrix<ls2g::cplx_type>& /*Ybus*/,
        ls2g::CplxVect& V,
        const ls2g::CplxVect& /*Sbus*/,
        Eigen::Ref<const ls2g::IntVect> /*slack_ids*/,
        const ls2g::RealVect& /*slack_weights*/,
        Eigen::Ref<const ls2g::IntVect> /*pv*/,
        Eigen::Ref<const ls2g::IntVect> /*pq*/,
        int /*max_iter*/,
        ls2g::real_type /*tol*/) override
    {
        // Store V unchanged and claim convergence.
        V_ = V;
        Va_ = V.array().arg();
        Vm_ = V.array().abs();
        n_  = static_cast<int>(V.size());
        nr_iter_ = 1;
        err_ = ls2g::ErrorType::NoError;
        return true;
    }
};

// ---------------------------------------------------------------------------
// Self-registration: fires when dlopen() / LoadLibrary() maps this .so.
// ---------------------------------------------------------------------------
namespace {
    ls2g::AlgorithmRegistrar _dummy_registrar(
        "DummyExternal",
        []{ return std::unique_ptr<ls2g::BaseAlgo>(new DummyExternalAlgo()); }
    );
}
