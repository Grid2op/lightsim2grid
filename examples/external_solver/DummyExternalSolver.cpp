// Minimal lightsim2grid solver plugin example.
//
// To register this solver with the C++ SolverRegistry at dlopen() time, we
// declare a static SolverRegistrar in an anonymous namespace.  No macro is
// needed — the object's constructor fires when the shared library is loaded,
// which registers "DummyExternal" into the singleton registry.
//
// Build:
//   mkdir build && cd build
//   cmake .. -DLIGHTSIM2GRID_SRC=<path/to/lightsim2grid/src>
//              -DEigen3_INCLUDE=<path/to/eigen>
//   make
//
// Python usage (from examples/external_solver/):
//   import lightsim2grid
//   lightsim2grid.load_solver_plugin("build/libdummy_solver.so")
//   grid.change_solver("DummyExternal")

#include <SolverRegistry.hpp>
#include <powerflow_algorithm/BaseAlgo.hpp>

// ---------------------------------------------------------------------------
// A trivial AC solver that always "converges" on the first call by returning
// the initial voltage vector unchanged.  Useful as a smoke-test for the
// plugin mechanism; not suitable for real power-flow calculations.
// ---------------------------------------------------------------------------
class DummyExternalSolver : public BaseAlgo {
public:
    DummyExternalSolver() : BaseAlgo(/*is_ac=*/true) {}

    bool compute_pf(
        const Eigen::SparseMatrix<cplx_type>& /*Ybus*/,
        CplxVect& V,
        const CplxVect& /*Sbus*/,
        Eigen::Ref<const IntVect> /*slack_ids*/,
        const RealVect& /*slack_weights*/,
        Eigen::Ref<const IntVect> /*pv*/,
        Eigen::Ref<const IntVect> /*pq*/,
        int /*max_iter*/,
        real_type /*tol*/) override
    {
        // Store V unchanged and claim convergence.
        V_ = V;
        Va_ = V.array().arg();
        Vm_ = V.array().abs();
        n_  = static_cast<int>(V.size());
        nr_iter_ = 1;
        err_ = ErrorType::NoError;
        return true;
    }
};

// ---------------------------------------------------------------------------
// Self-registration: fires when dlopen() / LoadLibrary() maps this .so.
// ---------------------------------------------------------------------------
namespace {
    SolverRegistrar _dummy_registrar(
        "DummyExternal",
        []{ return std::unique_ptr<BaseAlgo>(new DummyExternalSolver()); }
    );
}
