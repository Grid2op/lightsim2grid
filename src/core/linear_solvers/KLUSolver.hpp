// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifdef KLU_SOLVER_AVAILABLE
#ifndef KLSOLVER_H
#define KLSOLVER_H
#include <iostream>
#include <memory>

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Utils.hpp"

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

namespace ls2g {

// import klu package
extern "C" {
    #include "cs.h"
    #include "klu.h"
}

/**
class to handle the solver using newton-raphson method, using KLU algorithm and sparse matrices.

As long as the admittance matrix of the sytem does not change, you can reuse the same solver.
Reusing the same solver is possible, but "reset" method must be called.

Otherwise, unexpected behaviour might follow, including "segfault".

**/
class LS2G_API KLULinearSolver final
{
    public:
        KLULinearSolver() noexcept :
            common_(),
            symbolic_(nullptr, SymbolicDeleter{&common_}),
            numeric_(nullptr, NumericDeleter{&common_})
            { klu_defaults(&common_); }

        // symbolic_ / numeric_ are unique_ptr with custom deleters and free themselves.
        ~KLULinearSolver() noexcept = default;

        // public api
        ErrorType reset();
        ErrorType analyze(const EigenRefConstRealSpMat & J);   // reordering + symbolic factorization (structure only)
        ErrorType factorize(const EigenRefConstRealSpMat & J); // numeric factorization (requires values)
        ErrorType refactorize(const EigenRefConstRealSpMat & J);  // re-numeric factorization, reuses symbolic
        ErrorType solve(Eigen::Ref<RealVect> b);

        // can this linear solver solve problem where RHS is a matrix
        static const bool CAN_SOLVE_MAT;

    private:
        // KLU frees its symbolic / numeric handles through a pair of functions that
        // also need the `klu_common` control struct, so the deleters are stateful and
        // hold a pointer to `common_`. `common_` is therefore declared *before* the two
        // handles: members are destroyed in reverse declaration order, so the handles
        // (and their deleters, which dereference `common_`) are released while `common_`
        // is still alive.
        struct SymbolicDeleter {
            klu_common* common_ = nullptr;
            void operator()(klu_symbolic* ptr) const noexcept {
                if(ptr != nullptr){
                    klu_symbolic* tmp = ptr;  // klu_free_symbolic takes a klu_symbolic**
                    klu_free_symbolic(&tmp, common_);
                }
            }
        };
        struct NumericDeleter {
            klu_common* common_ = nullptr;
            void operator()(klu_numeric* ptr) const noexcept {
                if(ptr != nullptr){
                    klu_numeric* tmp = ptr;  // klu_free_numeric takes a klu_numeric**
                    klu_free_numeric(&tmp, common_);
                }
            }
        };

        // solver initialization (declaration order matters, see note above)
        klu_common common_;
        std::unique_ptr<klu_symbolic, SymbolicDeleter> symbolic_;
        std::unique_ptr<klu_numeric, NumericDeleter> numeric_;

        // no copy allowed
        KLULinearSolver(const KLULinearSolver&) = delete;
        KLULinearSolver(KLULinearSolver&&) = delete;
        KLULinearSolver & operator=(KLULinearSolver&&) = delete;
        KLULinearSolver & operator=(const KLULinearSolver&) = delete;
};

#endif // KLSOLVER_H

} // namespace ls2g

#elif defined(_READ_THE_DOCS)
#ifndef KLSOLVER_H
#define KLSOLVER_H
#include "Utils.hpp"

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

namespace ls2g {

/**
Doc-only stand-in for KLULinearSolver, used when the real KLU (SuiteSparse) library is
not available but _READ_THE_DOCS is set: it gives NR_KLU (and its DC_KLU / FDPF_*_KLU /
NRSing_KLU / NRRefactorRetry_KLU siblings, see Solvers.hpp) a genuine, distinct,
fully-functional C++ type -- with the same name and public interface as the real
KLULinearSolver -- so Sphinx can document them as real classes without requiring the
actual library. Every method is a trivial no-op: this is never meant to actually solve
a real system, only to be introspected.
**/
class LS2G_API KLULinearSolver final
{
    public:
        KLULinearSolver() noexcept = default;
        ~KLULinearSolver() noexcept = default;

        ErrorType reset() { return ErrorType::NoError; }
        ErrorType analyze(const EigenRefConstRealSpMat & /*J*/) { return ErrorType::NoError; }
        ErrorType factorize(const EigenRefConstRealSpMat & /*J*/) { return ErrorType::NoError; }
        ErrorType refactorize(const EigenRefConstRealSpMat & /*J*/) { return ErrorType::NoError; }
        ErrorType solve(Eigen::Ref<RealVect> /*b*/) { return ErrorType::NoError; }

        // can this linear solver solve problem where RHS is a matrix
        static constexpr bool CAN_SOLVE_MAT = false;

    private:
        // no copy allowed (matches the real KLULinearSolver)
        KLULinearSolver(const KLULinearSolver&) = delete;
        KLULinearSolver(KLULinearSolver&&) = delete;
        KLULinearSolver & operator=(KLULinearSolver&&) = delete;
        KLULinearSolver & operator=(const KLULinearSolver&) = delete;
};

} // namespace ls2g
#endif // KLSOLVER_H
#endif // KLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)