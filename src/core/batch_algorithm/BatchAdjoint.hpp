// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BATCH_ADJOINT_H
#define BATCH_ADJOINT_H

#include <cstddef>
#include <cstring>
#include <exception>   // std::exception_ptr: libc++ only forward-declares it in <vector>
#include <memory>
#include <vector>

#include "Eigen/Core"
#include "Eigen/SparseCore"

#include "Utils.hpp"
#include "linear_solvers/LinearSolverPolicy.hpp"
#include "linear_solvers/RefactorRetryLinearSolver.hpp"
#include "linear_solvers/SparseLUSolver.hpp"
#include "linear_solvers/KLUSolver.hpp"

namespace ls2g {

/**
Reverse-mode differentiation of a batch of powerflows: the adjoint system.

Every row of a batch solves G(x_i ; u_i) = 0 by Newton-Raphson and ends on a
Jacobian J_i. Differentiating a scalar loss L(V_1, ..., V_n) with respect to the
per-row inputs needs, for each row, the solution of

    J_i^T lambda_i = xbar_i

(the implicit function theorem), and nothing else: lambda_i does not depend on
which input is being differentiated, so ONE transposed solve per row yields the
gradient with respect to every injection at once. See docs/ for the full
derivation and for how lambda maps back onto the elements.

Two properties of this repository's batch algorithms are what make that cheap,
and they are why this class exists rather than a generic autodiff layer:

- **the Jacobian's sparsity pattern is the same for every row of a batch** (see
  BaseBatchSweep's class comment: a row may change a bus's role, but only by
  rewriting values inside a pattern reserved up front). So one pattern is
  snapshotted once and every row is just `nnz_J` reals;
- **KLU solves the transposed system out of the factorization of J itself**
  (`klu_tsolve`, see KLULinearSolver::solve_transpose). No transposed copy of J,
  no second symbolic analysis.

Memory is the price, and it is known before the batch runs: `nb_rows * nnz_J`
reals, allocated ONCE by allocate() before the row loop starts (see
memory_bytes(), which a caller can consult beforehand). Each worker thread then
memcpy's into its own rows -- J_values_ is row-major precisely so that a row is
one contiguous block and threads never share a cache line except at the two
boundaries of their range.

What is NOT stored is the numeric factorization: KLU's `klu_refactor` overwrites
its klu_numeric in place, and there is no way to clone one, so keeping a
factorization per row would mean a full `klu_factor` (more expensive than a
refactor) per row in the forward. backward therefore pays one refactorize plus
one transposed solve per row -- both timed separately in the returned stats, so
the trade can be revisited on measurements rather than on belief.

This class knows nothing about powerflows: it is a batch of "same pattern,
different values" linear systems and their transposed solves. Which rows carry
identity equations, and what lambda means once solved, are the caller's business
(see BaseBatchSweep::solve_JT).
**/
class LS2G_API BatchAdjoint
{
    public:
        // row-major: one row of J's values, or of xbar / lambda, is contiguous, so a
        // worker writes a block instead of striding across the whole matrix (and so
        // the numpy array a binding hands out is a view, not a transposing copy).
        using RealMatRM = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

        BatchAdjoint() = default;

        // --- filled during the forward -------------------------------------------

        /** Snapshot `J`'s pattern and size the value buffer for `nb_rows` rows. Called
        once, before the row loop, with the Jacobian of the batch's warm-up solve: by
        then the ledger is built and dim_J / nnz_J are known. Drops whatever a previous
        batch left. */
        void allocate(Eigen::Index nb_rows, const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & J);

        /** Store row `row`'s converged Jacobian. Thread-safe with respect to other rows
        (disjoint writes into a buffer that is never resized here). Throws if `J`'s
        pattern is not the snapshotted one -- that would silently produce wrong
        gradients, and the check costs a memcmp against a copy the row is about to
        make anyway. */
        void store_row(Eigen::Index row, const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & J);

        void clear();

        // --- consumed by the backward --------------------------------------------

        /** Solve `J_i^T lambda_i = xbar_i` for every row.

        `xbar` is `(nb_rows, nb_directions * dim_J)`: one row per batch row, holding
        `nb_directions` cotangents laid end to end (`dim_J` contiguous reals each --
        the layout `klu_tsolve`'s `nrhs` wants). Reverse-mode differentiation of a
        scalar loss uses one direction; several are only needed for a full Jacobian,
        and cost one triangular solve each, not one factorization each.

        `identity_rows[i]` lists the equations of row `i` that the batch pinned to the
        identity (a bus stranded by a contingency, a bus still PV in this row). Their
        multiplier is dropped from the result: it belongs to a frozen equation, and
        reporting it would put a gradient on an injection that this row's solution does
        not actually depend on. May be empty (no row pins anything), otherwise it must
        have one entry per row.

        Rows the forward never stored (they diverged, or were skipped) come back as
        zero, with `row_ok` 0. */
        RealMatRM solve_JT(const Eigen::Ref<const RealMatRM> & xbar,
                           const std::vector<std::vector<int> > & identity_rows,
                           int nb_thread,
                           std::vector<char> & row_ok) const;

        // --- introspection --------------------------------------------------------

        bool is_allocated() const noexcept { return dim_J() > 0 && J_values_.rows() > 0; }
        Eigen::Index nb_rows() const noexcept { return J_values_.rows(); }
        Eigen::Index dim_J() const noexcept { return pattern_.rows(); }
        Eigen::Index nnz_J() const noexcept { return pattern_.nonZeros(); }
        bool has_row(Eigen::Index row) const {
            return row >= 0 && row < static_cast<Eigen::Index>(row_has_J_.size()) && row_has_J_[static_cast<size_t>(row)] != 0;
        }

        /** Bytes the value buffer takes for `nb_rows` rows of this Jacobian. The only
        cost of keeping every row's J, and knowable before compute() runs. */
        static std::size_t memory_bytes(Eigen::Index nb_rows, Eigen::Index nnz) noexcept {
            return static_cast<std::size_t>(nb_rows) * static_cast<std::size_t>(nnz) * sizeof(real_type);
        }
        std::size_t memory_bytes() const noexcept { return memory_bytes(nb_rows(), nnz_J()); }

        /** Linear-solver counters and timings of the last solve_JT, summed over its
        worker threads: how much of a backward went into refactorizing versus into the
        transposed solves themselves. */
        const LinearSolverStats & get_linear_solver_stats() const noexcept { return stats_; }

    private:
        // The pattern every row shares, values zeroed: a worker copies it once to get
        // a working matrix it then only ever writes values into.
        Eigen::SparseMatrix<real_type> pattern_;
        RealMatRM J_values_;              // (nb_rows, nnz_J), row-major
        std::vector<char> row_has_J_;     // char, not bool: written concurrently per row
        mutable LinearSolverStats stats_;

        void _solve_range(Eigen::Index row_begin, Eigen::Index row_end,
                          const Eigen::Ref<const RealMatRM> & xbar,
                          const std::vector<std::vector<int> > & identity_rows,
                          Eigen::Index nb_directions,
                          RealMatRM & lambda, std::vector<char> & row_ok,
                          LinearSolverStats & stats, std::exception_ptr & err) const;

        BatchAdjoint(const BatchAdjoint&) = delete;
        BatchAdjoint & operator=(const BatchAdjoint&) = delete;
};


/**
The linear solver a BatchAdjoint worker owns, behind a virtual interface.

The algorithms pick their linear solver at compile time (LinearSolverPolicy is a
template parameter, never a base pointer -- see its class comment). A batch's
algorithm, though, is chosen at RUNTIME through AlgorithmSelector, so the adjoint
cannot know its solver type statically. Three virtual calls per row against one
refactorize and one triangular solve is not a measurable cost, and it keeps the
adjoint out of the template explosion.

Only solvers that answer the transposed system belong here -- hence the
static_assert: the alternative, forming J^T and factorizing it per row, is a
different (and much more expensive) algorithm, not a fallback to slip in
silently.
**/
class IAdjointLinearSolver
{
    public:
        virtual ~IAdjointLinearSolver() = default;
        virtual ErrorType analyze(const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & J) = 0;
        virtual ErrorType factorize(const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & J) = 0;
        virtual ErrorType refactorize(const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & J) = 0;
        virtual ErrorType solve_transpose(Eigen::Ref<RealVect> b) = 0;
        virtual const LinearSolverStats & get_linear_solver_stats() const = 0;
};

template<class LinearSolver>
class AdjointLinearSolver final : public IAdjointLinearSolver
{
    static_assert(LinearSolver::CAN_SOLVE_TRANSPOSE,
                  "a BatchAdjoint solver must answer J^T x = b out of the factorization of J");
    public:
        ErrorType analyze(const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & J) override { return inner_.analyze(J); }
        ErrorType factorize(const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & J) override { return inner_.factorize(J); }
        ErrorType refactorize(const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & J) override { return inner_.refactorize(J); }
        ErrorType solve_transpose(Eigen::Ref<RealVect> b) override { return inner_.solve_transpose(b); }
        const LinearSolverStats & get_linear_solver_stats() const override { return inner_.get_linear_solver_stats(); }
    private:
        // RefactorRetry, not the plain policy: a masked or pinned row rewrites values
        // where the base matrix had its pivot, which is exactly the case KLU's fixed
        // pivot sequence finds at zero. The batch algorithms turn the same fallback on
        // for their own solver for the same reason (BaseBatchSweep.cpp).
        RefactorRetryLinearSolver<LinearSolver> inner_;
};

/** The best available transposed-solve capable linear solver: KLU where it was
compiled in, Eigen's SparseLU otherwise (always available, same answers, slower).
*/
std::unique_ptr<IAdjointLinearSolver> LS2G_API make_adjoint_linear_solver();

} // namespace ls2g

#endif // BATCH_ADJOINT_H
