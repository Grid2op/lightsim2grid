// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BatchAdjoint.hpp"

#include <algorithm>
#include <exception>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace ls2g {

std::unique_ptr<IAdjointLinearSolver> make_adjoint_linear_solver()
{
    // the derived unique_ptr converts to the base one on return (IAdjointLinearSolver
    // has a virtual destructor), so no raw `new` has to be owned by hand here
    #ifdef KLU_SOLVER_AVAILABLE
    return std::make_unique<AdjointLinearSolver<KLULinearSolver> >();
    #else
    return std::make_unique<AdjointLinearSolver<SparseLULinearSolver> >();
    #endif
}

void BatchAdjoint::allocate(Eigen::Index nb_rows, const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & J)
{
    if(nb_rows <= 0){
        clear();
        return;
    }
    if(J.rows() != J.cols()){
        std::ostringstream exc_;
        exc_ << "BatchAdjoint::allocate: the Jacobian must be square, got " << J.rows()
             << " x " << J.cols() << ".";
        throw std::runtime_error(exc_.str());
    }
    // the pattern, values zeroed: a worker copies this to get a working matrix whose
    // index arrays it never has to rebuild
    pattern_ = J;
    real_type * val = pattern_.valuePtr();
    std::fill(val, val + pattern_.nonZeros(), static_cast<real_type>(0.));

    J_values_ = RealMatRM::Zero(nb_rows, pattern_.nonZeros());
    row_has_J_.assign(static_cast<size_t>(nb_rows), 0);
    stats_ = LinearSolverStats();
}

void BatchAdjoint::clear()
{
    pattern_ = Eigen::SparseMatrix<real_type>();
    J_values_ = RealMatRM();
    row_has_J_.clear();
    stats_ = LinearSolverStats();
}

void BatchAdjoint::store_row(Eigen::Index row, const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & J)
{
    if(row < 0 || row >= J_values_.rows()){
        std::ostringstream exc_;
        exc_ << "BatchAdjoint::store_row: row " << row << " is out of the [0, "
             << J_values_.rows() << "[ range this adjoint was allocated for.";
        throw std::runtime_error(exc_.str());
    }
    // The whole premise of a batch is that this never fires (see the class comment).
    // It is checked anyway, on every row: a pattern that drifted would not throw, it
    // would quietly return a gradient computed from somebody else's matrix -- and the
    // comparison costs a pass over arrays the memcpy below walks regardless.
    const bool same_shape = (J.rows() == pattern_.rows()) && (J.cols() == pattern_.cols()) &&
                            (J.nonZeros() == pattern_.nonZeros());
    const bool same_pattern = same_shape &&
        std::memcmp(J.outerIndexPtr(), pattern_.outerIndexPtr(),
                    static_cast<size_t>(pattern_.outerSize() + 1) * sizeof(Eigen::SparseMatrix<real_type>::StorageIndex)) == 0 &&
        std::memcmp(J.innerIndexPtr(), pattern_.innerIndexPtr(),
                    static_cast<size_t>(pattern_.nonZeros()) * sizeof(Eigen::SparseMatrix<real_type>::StorageIndex)) == 0;
    if(!same_pattern){
        std::ostringstream exc_;
        exc_ << "BatchAdjoint::store_row: row " << row << " has a Jacobian whose sparsity "
                "pattern differs from the one the batch was allocated with (" << J.rows()
             << " x " << J.cols() << ", " << J.nonZeros() << " nonzeros, against "
             << pattern_.rows() << " x " << pattern_.cols() << ", " << pattern_.nonZeros()
             << "). A batch algorithm must keep one pattern for every row -- a row that "
                "changes a bus's role has to do it by rewriting values inside the pattern "
                "reserved up front, never by re-deriving the pattern.";
        throw std::runtime_error(exc_.str());
    }
    std::memcpy(J_values_.row(row).data(), J.valuePtr(),
                static_cast<size_t>(pattern_.nonZeros()) * sizeof(real_type));
    row_has_J_[static_cast<size_t>(row)] = 1;
}

void BatchAdjoint::_solve_range(Eigen::Index row_begin, Eigen::Index row_end,
                                const Eigen::Ref<const RealMatRM> & xbar,
                                const std::vector<std::vector<int> > & identity_rows,
                                Eigen::Index nb_directions,
                                RealMatRM & lambda, std::vector<char> & row_ok,
                                LinearSolverStats & stats, std::exception_ptr & err) const
{
    try {
        const Eigen::Index dim = dim_J();
        std::unique_ptr<IAdjointLinearSolver> solver = make_adjoint_linear_solver();
        // this worker's own copy of the shared pattern: from here on only its values
        // are ever written, which is what lets every row after the first refactorize
        Eigen::SparseMatrix<real_type> work = pattern_;
        RealVect rhs(dim);
        bool analyzed = false;

        for(Eigen::Index i = row_begin; i < row_end; ++i){
            if(!has_row(i)) continue;   // never solved in the forward: lambda stays 0

            std::memcpy(work.valuePtr(), J_values_.row(i).data(),
                        static_cast<size_t>(pattern_.nonZeros()) * sizeof(real_type));

            ErrorType err_solver;
            if(!analyzed){
                // one symbolic analysis per worker, for the whole range it owns
                err_solver = solver->analyze(work);
                if(err_solver == ErrorType::NoError) err_solver = solver->factorize(work);
                analyzed = (err_solver == ErrorType::NoError);
            } else {
                err_solver = solver->refactorize(work);
            }
            if(err_solver != ErrorType::NoError) continue;   // row_ok stays 0

            bool all_ok = true;
            for(Eigen::Index d = 0; d < nb_directions; ++d){
                rhs = xbar.row(i).segment(d * dim, dim);
                if(solver->solve_transpose(rhs) != ErrorType::NoError){
                    all_ok = false;
                    break;
                }
                lambda.row(i).segment(d * dim, dim) = rhs;
            }
            if(!all_ok) continue;

            // Drop the multipliers of the equations this row froze to the identity.
            // Such a row of J is a single 1: it feeds no other equation (every other
            // entry of that row is zero), so the rest of lambda is exactly what it
            // would have been without the frozen equation -- but the multiplier
            // itself belongs to an equation the row does not really solve, and
            // handing it back would put a gradient on an injection the solution does
            // not depend on.
            if(!identity_rows.empty()){
                for(int eq : identity_rows[static_cast<size_t>(i)]){
                    if(eq < 0 || eq >= dim) continue;
                    for(Eigen::Index d = 0; d < nb_directions; ++d) lambda(i, d * dim + eq) = 0.;
                }
            }
            row_ok[static_cast<size_t>(i)] = 1;
        }
        stats = solver->get_linear_solver_stats();
    } catch(...) {
        err = std::current_exception();
    }
}

BatchAdjoint::RealMatRM BatchAdjoint::solve_JT(
    const Eigen::Ref<const RealMatRM> & xbar,
    const std::vector<std::vector<int> > & identity_rows,
    int nb_thread,
    std::vector<char> & row_ok) const
{
    if(!is_allocated()){
        throw std::runtime_error("BatchAdjoint::solve_JT: no Jacobian was kept. Set "
                                 "`keep_jacobian` to True before running the batch.");
    }
    const Eigen::Index rows = nb_rows();
    const Eigen::Index dim = dim_J();
    if(xbar.rows() != rows){
        std::ostringstream exc_;
        exc_ << "BatchAdjoint::solve_JT: xbar has " << xbar.rows() << " rows, but the batch "
                "has " << rows << ".";
        throw std::runtime_error(exc_.str());
    }
    if(xbar.cols() <= 0 || (xbar.cols() % dim) != 0){
        std::ostringstream exc_;
        exc_ << "BatchAdjoint::solve_JT: xbar has " << xbar.cols() << " columns, which is not "
                "a positive multiple of the Jacobian's dimension (" << dim << "). Each row "
                "holds one or more cotangents of " << dim << " coefficients, laid end to end.";
        throw std::runtime_error(exc_.str());
    }
    if(!identity_rows.empty() && static_cast<Eigen::Index>(identity_rows.size()) != rows){
        std::ostringstream exc_;
        exc_ << "BatchAdjoint::solve_JT: identity_rows has " << identity_rows.size()
             << " entries for a batch of " << rows << " rows (it must be empty, or have one "
                "entry per row).";
        throw std::runtime_error(exc_.str());
    }
    const Eigen::Index nb_directions = xbar.cols() / dim;

    RealMatRM lambda = RealMatRM::Zero(rows, xbar.cols());
    row_ok.assign(static_cast<size_t>(rows), 0);

    const int nb_th = std::min(static_cast<int>(rows), std::max(1, nb_thread));
    std::vector<LinearSolverStats> th_stats(nb_th);
    std::vector<std::exception_ptr> th_err(nb_th);

    if(nb_th <= 1){
        _solve_range(0, rows, xbar, identity_rows, nb_directions, lambda, row_ok,
                     th_stats[0], th_err[0]);
    } else {
        // contiguous ranges, exactly like the forward's row loop: disjoint rows of
        // lambda / row_ok, one solver and one working matrix per worker
        std::vector<std::thread> threads;
        threads.reserve(static_cast<size_t>(nb_th));
        const Eigen::Index per_thread = rows / nb_th;
        const Eigen::Index remainder = rows % nb_th;
        Eigen::Index begin = 0;
        for(int t = 0; t < nb_th; ++t){
            const Eigen::Index end = begin + per_thread + (t < remainder ? 1 : 0);
            threads.emplace_back([this, begin, end, &xbar, &identity_rows, nb_directions,
                                  &lambda, &row_ok, &th_stats, &th_err, t](){
                this->_solve_range(begin, end, xbar, identity_rows, nb_directions,
                                   lambda, row_ok, th_stats[t], th_err[t]);
            });
            begin = end;
        }
        for(auto & th : threads) th.join();
    }

    for(int t = 0; t < nb_th; ++t){
        if(th_err[t]) std::rethrow_exception(th_err[t]);
    }

    stats_ = LinearSolverStats();
    for(const auto & s : th_stats){
        stats_.nb_analyze += s.nb_analyze;
        stats_.nb_factorize += s.nb_factorize;
        stats_.nb_refactorize += s.nb_refactorize;
        stats_.nb_refactorize_failed += s.nb_refactorize_failed;
        stats_.nb_fallback_factorize += s.nb_fallback_factorize;
        stats_.nb_fallback_factorize_failed += s.nb_fallback_factorize_failed;
        stats_.nb_solve_transpose += s.nb_solve_transpose;
        stats_.timer_initialize_ += s.timer_initialize_;
        stats_.timer_factor_ += s.timer_factor_;
        stats_.timer_refactor_ += s.timer_refactor_;
        stats_.timer_solve_transpose_ += s.timer_solve_transpose_;
    }
    return lambda;
}

} // namespace ls2g
