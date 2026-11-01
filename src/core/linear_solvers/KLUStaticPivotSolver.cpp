// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "KLUStaticPivotSolver.hpp"

#include <algorithm>
#include <cmath>

namespace ls2g {

void KLUStaticPivotLinearSolver::set_defaults(){
    klu_defaults(&common_);
    // one block: the frozen plan does not follow BTF's blocks and off-diagonal part
    common_.btf = 0;
}

bool constrained_amd_order(int n, const int * Ap, const int * Ai,
                           const OrderingConstraints & constraints,
                           std::vector<int> & perm){
    perm.resize(n);
    if(n == 0) return true;
    double control[AMD_CONTROL];
    double info[AMD_INFO];
    amd_defaults(control);
    const int status = amd_order(n, Ap, Ai, perm.data(), control, info);
    if(status != AMD_OK && status != AMD_OK_BUT_JUMBLED) return false;
    if(constraints.empty()) return true;

    std::vector<int> pos(n);
    for(int k = 0; k < n; ++k) pos[perm[k]] = k;
    // sort key: twice the AMD position; a follow index that comes before its last lead one
    // gets the odd key right after it. Leads keep their place, so the keys of all the
    // constraints are read off the same AMD order and one pass is enough.
    // An `in_order` constraint gives all its follow indices the odd key after the last of
    // its leads and follows, and their rank breaks the tie. Its follows are never another
    // constraint's leads, so the keys stay read off the AMD order alone.
    std::vector<long> key(n);
    std::vector<int> rank(n, 0);
    for(int i = 0; i < n; ++i) key[i] = 2L * pos[i];
    for(const OrderingConstraint & c : constraints){
        int last_lead = -1;
        for(int i : c.lead){
            if(i >= 0 && i < n && pos[i] > last_lead) last_lead = pos[i];
        }
        if(last_lead < 0) continue;
        if(c.in_order){
            int last = last_lead;
            for(int i : c.follow){
                if(i >= 0 && i < n) last = std::max(last, pos[i]);
            }
            int r = 0;
            for(int i : c.follow){
                if(i < 0 || i >= n) continue;
                key[i] = 2L * last + 1;
                rank[i] = r++;
            }
            continue;
        }
        for(int i : c.follow){
            if(i >= 0 && i < n && pos[i] < last_lead) key[i] = std::max(key[i], 2L * last_lead + 1);
        }
    }
    std::stable_sort(perm.begin(), perm.end(), [&key, &rank](int a, int b){
        return key[a] < key[b] || (key[a] == key[b] && rank[a] < rank[b]);
    });
    return true;
}

ErrorType KLUStaticPivotLinearSolver::reset(){
    numeric_.reset();
    symbolic_.reset();
    common_ = klu_common();
    set_defaults();
    plan_valid_ = false;
    factors_ = Factors::None;
    ordering_constraints_.clear();
    return ErrorType::NoError;
}

ErrorType KLUStaticPivotLinearSolver::analyze(const EigenRefConstRealSpMat & J){
    // J is const here, but `klu_analyze` signature expects non-const arrays, so I const_cast
    const auto n = J.cols();
    numeric_.reset();
    symbolic_.reset();
    common_ = klu_common();
    set_defaults();
    plan_valid_ = false;
    factors_ = Factors::None;
    auto * Ap = const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.outerIndexPtr());
    auto * Ai = const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.innerIndexPtr());
    if(ordering_constraints_.empty()){
        symbolic_.reset(klu_analyze(n, Ap, Ai, &common_));
    } else {
        // symmetric: the rows the constraints move sit on the diagonal of their pivot column
        std::vector<int> perm;
        if(!constrained_amd_order(static_cast<int>(n), Ap, Ai, ordering_constraints_, perm)) return ErrorType::SolverAnalyze;
        symbolic_.reset(klu_analyze_given(n, Ap, Ai, perm.data(), perm.data(), &common_));
    }
    if(common_.status != KLU_OK) return ErrorType::SolverAnalyze;
    return ErrorType::NoError;
}

ErrorType KLUStaticPivotLinearSolver::factorize(const EigenRefConstRealSpMat & J){
    // J is const here, but `klu_factor` signature expects non-const arrays, so I const_cast
    factors_ = Factors::None;
    numeric_.reset(klu_factor(const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.outerIndexPtr()),
                              const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.innerIndexPtr()),
                              const_cast<real_type*>(J.valuePtr()),
                              symbolic_.get(), &common_));
    if(common_.status != KLU_OK) return ErrorType::SolverFactor;
    // The patterns of L and U follow from the pattern of J and the pivot sequence alone: a
    // factorization that chose the frozen pivots again (the column order only changes at
    // analyze) keeps the plan, and the next refactorization does not extract anything.
    if(plan_valid_ && !same_pivots_as_plan()) plan_valid_ = false;
    factors_ = Factors::Klu;
    return ErrorType::NoError;
}

bool KLUStaticPivotLinearSolver::same_pivots_as_plan() const {
    const int * Pnum = numeric_->Pnum;
    for(int k = 0; k < n_; ++k){
        if(Pnum[k] != P_[k]) return false;
    }
    return true;
}

bool KLUStaticPivotLinearSolver::freeze_plan(const EigenRefConstRealSpMat & J){
    if(symbolic_ == nullptr || numeric_ == nullptr) return false;
    if(symbolic_->nblocks != 1 || numeric_->nzoff != 0) return false;
    const int n = static_cast<int>(symbolic_->n);
    const int lnz = numeric_->lnz;  // the diagonals included
    const int unz = numeric_->unz;

    // only the patterns and the permutations are kept (the values are what refactorize
    // computes), but klu_extract fills a factor's pattern only along with its values
    ext_Lp_.resize(n + 1);
    ext_Li_.resize(lnz);
    ext_Lx_.resize(lnz);
    ext_Up_.resize(n + 1);
    ext_Ui_.resize(unz);
    ext_Ux_.resize(unz);
    P_.resize(n);
    Q_.resize(n);
    const int ok = klu_extract(numeric_.get(), symbolic_.get(),
                               ext_Lp_.data(), ext_Li_.data(), ext_Lx_.data(),
                               ext_Up_.data(), ext_Ui_.data(), ext_Ux_.data(),
                               nullptr, nullptr, nullptr,
                               P_.data(), Q_.data(), nullptr, nullptr,
                               &common_);
    if(ok != 1) return false;

    // L without its unit diagonal, U without its diagonal and with ascending rows
    Lp_.assign(n + 1, 0);
    Li_.clear();
    Li_.reserve(std::max(lnz - n, 0));
    Up_.assign(n + 1, 0);
    Ui_.clear();
    Ui_.reserve(std::max(unz - n, 0));
    for(int k = 0; k < n; ++k){
        for(int p = ext_Lp_[k]; p < ext_Lp_[k + 1]; ++p){
            if(ext_Li_[p] != k) Li_.push_back(ext_Li_[p]);
        }
        Lp_[k + 1] = static_cast<int>(Li_.size());
        for(int p = ext_Up_[k]; p < ext_Up_[k + 1]; ++p){
            if(ext_Ui_[p] != k) Ui_.push_back(ext_Ui_[p]);
        }
        Up_[k + 1] = static_cast<int>(Ui_.size());
        std::sort(Ui_.begin() + Up_[k], Ui_.end());
    }
    Lx_.resize(Li_.size());
    Ux_.resize(Ui_.size());
    Udiag_.resize(n);

    // where each nonzero of J lands in B = J(P, Q)
    std::vector<int> Pinv(n);
    for(int k = 0; k < n; ++k) Pinv[P_[k]] = k;
    const int nnz = static_cast<int>(J.nonZeros());
    const auto * Ai = J.innerIndexPtr();
    Brow_.resize(nnz);
    for(int p = 0; p < nnz; ++p) Brow_[p] = Pinv[Ai[p]];

    x_.assign(n, 0.);
    y_.resize(n);
    n_ = n;
    nnz_ = nnz;
    plan_valid_ = true;
    return true;
}

ErrorType KLUStaticPivotLinearSolver::refactorize(const EigenRefConstRealSpMat & J){
    if(symbolic_ == nullptr || numeric_ == nullptr) return ErrorType::SolverReFactor;
    if(!plan_valid_ && !freeze_plan(J)) return ErrorType::SolverReFactor;
    // the pattern is the one of the analysis: only its size can be checked cheaply
    if(J.cols() != n_ || J.nonZeros() != nnz_) return ErrorType::SolverReFactor;
    factors_ = Factors::None;  // until this one has gone through

    const auto * Ap = J.outerIndexPtr();
    const real_type * Ax = J.valuePtr();
    const int * Brow = Brow_.data();
    const int * Lp = Lp_.data();
    const int * Li = Li_.data();
    const int * Up = Up_.data();
    const int * Ui = Ui_.data();
    real_type * Lx = Lx_.data();
    real_type * Ux = Ux_.data();
    real_type * x = x_.data();

    // The reciprocal pivot growth of this factorization: min over the columns of
    // max |B(:, k)| / max |U(:, k)|. Computed on B itself, the matrix factorized here (KLU
    // computes it on R \ B, the matrix it factorizes), it costs a max per entry read anyway.
    const bool check_growth = rgrowth_threshold_ > 0.;
    real_type rgrowth = 1.;

    for(int k = 0; k < n_; ++k){
        // x = B(:, k)
        const int col = Q_[k];
        real_type a_max = 0.;
        for(int p = Ap[col]; p < Ap[col + 1]; ++p){
            x[Brow[p]] = Ax[p];
            a_max = std::max(a_max, std::abs(Ax[p]));
        }
        // x -= L(:, j) * x(j) for the j < k with U(j, k) != 0, ascending; U(j, k) = x(j)
        real_type u_max = 0.;
        for(int q = Up[k]; q < Up[k + 1]; ++q){
            const int j = Ui[q];
            const real_type xj = x[j];
            x[j] = 0.;
            Ux[q] = xj;
            if(xj != 0.){
                for(int r = Lp[j]; r < Lp[j + 1]; ++r) x[Li[r]] -= Lx[r] * xj;
            }
            u_max = std::max(u_max, std::abs(xj));
        }
        const real_type pivot = x[k];
        x[k] = 0.;
        Udiag_[k] = pivot;
        if(!(std::abs(pivot) > 0.) || !std::isfinite(pivot)){
            // the frozen pivot is no longer usable: leave x clean for the next call
            for(int r = Lp[k]; r < Lp[k + 1]; ++r) x[Li[r]] = 0.;
            last_rgrowth_ = 0.;
            return ErrorType::SolverReFactor;
        }
        // L(:, k) = x(rows > k) / pivot
        for(int r = Lp[k]; r < Lp[k + 1]; ++r){
            const int i = Li[r];
            Lx[r] = x[i] / pivot;
            x[i] = 0.;
        }
        u_max = std::max(u_max, std::abs(pivot));
        rgrowth = std::min(rgrowth, a_max / u_max);
    }
    last_rgrowth_ = rgrowth;
    if(check_growth && !(rgrowth >= rgrowth_threshold_)) return ErrorType::SolverReFactor;
    factors_ = Factors::Own;
    return ErrorType::NoError;
}

ErrorType KLUStaticPivotLinearSolver::solve(Eigen::Ref<RealVect> b){
    if(factors_ == Factors::Klu){
        // right after a factorize: KLU's own factors
        const int n = static_cast<int>(b.size());
        int ok = klu_solve(symbolic_.get(), numeric_.get(), n, 1, &b(0), &common_);
        if(ok != 1) return ErrorType::SolverSolve;
        return ErrorType::NoError;
    }
    if(factors_ != Factors::Own || b.size() != n_) return ErrorType::SolverSolve;

    // J x = b  <=>  L U y = c, with c(k) = b(P(k)) and x(Q(k)) = y(k)
    real_type * y = y_.data();
    for(int k = 0; k < n_; ++k) y[k] = b(P_[k]);
    for(int j = 0; j < n_; ++j){
        const real_type yj = y[j];
        if(yj != 0.){
            for(int r = Lp_[j]; r < Lp_[j + 1]; ++r) y[Li_[r]] -= Lx_[r] * yj;
        }
    }
    for(int j = n_ - 1; j >= 0; --j){
        const real_type yj = y[j] / Udiag_[j];
        y[j] = yj;
        if(yj != 0.){
            for(int q = Up_[j]; q < Up_[j + 1]; ++q) y[Ui_[q]] -= Ux_[q] * yj;
        }
    }
    for(int k = 0; k < n_; ++k) b(Q_[k]) = y[k];
    return ErrorType::NoError;
}

ErrorType KLUStaticPivotLinearSolver::solve_transpose(Eigen::Ref<RealVect> b){
    if(factors_ == Factors::Klu){
        const int n = static_cast<int>(b.size());
        int ok = klu_tsolve(symbolic_.get(), numeric_.get(), n, 1, &b(0), &common_);
        if(ok != 1) return ErrorType::SolverSolve;
        return ErrorType::NoError;
    }
    if(factors_ != Factors::Own || b.size() != n_) return ErrorType::SolverSolve;

    // J^T x = b  <=>  U^T L^T z = d, with d(k) = b(Q(k)) and x(P(k)) = z(k)
    real_type * y = y_.data();
    for(int k = 0; k < n_; ++k) y[k] = b(Q_[k]);
    for(int j = 0; j < n_; ++j){
        // row j of U^T is column j of U
        real_type s = y[j];
        for(int q = Up_[j]; q < Up_[j + 1]; ++q) s -= Ux_[q] * y[Ui_[q]];
        y[j] = s / Udiag_[j];
    }
    for(int j = n_ - 1; j >= 0; --j){
        real_type s = y[j];
        for(int r = Lp_[j]; r < Lp_[j + 1]; ++r) s -= Lx_[r] * y[Li_[r]];
        y[j] = s;
    }
    for(int k = 0; k < n_; ++k) b(P_[k]) = y[k];
    return ErrorType::NoError;
}

} // namespace ls2g
