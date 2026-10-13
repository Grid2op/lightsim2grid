// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifdef KLU_SOLVER_AVAILABLE
#ifndef KLU_STATIC_PIVOT_SOLVER_H
#define KLU_STATIC_PIVOT_SOLVER_H
#include <memory>
#include <vector>

#include "Utils.hpp"
#include "linear_solvers/OrderingConstraints.hpp"

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

namespace ls2g {

extern "C" {
    #include "cs.h"
    #include "klu.h"
}

/**
KLU with a static-pivot refactorization.

`factorize` is KLU's own (`klu_factor`: it chooses the pivots). `refactorize` is not
`klu_refactor`: it keeps the row and column permutations and the nonzero patterns of L and U
of the last `klu_factor`, frozen in plain compressed-column arrays, and recomputes the values
with the left-looking column algorithm of Gilbert and Peierls on that frozen order, with no
pivot search, no symbolic step and no row scaling. For column k of B = A(P, Q): scatter
B(:, k) into a dense work vector x; for every j < k with U(j, k) != 0, in ascending j (a
topological order, L being lower triangular), x -= L(:, j) * x(j); then U(:, k) = x(rows <=
k) and L(:, k) = x(rows > k) / x(k). That is the same numerical work as `klu_refactor`,
without the scaling of the rows and the packed storage it goes through on every call.

A frozen pivot order can turn bad when the values change. A zero or non-finite pivot always
rejects the refactorization; above a threshold (`set_rgrowth_threshold`, 0 by default), so
does a reciprocal pivot growth below it, computed on the new factors as the KLU thesis
defines it: min over the columns j of max |B(:, j)| / max |U(:, j)|. B is the matrix
factorized here; KLU computes the same on R \ B, the matrix it factorizes after scaling its
rows. A rejection returns ErrorType::SolverReFactor, and the fallback of
LinearSolverPolicy (set_refactor_fallback) runs a fresh `klu_factor`, which chooses new pivots.

`analyze` turns KLU's block-triangular preordering off, so that the factorization is one
block (a connected power-flow Jacobian is irreducible anyway) and the frozen plan needs no
block structure. Given ordering constraints (`set_ordering_constraints`, which the Newton
system hands over before each analysis), it orders with AMD then moves each constrained index
just after the last index it must follow (`constrained_amd_order`), instead of KLU's own AMD:
the rows the outer loops switch by value then keep a usable pivot in both their forms.

References:

- J. R. Gilbert, T. Peierls, "Sparse Partial Pivoting in Time Proportional to Arithmetic
  Operations", SIAM Journal on Scientific and Statistical Computing 9(5):862-874, 1988.
  https://doi.org/10.1137/0909058 -- the left-looking column algorithm.
- T. A. Davis, E. Palamadai Natarajan, "Algorithm 907: KLU, A Direct Sparse Solver for
  Circuit Simulation Problems", ACM Transactions on Mathematical Software 37(3), article 36,
  2010. https://doi.org/10.1145/1824801.1824814 -- the ordering and the first factorization.
- E. Palamadai Natarajan, "KLU -- A High Performance Sparse Linear Solver for Circuit
  Simulation Problems", MS thesis, University of Florida, 2005.
  https://ufdcimages.uflib.ufl.edu/UF/E0/01/17/21/00001/palamadai_e.pdf -- refactorization
  on the same permutations and nonzero pattern of L and U, and the reciprocal pivot growth.
- X. Chen, L. Xia, Y. Wang, H. Yang, "Sparsity-Oriented Sparse Solver Design for Circuit
  Simulation", DATE 2016, pp. 1580-1585. https://doi.org/10.3850/9783981537079_0839 --
  refactorization without pivoting, column k depending on column j iff U(j, k) != 0.
- X. Chen, Y. Wang, H. Yang, "NICSLU: An Adaptive Sparse Matrix Solver for Parallel Circuit
  Simulation", IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems
  32(2):261-274, 2013. https://doi.org/10.1109/TCAD.2012.2217964 -- refactorization without
  pivoting in a Gilbert-Peierls solver.
- B. Bruned, J. Mahseredjian et al., "Sparse solver application for parallel real-time
  electromagnetic transient simulations", Electric Power Systems Research 223:109585, 2023.
  https://doi.org/10.1016/j.epsr.2023.109585 -- a refactorization stopped on an invalid
  pivot, followed by a full factorization.
**/
class LS2G_API KLUStaticPivotLinearSolver final
{
    public:
        KLUStaticPivotLinearSolver() noexcept :
            common_(),
            symbolic_(nullptr, SymbolicDeleter{&common_}),
            numeric_(nullptr, NumericDeleter{&common_})
            { set_defaults(); }

        ~KLUStaticPivotLinearSolver() noexcept = default;

        // public api
        ErrorType reset();
        ErrorType analyze(const EigenRefConstRealSpMat & J);      // KLU ordering, no BTF
        ErrorType factorize(const EigenRefConstRealSpMat & J);    // klu_factor: new pivots
        ErrorType refactorize(const EigenRefConstRealSpMat & J);  // static pivots, see above
        ErrorType solve(Eigen::Ref<RealVect> b);
        ErrorType solve_transpose(Eigen::Ref<RealVect> b);  // J^T x = b, out of the factorization of J

        // can this linear solver solve problem where RHS is a matrix
        static constexpr bool CAN_SOLVE_MAT = false;

        // can this linear solver solve J^T x = b out of the factorization of J
        static constexpr bool CAN_SOLVE_TRANSPOSE = true;

        // Above 0, a refactorization whose reciprocal pivot growth falls below this threshold
        // is rejected, so that the caller factorizes again (choosing new pivots); same meaning
        // as KLULinearSolver::set_rgrowth_threshold. 0 (the default): only a zero or
        // non-finite pivot rejects it.
        void set_rgrowth_threshold(real_type val) noexcept { rgrowth_threshold_ = val; }
        real_type rgrowth_threshold() const noexcept { return rgrowth_threshold_; }

        // the reciprocal pivot growth of the last refactorization (1 when it was not computed)
        real_type last_rgrowth() const noexcept { return last_rgrowth_; }

        // Indices to eliminate after others, for the next analyze() (see OrderingConstraints).
        // Empty (the default): KLU's own AMD ordering.
        void set_ordering_constraints(OrderingConstraints constraints) { ordering_constraints_ = std::move(constraints); }
        const OrderingConstraints & ordering_constraints() const noexcept { return ordering_constraints_; }

    private:
        // see KLULinearSolver: the deleters need `common_`, declared before the handles
        struct SymbolicDeleter {
            klu_common* common_ = nullptr;
            void operator()(klu_symbolic* ptr) const noexcept {
                if(ptr != nullptr){
                    klu_symbolic* tmp = ptr;
                    klu_free_symbolic(&tmp, common_);
                }
            }
        };
        struct NumericDeleter {
            klu_common* common_ = nullptr;
            void operator()(klu_numeric* ptr) const noexcept {
                if(ptr != nullptr){
                    klu_numeric* tmp = ptr;
                    klu_free_numeric(&tmp, common_);
                }
            }
        };

        // which factors a solve uses
        enum class Factors { None, Klu, Own };

        void set_defaults();
        // freezes the permutations and the patterns of L and U of the current KLU numeric
        // factorization; false if it cannot (more than one block, an extract that fails)
        bool freeze_plan(const EigenRefConstRealSpMat & J);
        // whether the current KLU numeric factorization pivoted on the frozen rows
        bool same_pivots_as_plan() const;

        klu_common common_;
        std::unique_ptr<klu_symbolic, SymbolicDeleter> symbolic_;
        std::unique_ptr<klu_numeric, NumericDeleter> numeric_;
        real_type rgrowth_threshold_ = 0.;
        real_type last_rgrowth_ = 1.;
        Factors factors_ = Factors::None;
        OrderingConstraints ordering_constraints_;

        // the frozen plan, valid while `plan_valid_`
        bool plan_valid_ = false;
        int n_ = 0;
        int nnz_ = 0;
        std::vector<int> P_;      // row k of B is row P_[k] of A
        std::vector<int> Q_;      // column k of B is column Q_[k] of A
        std::vector<int> Brow_;   // for each nonzero of A (its storage order): its row in B
        // strictly lower part of L (unit diagonal not stored), compressed by column
        std::vector<int> Lp_, Li_;
        std::vector<real_type> Lx_;
        // strictly upper part of U, compressed by column, rows ascending; diagonal apart
        std::vector<int> Up_, Ui_;
        std::vector<real_type> Ux_, Udiag_;

        // work vectors: x_ is all zeros between two refactorizations
        std::vector<real_type> x_;
        std::vector<real_type> y_;
        // klu_extract's output, kept to not reallocate at every freeze
        std::vector<int> ext_Lp_, ext_Li_, ext_Up_, ext_Ui_;
        std::vector<real_type> ext_Lx_, ext_Ux_;

        // no copy allowed
        KLUStaticPivotLinearSolver(const KLUStaticPivotLinearSolver&) = delete;
        KLUStaticPivotLinearSolver(KLUStaticPivotLinearSolver&&) = delete;
        KLUStaticPivotLinearSolver & operator=(KLUStaticPivotLinearSolver&&) = delete;
        KLUStaticPivotLinearSolver & operator=(const KLUStaticPivotLinearSolver&) = delete;
};

/**
A fill-reducing ordering of the n-by-n pattern (Ap, Ai) that honours `constraints`: AMD on the
pattern of A + A^T, then every `follow` index placed after the last of its `lead` ones, by a
stable sort on the position (an index only ever moves down, to just after its last lead, and
keeps its relative order with what it passes). Moving the follow indices down rather than the
lead ones up is what keeps the fill local: a lead index moved up drags its whole neighbourhood
forward. Returns false if AMD fails; `perm` (size n) is the elimination order, perm[k] being
the index eliminated k-th.
**/
LS2G_API bool constrained_amd_order(int n, const int * Ap, const int * Ai,
                                    const OrderingConstraints & constraints,
                                    std::vector<int> & perm);

} // namespace ls2g

#endif // KLU_STATIC_PIVOT_SOLVER_H

#elif defined(_READ_THE_DOCS)
#ifndef KLU_STATIC_PIVOT_SOLVER_H
#define KLU_STATIC_PIVOT_SOLVER_H
#include "Utils.hpp"

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

namespace ls2g {

/**
Doc-only stand-in for KLUStaticPivotLinearSolver, see the KLULinearSolver one: same name and
public interface, every method a no-op.
**/
class LS2G_API KLUStaticPivotLinearSolver final
{
    public:
        KLUStaticPivotLinearSolver() noexcept = default;
        ~KLUStaticPivotLinearSolver() noexcept = default;

        ErrorType reset() { return ErrorType::NoError; }
        ErrorType analyze(const EigenRefConstRealSpMat & /*J*/) { return ErrorType::NoError; }
        ErrorType factorize(const EigenRefConstRealSpMat & /*J*/) { return ErrorType::NoError; }
        ErrorType refactorize(const EigenRefConstRealSpMat & /*J*/) { return ErrorType::NoError; }
        ErrorType solve(Eigen::Ref<RealVect> /*b*/) { return ErrorType::NoError; }
        ErrorType solve_transpose(Eigen::Ref<RealVect> /*b*/) { return ErrorType::NoError; }

        static constexpr bool CAN_SOLVE_MAT = false;
        static constexpr bool CAN_SOLVE_TRANSPOSE = true;

    private:
        KLUStaticPivotLinearSolver(const KLUStaticPivotLinearSolver&) = delete;
        KLUStaticPivotLinearSolver(KLUStaticPivotLinearSolver&&) = delete;
        KLUStaticPivotLinearSolver & operator=(KLUStaticPivotLinearSolver&&) = delete;
        KLUStaticPivotLinearSolver & operator=(const KLUStaticPivotLinearSolver&) = delete;
};

} // namespace ls2g
#endif // KLU_STATIC_PIVOT_SOLVER_H
#endif // KLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)
