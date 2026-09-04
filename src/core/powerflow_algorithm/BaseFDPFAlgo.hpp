// Copyright (c) 2023-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASEFDPFALGO_H
#define BASEFDPFALGO_H

#include "BaseAlgo.hpp"
#include "linear_solvers/LinearSolverStats.hpp"

namespace ls2g {

/**
Base class for Fast Decoupled Powerflow based solver
**/
template<class LinearSolver, FDPFMethod XB_BX>
class BaseFDPFAlgo final: public BaseAlgo
{
    public:

        BaseFDPFAlgo() noexcept :BaseAlgo(true), need_factorize_(true) {}
        ~BaseFDPFAlgo() noexcept override = default;


        static constexpr bool IS_FDPF = true;
        bool is_fdpf() const noexcept override { return IS_FDPF; }

        bool compute_pf(const EigenRefConstCplxSpMat     & Ybus,
                        const Eigen::Ref<const CplxVect> & V,
                        const Eigen::Ref<const CplxVect> & Sbus,
                        const Eigen::Ref<const IntVect>  & slack_ids,
                        const Eigen::Ref<const RealVect> & slack_weights,
                        const Eigen::Ref<const IntVect>  & pv,
                        const Eigen::Ref<const IntVect>  & pq,
                        int                              max_iter,
                        real_type                        tol
                        ) override;  // requires a gridmodel !

        void reset() override
        {   
            BaseAlgo::reset();
            // solution of the problem
            Bp_ = Eigen::SparseMatrix<real_type> ();  // the B prime matrix (size n_pvpq)
            Bpp_ = Eigen::SparseMatrix<real_type>();  // the B double prime matrix  (size n_pq)
            grid_Bp_ = Eigen::SparseMatrix<real_type> ();  // the B prime matrix (size n_pvpq)
            grid_Bpp_ = Eigen::SparseMatrix<real_type>();  // the B double prime matrix  (size n_pq)
            p_ = RealVect();
            q_ = RealVect();
            mis_ = CplxVect();
            ybus_v_ = CplxVect();
            need_factorize_ = true;

            // reset linear solvers
            ErrorType reset_status = _linear_solver_Bp.reset();
            if(reset_status != ErrorType::NoError) err_ = reset_status;
            reset_status = _linear_solver_Bpp.reset();
            if((reset_status != ErrorType::NoError) && (err_ != ErrorType::NotInitError)) err_ = reset_status;
        }

        // bare (non-Ref) return type: these are bound directly to python
        // (binding_solvers.cpp), and pybind11's Eigen support has no type_caster for
        // Eigen::Ref<const SparseMatrix<...>> (only for a concrete, owning
        // SparseMatrix -- see TODO in CHANGELOG.rst), same reason as get_J_python.
        Eigen::SparseMatrix<real_type> debug_get_Bp_python()  const { return Bp_;}
        Eigen::SparseMatrix<real_type> debug_get_Bpp_python() const { return Bpp_;}

        // Two linear solvers (B' and B''): reported separately rather than merged, so
        // no information is lost (unlike timer_initialize_ below, which historically
        // combines both for backward compatibility with get_timers_jacobian()).
        LinearSolverStats get_linear_solver_stats_bp() const {
            return detail::get_stats_impl(_linear_solver_Bp, 0);
        }
        LinearSolverStats get_linear_solver_stats_bpp() const {
            return detail::get_stats_impl(_linear_solver_Bpp, 0);
        }

        TimerJac get_timers_jacobian() const override
        {
            const LinearSolverStats bp  = detail::get_stats_impl(_linear_solver_Bp, 0);
            const LinearSolverStats bpp = detail::get_stats_impl(_linear_solver_Bpp, 0);
            TimerJac res;
            res.timer_Fx_         = timer_Fx_;
            res.timer_solve_      = bp.timer_solve_ + bpp.timer_solve_;
            // FDPF never refactorizes (B'/B'' are fixed matrices, factorized once), and
            // historically combined analyze+factor into a single timer_initialize_ for
            // both solvers -- preserved here for backward compatibility.
            res.timer_initialize_ = bp.timer_initialize_ + bp.timer_factor_
                                   + bpp.timer_initialize_ + bpp.timer_factor_;
            res.timer_check_      = timer_check_;
            res.timer_total_nr_   = timer_total_nr_;
            return res;
        }

    protected:
        void reset_timer() override {
            BaseAlgo::reset_timer();
            detail::reset_stats_timers_impl(_linear_solver_Bp, 0);
            detail::reset_stats_timers_impl(_linear_solver_Bpp, 0);
        }

        // (V * conj(Ybus * V) - Sbus + slack_absorbed * slack_weights), written
        // into the persistent mis_ buffer rather than returned by value: this
        // runs once or twice per FDPF iteration, and returning an nb_bus
        // complex vector meant heap-allocating and freeing one every time.
        // Ybus * V lands in its own buffer first (`noalias`): a sparse * dense
        // product left inside the coefficient-wise expression is one more
        // temporary Eigen has to allocate to evaluate it.
        void evaluate_mismatch_into(const EigenRefConstCplxSpMat     &  Ybus,
                                    const Eigen::Ref<const CplxVect> & V,
                                    const Eigen::Ref<const CplxVect> & Sbus,
                                    size_t /*slack_id*/,  // id of the ref slack bus
                                    real_type slack_absorbed,
                                    const Eigen::Ref<const RealVect> & slack_weights)
        {
            ybus_v_.noalias() = Ybus * V;
            mis_ = V.array() * ybus_v_.array().conjugate() - Sbus.array() + slack_absorbed * slack_weights.array();
        }

        // value-returning form, kept for out-of-tree derived algorithms
        CplxVect evaluate_mismatch(const EigenRefConstCplxSpMat     &  Ybus,
                                   const Eigen::Ref<const CplxVect> & V,
                                   const Eigen::Ref<const CplxVect> & Sbus,
                                   size_t slack_id,  // id of the ref slack bus
                                   real_type slack_absorbed,
                                   const Eigen::Ref<const RealVect> & slack_weights)
        {
            evaluate_mismatch_into(Ybus, V, Sbus, slack_id, slack_absorbed, slack_weights);
            return mis_;
        }
        
        void fillBp_Bpp(
            Eigen::SparseMatrix<real_type> & Bp, 
            Eigen::SparseMatrix<real_type> & Bpp) const;  // defined in Solvers.cpp !

        void initialize() {
            err_ = ErrorType::NoError; // reset error message
            // analyze (structure) + factorize (values) for Bp solver
            ErrorType init_status = _linear_solver_Bp.analyze(Bp_);
            if(init_status == ErrorType::NoError)
                init_status = _linear_solver_Bp.factorize(Bp_);
            if(init_status != ErrorType::NoError){
                _linear_solver_Bp.reset();
                _linear_solver_Bpp.reset();
                err_ = init_status;
                need_factorize_ = true;
                return;
            }
            // analyze (structure) + factorize (values) for Bpp solver (if Bp succeeded)
            init_status = _linear_solver_Bpp.analyze(Bpp_);
            if(init_status == ErrorType::NoError)
                init_status = _linear_solver_Bpp.factorize(Bpp_);
            if(init_status != ErrorType::NoError){
                _linear_solver_Bp.reset();
                _linear_solver_Bpp.reset();
                err_ = init_status;
                need_factorize_ = true;
                return;
            }

            // everything went well
            need_factorize_ = false;
        }

        void solve(LinearSolver& linear_solver,
                   Eigen::Ref<RealVect> b){
            const ErrorType solve_status = linear_solver.solve(b);
            if(solve_status != ErrorType::NoError){
                // std::cout << "solve error: " << solve_status << std::endl;
                err_ = solve_status;
            }
        }

        bool has_converged(
            const Eigen::Ref<const CplxVect >       & tmp_va,
            const EigenRefConstCplxSpMat            & Ybus,
            const Eigen::Ref<const CplxVect>        & Sbus,
            size_t                                  slack_bus_id,
            real_type                               & slack_absorbed,
            const Eigen::Ref<const RealVect>        & slack_weights,
            const Eigen::Ref<const Eigen::VectorXi> & pvpq,
            const Eigen::Ref<const Eigen::VectorXi> & pq,
            real_type                               tol)
        {
            /**
            It is suppose to implement the python code bellow (so it updates p_, q_, v_, va_ and vm_):

            .. code-block:: python

                V = Vm * exp(1j * Va)

                ## evalute mismatch
                mis = (V * conj(Ybus * V) - Sbus) / Vm
                P = mis[pvpq].real
                Q = mis[pq].imag

                ## check tolerance
                normP = linalg.norm(P, Inf)
                normQ = linalg.norm(Q, Inf)
                if verbose > 1:
                    sys.stdout.write('\n  Q  %3d   %10.3e   %10.3e' % (i, normP, normQ))
                if normP < tol and normQ < tol:
                    converged = 1
                    if verbose:
                        sys.stdout.write('\nFast-decoupled power flow converged in %d '
                            'P-iterations and %d Q-iterations.\n' % (i, i))
                    break
            **/

            // V = Vm * exp(1j * Va) : 
            V_ = Vm_.array() * tmp_va.array();
            // ... then repair a magnitude the Q iteration drove past zero, because
            // the division by Vm_ two lines below would otherwise flip that bus's P
            // and Q. Guarded the way the Newton-Raphson guards the same repair (see
            // NRSystem::apply_step): a converging trajectory never goes negative, so
            // the ordinary solve pays one reduction instead of two full-length
            // passes, twice per iteration. V_ is built above and is unchanged either
            // way. The angle wrap that used to run here with it has moved to the end
            // of compute_pf -- see wrap_va for why once is enough.
            if (Vm_.minCoeff() < my_zero_) fix_negative_vm(Vm_, Va_);

            evaluate_mismatch_into(Ybus, V_, Sbus, slack_bus_id, slack_absorbed, slack_weights);  // mis_ = V * conj(Ybus * V) - Sbus
            mis_.array() /= Vm_.array();  // mis = (V * conj(Ybus * V) - Sbus) / Vm (do not forget the / Vm !)

            bool tmp = mis_.allFinite();
            if(!tmp){
                err_ = ErrorType::InifiniteValue;
                return false; // divergence due to Nans
            }
            p_ = mis_(pvpq).real();  // P = mis[pvpq].real
            q_ = mis_(pq).imag();  // Q = mis[pq].imag
            return _check_for_convergence(p_, q_, tol);
        }

        void fill_sparse_matrices(
            const EigenRefConstRealSpMat & grid_Bp,
            const EigenRefConstRealSpMat & grid_Bpp,
            const std::vector<int>       & pvpq_inv,
            const std::vector<int>       & pq_inv,
            size_t                       n_pvpq,
            size_t                       n_pq);

        void aux_fill_sparse_matrices(
            const EigenRefConstRealSpMat   & grid_Bp_Bpp,
            const std::vector<int>         & ind_inv,
            size_t                         mat_dim,
            Eigen::SparseMatrix<real_type> & res);

    protected:
        // use 2 linear solvers
        LinearSolver _linear_solver_Bp;
        LinearSolver _linear_solver_Bpp;

        // solution of the problem
        Eigen::SparseMatrix<real_type> grid_Bp_;
        Eigen::SparseMatrix<real_type> grid_Bpp_;
        Eigen::SparseMatrix<real_type> Bp_;  // the B prime matrix (size n_pvpq)
        Eigen::SparseMatrix<real_type> Bpp_;  // the B double prime matrix  (size n_pq)
        RealVect p_;  // (size n_pvpq)
        RealVect q_;  // (size n_pq)
        // per-iteration scratch, kept across iterations (see evaluate_mismatch_into)
        CplxVect mis_;     // per-bus complex power mismatch (size nb_bus)
        CplxVect ybus_v_;  // Ybus * V (size nb_bus)
        bool need_factorize_;

    private:
        // no copy allowed
        BaseFDPFAlgo(const BaseFDPFAlgo&) = delete;
        BaseFDPFAlgo(BaseFDPFAlgo&&) = delete;
        BaseFDPFAlgo & operator=(BaseFDPFAlgo&&) = delete;
        BaseFDPFAlgo & operator=(const BaseFDPFAlgo&) = delete;
};

#include "BaseFDPFAlgo.tpp"


} // namespace ls2g

#endif // BASEFDPFALGO_H