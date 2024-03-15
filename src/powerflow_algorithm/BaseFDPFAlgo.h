// Copyright (c) 2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASEFDPFALGO_H
#define BASEFDPFALGO_H

#include "BaseAlgo.h"

/**
Base class for Fast Decoupled Powerflow based solver
**/
template<class LinearSolver, FDPFMethod XB_BX>
class BaseFDPFAlgo: public BaseAlgo
{
    public:
        BaseFDPFAlgo():BaseAlgo(true), need_factorize_(true) {}

        virtual
        bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                        CplxVect & V,
                        const CplxVect & Sbus,
                        const Eigen::VectorXi & slack_ids,
                        const RealVect & slack_weights,
                        const Eigen::VectorXi & pv,
                        const Eigen::VectorXi & pq,
                        int max_iter,
                        real_type tol
                        ) ;  // requires a gridmodel !

        // bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
        //                 const Eigen::SparseMatrix<cplx_type> & Bp,
        //                 const Eigen::SparseMatrix<cplx_type> & Bpp,
        //                 CplxVect & V,
        //                 const CplxVect & Sbus,
        //                 const Eigen::VectorXi & slack_ids,
        //                 const RealVect & slack_weights,
        //                 const Eigen::VectorXi & pv,
        //                 const Eigen::VectorXi & pq,
        //                 int max_iter,
        //                 real_type tol
        //                 ) ;  // TODO add Bp and Bpp as argument for use in python directly !

        virtual void reset()
        {   
            BaseAlgo::reset();
            // solution of the problem
            Bp_ = Eigen::SparseMatrix<real_type> ();  // the B prime matrix (size n_pvpq)
            Bpp_ = Eigen::SparseMatrix<real_type>();  // the B double prime matrix  (size n_pq)
            grid_Bp_ = Eigen::SparseMatrix<real_type> ();  // the B prime matrix (size n_pvpq)
            grid_Bpp_ = Eigen::SparseMatrix<real_type>();  // the B double prime matrix  (size n_pq)
            p_ = RealVect();
            q_ = RealVect();
            need_factorize_ = true;

            // reset linear solvers
            ErrorType reset_status = _linear_solver_Bp.reset();
            if(reset_status != ErrorType::NoError) err_ = reset_status;
            reset_status = _linear_solver_Bpp.reset();
            if((reset_status != ErrorType::NoError) && (err_ != ErrorType::NotInitError)) err_ = reset_status;
        }

        Eigen::SparseMatrix<real_type> debug_get_Bp_python() { return Bp_;}
        Eigen::SparseMatrix<real_type> debug_get_Bpp_python() { return Bpp_;}

    protected:
        virtual void reset_timer(){
            BaseAlgo::reset_timer();
        }

        CplxVect evaluate_mismatch(const Eigen::SparseMatrix<cplx_type> &  Ybus,
                                   const CplxVect & V,
                                   const CplxVect & Sbus,
                                   Eigen::Index slack_id,  // id of the ref slack bus
                                   real_type slack_absorbed,
                                   const RealVect & slack_weights)
        {
            CplxVect tmp = Ybus * V;  // this is a vector
            tmp = tmp.array().conjugate();  // i take the conjugate
            auto mis = V.array() * tmp.array() - Sbus.array() + slack_absorbed * slack_weights.array();
            return mis;
        }
        
        void fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, Eigen::SparseMatrix<real_type> & Bpp) const;  // defined in Solvers.cpp !

        virtual
        void initialize(){
            auto timer = CustTimer();
            err_ = ErrorType::NoError; // reset error message
            // init Bp solver
            ErrorType init_status = _linear_solver_Bp.initialize(Bp_);
            if(init_status != ErrorType::NoError){
                _linear_solver_Bp.reset();
                _linear_solver_Bpp.reset();
                err_ = init_status;
                need_factorize_ = true;
                timer_initialize_ += timer.duration();
                return;
            }
            // init Bpp solver (if Bp is sucesfull of course)
            init_status = _linear_solver_Bpp.initialize(Bpp_);
            if(init_status != ErrorType::NoError){
                _linear_solver_Bp.reset();
                _linear_solver_Bpp.reset();
                err_ = init_status;
                need_factorize_ = true;
                timer_initialize_ += timer.duration();
                return;
            }

            // everything went well
            need_factorize_ = false;
            timer_initialize_ += timer.duration();
        }

        virtual
        void solve(LinearSolver& linear_solver,
                   Eigen::SparseMatrix<real_type>& mat,
                   RealVect & b,
                   bool has_just_been_inialized){
            auto timer = CustTimer();
            // const ErrorType solve_status = linear_solver.solve(mat, b, has_just_been_inialized);
            const ErrorType solve_status = linear_solver.solve(mat, b, true);  // true because i don't need to refactorize the matrix
            if(solve_status != ErrorType::NoError){
                // std::cout << "solve error: " << solve_status << std::endl;
                err_ = solve_status;
            }
            timer_solve_ += timer.duration();
        }

        bool has_converged(const Eigen::Ref<const CplxVect > & tmp_va,
                           const Eigen::SparseMatrix<cplx_type> & Ybus,
                           const CplxVect & Sbus,
                           Eigen::Index slack_bus_id,
                           real_type & slack_absorbed,
                           const RealVect & slack_weights,
                           const Eigen::Ref<const Eigen::VectorXi> & pvpq,
                           const Eigen::Ref<const Eigen::VectorXi> & pq,
                           real_type tol)
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
            Vm_ = V_.array().abs();
            Va_ = V_.array().arg();
            auto mis = evaluate_mismatch(Ybus, V_, Sbus, slack_bus_id, slack_absorbed, slack_weights);  // mis = (V * conj(Ybus * V) - Sbus) / Vm
            mis.array() /= Vm_.array();  // mis = (V * conj(Ybus * V) - Sbus) / Vm (do not forget the / Vm !)

            bool tmp = mis.allFinite();
            if(!tmp){
                err_ = ErrorType::InifiniteValue;
                return false; // divergence due to Nans
            }
            p_ = mis(pvpq).real();  // P = mis[pvpq].real
            q_ = mis(pq).imag();  // Q = mis[pq].imag
            return _check_for_convergence(p_, q_, tol);
        }

        void fill_sparse_matrices(const Eigen::SparseMatrix<real_type> & grid_Bp,
                                  const Eigen::SparseMatrix<real_type> & grid_Bpp,
                                  const std::vector<int> & pvpq_inv,
                                  const std::vector<int> & pq_inv,
                                  Eigen::Index n_pvpq,
                                  Eigen::Index n_pq);

        void aux_fill_sparse_matrices(const Eigen::SparseMatrix<real_type> & grid_Bp_Bpp,
                                      const std::vector<int> & ind_inv,
                                      Eigen::Index mat_dim,
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
        bool need_factorize_;

        // timers
        double timer_initialize_;
        // double timer_dSbus_;
        // double timer_fillJ_;

    private:
        // no copy allowed
        BaseFDPFAlgo( const BaseFDPFAlgo & ) =delete ;
        BaseFDPFAlgo & operator=( const BaseFDPFAlgo & ) =delete ;
};

#include "BaseFDPFAlgo.tpp"

#endif // BASEFDPFALGO_H
