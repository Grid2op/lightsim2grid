// Copyright (c) 2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASEFDPFSOLVER_H
#define BASEFDPFSOLVER_H

#include "BaseSolver.h"

/**
Base class for Fast Decoupled Powerflow based solver
**/
template<class LinearSolver, FDPFMethod XB_BX>
class BaseFDPFSolver : public BaseSolver
{
    public:
        BaseFDPFSolver():BaseSolver(true), need_factorize_(true) {}

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
            BaseSolver::reset();
            // solution of the problem
            Bp_ = Eigen::SparseMatrix<real_type> ();  // the B prime matrix (size n_pvpq)
            Bpp_ = Eigen::SparseMatrix<real_type>();  // the B double prime matrix  (size n_pq)
            p_ = RealVect();
            q_ = RealVect();
            need_factorize_ = true;

            // reset linear solvers
            ErrorType reset_status = _linear_solver_Bp.reset();
            if(reset_status != ErrorType::NoError) err_ = reset_status;
            reset_status = _linear_solver_Bpp.reset();
            if((reset_status != ErrorType::NoError) && (err_ != ErrorType::NotInitError)) err_ = reset_status;
        }

    protected:
        virtual void reset_timer(){
            BaseSolver::reset_timer();
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
        
        void fillBp(Eigen::SparseMatrix<real_type> & res) const;  // defined in Solvers.cpp !
        void fillBpp(Eigen::SparseMatrix<real_type> & res) const;  // defined in Solvers.cpp !

        // TODO !!!
        // virtual
        // void initialize(){
        //     auto timer = CustTimer();
        //     n_ = static_cast<int>(J_.cols()); // should be equal to J_.nrows()
        //     err_ = ErrorType::NoError; // reset error message
        //     const ErrorType init_status = _linear_solver.initialize(J_);
        //     if(init_status != ErrorType::NoError){
        //         // std::cout << "init_ok " << init_ok << std::endl;
        //         err_ = init_status;
        //     }
        //     need_factorize_ = false;
        //     timer_initialize_ += timer.duration();
        // }

        virtual
        void solve(LinearSolver& linear_solver,
                   Eigen::SparseMatrix<real_type>& mat,
                   RealVect & b,
                   bool has_just_been_inialized){
            auto timer = CustTimer();
            const ErrorType solve_status = linear_solver.solve(mat, b, has_just_been_inialized);
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
            bool tmp = mis.allFinite();
            if(!tmp){
                err_ = ErrorType::InifiniteValue;
                return false; // divergence due to Nans
            }
            p_ = mis.real()(pvpq);  // P = mis[pvpq].real
            q_ = mis.imag()(pq);  // Q = mis[pq].imag
            return _check_for_convergence(p_, q_, tol);
        }

        void _dSbus_dV(const Eigen::Ref<const Eigen::SparseMatrix<cplx_type> > & Ybus,
                       const Eigen::Ref<const CplxVect > & V);

        void _get_values_J(int & nb_obj_this_col,
                           std::vector<Eigen::Index> & inner_index,
                           std::vector<real_type> & values,
                           const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & mat,  // ex. dS_dVa_r
                           const std::vector<int> & index_row_inv, // ex. pvpq_inv
                           const Eigen::VectorXi & index_col, // ex. pvpq
                           Eigen::Index col_id,
                           Eigen::Index row_lag,  // 0 for J11 for example, n_pvpq for J12
                           Eigen::Index col_lag
                           );
        void _get_values_J(int & nb_obj_this_col,
                           std::vector<Eigen::Index> & inner_index,
                           std::vector<real_type> & values,
                           const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & mat,  // ex. dS_dVa_r
                           const std::vector<int> & index_row_inv, // ex. pvpq_inv
                           Eigen::Index col_id_mat, // ex. pvpq(col_id)
                           Eigen::Index row_lag,  // 0 for J11 for example, n_pvpq for J12
                           Eigen::Index col_lag  // to remove the ref slack bus from this
                           );

        void fill_jacobian_matrix(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                  const CplxVect & V,
                                  Eigen::Index slack_bus_id,
                                  const RealVect & slack_weights,
                                  const Eigen::VectorXi & pq,
                                  const Eigen::VectorXi & pvpq,
                                  const std::vector<int> & pq_inv,
                                  const std::vector<int> & pvpq_inv
                                  );
        void fill_jacobian_matrix_kown_sparsity_pattern(
                 Eigen::Index slack_bus_id,
                 const Eigen::VectorXi & pq,
                 const Eigen::VectorXi & pvpq
                 );
        void fill_jacobian_matrix_unkown_sparsity_pattern(
                 const Eigen::SparseMatrix<cplx_type> & Ybus,
                 const CplxVect & V,
                 Eigen::Index slack_bus_id,
                 const RealVect & slack_weights,
                 const Eigen::VectorXi & pq,
                 const Eigen::VectorXi & pvpq,
                 const std::vector<int> & pq_inv,
                 const std::vector<int> & pvpq_inv
                 );

        void fill_value_map(Eigen::Index slack_bus_id,
                            const Eigen::VectorXi & pq,
                            const Eigen::VectorXi & pvpq);

    protected:
        // used linear solver
        LinearSolver _linear_solver_Bp;
        LinearSolver _linear_solver_Bpp;

        // solution of the problem
        Eigen::SparseMatrix<real_type> Bp_;  // the B prime matrix (size n_pvpq)
        Eigen::SparseMatrix<real_type> Bpp_;  // the B double prime matrix  (size n_pq)
        RealVect p_;  // (size n_pvpq)
        RealVect q_;  // (size n_pq)
        bool need_factorize_;

        // to store the mapping from the element of J_ in dS_dVm_ and dS_dVa_
        // it does not own any memory at all !
        // std::vector<cplx_type*> value_map_;

        // timers
        // double timer_initialize_;
        // double timer_dSbus_;
        // double timer_fillJ_;

    private:
        // no copy allowed
        BaseFDPFSolver( const BaseFDPFSolver & ) =delete ;
        BaseFDPFSolver & operator=( const BaseFDPFSolver & ) =delete ;
};

#include "BaseFDPFSolver.tpp"

#endif // BASEFDPFSOLVER_H
