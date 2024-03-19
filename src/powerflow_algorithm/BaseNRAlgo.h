// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASE_NR_ALGO_H
#define BASE_NR_ALGO_H

#include "BaseAlgo.h"

/**
Base class for Newton Raphson based solver
**/
template<class LinearSolver>
class BaseNRAlgo : public BaseAlgo
{
    public:
        BaseNRAlgo():
            BaseAlgo(true),
            need_factorize_(true),
            timer_initialize_(0.),
            timer_dSbus_(0.),
            timer_fillJ_(0.),
            timer_Va_Vm_(0.),
            timer_pre_proc_(0.){}

        virtual
        Eigen::Ref<const Eigen::SparseMatrix<real_type> > get_J() const {
            return J_;
        }
        
        virtual
        Eigen::SparseMatrix<real_type> get_J_python() const {
            Eigen::SparseMatrix<real_type> res = get_J();
            return res;
        }

        virtual
        TimerJacType get_timers_jacobian() const
        {
            // TODO refacto that, and change the order
            auto res = TimerJacType(timer_Fx_,
                                    timer_solve_,
                                    timer_initialize_,
                                    timer_check_,
                                    timer_dSbus_,
                                    timer_fillJ_,
                                    timer_Va_Vm_,
                                    timer_pre_proc_,
                                    timer_total_nr_);
            return res;
        }

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
                        ) ;

        virtual void reset();

    protected:
        virtual void reset_timer(){
            BaseAlgo::reset_timer();
            timer_dSbus_ = 0.;
            timer_fillJ_ = 0.;
            timer_Va_Vm_ = 0.;
            timer_pre_proc_ = 0.;
            timer_initialize_ = 0.;
        }
        virtual
        void initialize(){
            auto timer = CustTimer();
            n_ = static_cast<int>(J_.cols()); // should be equal to J_.nrows()
            err_ = ErrorType::NoError; // reset error message
            const ErrorType init_status = _linear_solver.initialize(J_);
            if(init_status != ErrorType::NoError){
                // std::cout << "init_ok " << init_ok << std::endl;
                err_ = init_status;
            }
            need_factorize_ = false;
            timer_initialize_ += timer.duration();
        }

        virtual
        void solve(RealVect & b, bool has_just_been_inialized){
            auto timer = CustTimer();
            const ErrorType solve_status = _linear_solver.solve(J_, b, has_just_been_inialized);
            if(solve_status != ErrorType::NoError){
                // std::cout << "solve error: " << solve_status << std::endl;
                err_ = solve_status;
            }
            timer_solve_ += timer.duration();
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
                            const Eigen::VectorXi & pvpq,
                            bool reset_J);

        void reset_if_needed(){
            if(_solver_control.need_reset_solver() || 
               _solver_control.has_dimension_changed() ||
               _solver_control.ybus_change_sparsity_pattern() ||
               _solver_control.has_ybus_some_coeffs_zero() ||
               _solver_control.has_slack_participate_changed() ||
               _solver_control.has_pv_changed() ||
               _solver_control.has_pq_changed()
               ){
               reset();
            }
        }
    protected:
        // used linear solver
        LinearSolver _linear_solver;

        // solution of the problem
        Eigen::SparseMatrix<real_type> J_;  // the jacobian matrix
        Eigen::SparseMatrix<cplx_type> dS_dVm_;
        Eigen::SparseMatrix<cplx_type> dS_dVa_;
        bool need_factorize_;

        // to store the mapping from the element of J_ in dS_dVm_ and dS_dVa_
        // it does not own any memory at all !
        std::vector<cplx_type*> value_map_;
        // std::vector<int> col_map_;
        // std::vector<int> row_map_;

        // timers
        double timer_initialize_;
        double timer_dSbus_;
        double timer_fillJ_;
        double timer_Va_Vm_;
        double timer_pre_proc_;


    Eigen::SparseMatrix<real_type>
        create_jacobian_matrix_test(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                    const CplxVect & V,
                                    const RealVect & slack_weights,
                                    const Eigen::VectorXi & pq,
                                    const Eigen::VectorXi & pvpq
                                    ){
            // DO NOT USE, FOR DEBUG ONLY (especially for multiple slacks)
            const auto & n_pvpq = pvpq.size();
            const auto & n_pq = pvpq.size();
            std::vector<int> pvpq_inv(V.size(), -1);
            for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
            std::vector<int> pq_inv(V.size(), -1);
            for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;
            // TODO if bug when using it, check the "pvpq" below, 
            // in theory its "pv" !
            int slack_bus_id = extract_slack_bus_id(pvpq, pq,
                                                    static_cast<unsigned int>(V.size())
                                                    );
            fill_jacobian_matrix(Ybus, V, static_cast<Eigen::Index>(slack_bus_id),
                                 slack_weights, pq, pvpq, pq_inv, pvpq_inv);
            return J_;
        }

    private:
        // no copy allowed
        BaseNRAlgo( const BaseNRAlgo & ) =delete ;
        BaseNRAlgo & operator=( const BaseNRAlgo & ) =delete ;

        /** helper function to print the max_col left most columns of the J matrix **/
        void print_J(int min_col=-1, int max_col=-1) const{
            auto size_J = J_.cols();
            if(max_col == -1) max_col = static_cast<int>(size_J);
            if(min_col == -1) min_col = 0;
            for (int col_id=min_col; col_id < max_col; ++col_id){
                for (Eigen::SparseMatrix<real_type>::InnerIterator it(J_, col_id); it; ++it)
                {
                    std::cout << it.row() << ", " << it.col() << ": " << it.value() << std::endl;
                }
                std::cout << std::endl;
            }
        }


    private:
        Eigen::SparseMatrix<cplx_type>
            _make_diagonal_matrix(const Eigen::Ref<const CplxVect > & diag_val){
            // TODO their might be a more efficient way to do that
            auto n = diag_val.size();
            Eigen::SparseMatrix<cplx_type> res(n,n);
            // first method, without a loop of mine
            res.setIdentity();  // segfault if attempt to use this function without this
            res.diagonal() = diag_val;

            // second method, with an "optimized" loop
            // res.reserve(Eigen::VectorXi::Constant(n,1)); // i reserve one number per columns (speed optim)
            // for(unsigned int i = 0; i < n; ++i){
            //    res.insert(i,i) = diag_val(i);
            //}
            return res;
        }
};

#include "BaseNRAlgo.tpp"

#endif // BASE_NR_ALGO_H
