// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASE_DC_ALGO_H
#define BASE_DC_ALGO_H

#include "BaseAlgo.hpp"

namespace ls2g {

template<class LinearSolver>
class BaseDCAlgo final: public BaseAlgo
{
    public:
        BaseDCAlgo() noexcept :
            BaseAlgo(false),
            _linear_solver(),
            need_factorize_(true),
            need_refactor_(true),
            timer_factor_(0.),
            timer_refactor_(0.),
            timer_initialize_(0.),
            timer_pre_proc_(0.),
            timer_mismatch_(0.),  // used for all the post processing
            timer_ptdf_(0.),
            timer_lodf_(0.),
            sizeYbus_with_slack_(0),
            sizeYbus_without_slack_(0){};

        virtual ~BaseDCAlgo() noexcept = default;

        virtual void reset() final;
        virtual void reset_timer() final{
            BaseAlgo::reset_timer();
            timer_refactor_ = 0.;
            timer_factor_ = 0.;
            timer_initialize_  = 0.;
            timer_pre_proc_  = 0.;
            timer_mismatch_  = 0.;
            timer_solve_ = 0.;
            timer_ptdf_ = 0.;
            timer_lodf_ = 0.;
        }

        virtual TimerJac get_timers_jacobian() const final
        {
            TimerJac res;
            res.timer_Fx_         = timer_Fx_;
            res.timer_solve_      = timer_solve_;
            res.timer_factor_     = timer_factor_;
            res.timer_refactor_   = timer_refactor_;
            res.timer_initialize_ = timer_initialize_;
            res.timer_check_      = timer_check_;
            res.timer_total_nr_   = timer_total_nr_;
            res.timer_pre_proc_   = timer_pre_proc_;
            res.timer_mismatch_   = timer_mismatch_;
            return res;
        }

        virtual TimerPTDFLODFType get_timers_ptdf_lodf() const final
        {
            TimerPTDFLODFType res = {
                timer_ptdf_,  
                timer_lodf_ - timer_ptdf_,
                -1.,  // not available yet so I put -1
            };
            return res;
        }

        // TODO SLACK : this should be handled in Sbus by the gridmodel maybe ?
        virtual
        bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                        CplxVect & V,
                        const CplxVect & Sbus,
                        Eigen::Ref<const IntVect> slack_ids,
                        const RealVect & slack_weights,  // currently unused
                        Eigen::Ref<const IntVect> pv,
                        Eigen::Ref<const IntVect> pq,
                        int max_iter,
                        real_type tol
                        ) final;

        virtual RealMat get_ptdf() final;
        virtual RealMat get_lodf(const IntVect & from_bus,
                                 const IntVect & to_bus) final;
        virtual Eigen::SparseMatrix<real_type> get_bsdf() final;  // TODO BSDF
        
        virtual void update_internal_Ybus(const Coeff & coeff, bool add) final{
            int row_res = static_cast<int>(coeff.row_id);
            row_res = mat_bus_id_(row_res);
            if(row_res == -1) return;
            int col_res = static_cast<int>(coeff.col_id);
            col_res = mat_bus_id_(col_res);
            if(col_res == -1) return;
            real_type val = add ? std::real(coeff.value) : - std::real(coeff.value);
            dcYbus_noslack_.coeffRef(row_res, col_res) += val;

            // need to refactor the linear solver (Ybus changed)
            if(!add) need_refactor_ = true;
        }

    private:
        // no copy allowed
        BaseDCAlgo(const BaseDCAlgo&) = delete;
        BaseDCAlgo(BaseDCAlgo&&) = delete;
        BaseDCAlgo & operator=(BaseDCAlgo&&) = delete;
        BaseDCAlgo & operator=(const BaseDCAlgo&) = delete;

    protected:
        void fill_mat_bus_id(int nb_bus_solver);
        void fill_dcYbus_noslack(int nb_bus_solver, const Eigen::SparseMatrix<cplx_type> & ref_mat);

        // remove_slack_buses: res_mat is initialized and make_compressed in this function
        template<typename ref_mat_type>  // ref_mat_type should be `real_type` or `cplx_type`
        void remove_slack_buses(int nb_bus_solver, const Eigen::SparseMatrix<ref_mat_type> & ref_mat, Eigen::SparseMatrix<real_type> & res_mat);

    protected:
        LinearSolver  _linear_solver;
        bool need_factorize_;
        bool need_refactor_;

        double timer_factor_;
        double timer_refactor_;
        double timer_initialize_;
        double timer_pre_proc_;
        double timer_mismatch_;  // used for all the post processing

        double timer_ptdf_;
        double timer_lodf_;

        // save this not to recompute them when not needed
        int sizeYbus_with_slack_;
        int sizeYbus_without_slack_;
        RealVect dcSbus_noslack_;
        Eigen::SparseMatrix<real_type> dcYbus_noslack_;
        Eigen::VectorXi my_pv_;
        Eigen::VectorXi slack_buses_ids_solver_;
        // -1 if bus is slack , else the id of the row / column used in the linear solver representing this bus
        Eigen::VectorXi mat_bus_id_;   // formerly `ybus_to_me`
        // ybus_ids of non-slack buses in solver order (mat_bus_id_ inverse for non-slack); precomputed in fill_mat_bus_id
        Eigen::VectorXi nonslack_ybus_ids_;

};

#include "BaseDCAlgo.tpp"


} // namespace ls2g

#endif // BASE_DC_ALGO_H