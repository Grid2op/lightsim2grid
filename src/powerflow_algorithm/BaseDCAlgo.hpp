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

template<class LinearSolver>
class BaseDCAlgo final: public BaseAlgo
{
    public:
        BaseDCAlgo() noexcept :
            BaseAlgo(false),
            _linear_solver(),
            need_factorize_(true),
            need_refactor_(true),
            timer_refactor_(0.),
            timer_ptdf_(0.),
            timer_lodf_(0.),
            sizeYbus_with_slack_(0),
            sizeYbus_without_slack_(0){};

        virtual ~BaseDCAlgo() noexcept = default;

        virtual void reset();
        virtual void reset_timer(){
            BaseAlgo::reset_timer();
            timer_refactor_ = 0.;
            timer_ptdf_ = 0.;
            timer_lodf_ = 0.;
        }

        virtual TimerJacType get_timers_jacobian() const
        {
            TimerJacType res = {
                timer_Fx_,
                timer_solve_,
                timer_refactor_,
                -1.,  // timer_initialize_: not applicable to DC solver
                timer_check_,
                -1.,  // timer_dSbus_: not applicable to DC solver
                -1.,  // timer_fillJ_: not applicable to DC solver
                -1.,  // timer_Va_Vm_: not applicable to DC solver
                -1.,  // timer_pre_proc_: not applicable to DC solver
                timer_total_nr_
            };
            return res;
        }

        virtual TimerPTDFLODFType get_timers_ptdf_lodf() const
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
                        );

        virtual RealMat get_ptdf();
        virtual RealMat get_lodf(const IntVect & from_bus,
                                 const IntVect & to_bus);
        virtual Eigen::SparseMatrix<real_type> get_bsdf();  // TODO BSDF
        
        virtual void update_internal_Ybus(const Coeff & coeff, bool add){
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

        double timer_refactor_;
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

};

#include "BaseDCAlgo.tpp"

#endif // BASE_DC_ALGO_H
