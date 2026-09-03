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
#include "HvdcDroopData.hpp"
#include "linear_solvers/LinearSolverStats.hpp"

namespace ls2g {

template<class LinearSolver>
class BaseDCAlgo final: public BaseAlgo
{
    public:
        // above this (absolute) voltage angle [rad] the DC solution is considered non physical
        // (used as a divergence criterion together with the "all finite" check)
        static constexpr real_type _max_dc_voltage_angle = 1e6;

    public:
        BaseDCAlgo() noexcept :
            BaseAlgo(false),
            _linear_solver(),
            need_factorize_(true),
            need_refactor_(true),
            timer_pre_proc_(0.),
            timer_mismatch_(0.),  // used for all the post processing
            timer_ptdf_(0.),
            timer_lodf_(0.),
            sizeYbus_with_slack_(0),
            sizeYbus_without_slack_(0),
            _lazy_v_(false){};

        ~BaseDCAlgo() noexcept override = default;

        static constexpr bool IS_DC = true;
        bool is_dc() const noexcept override { return IS_DC; }

        void reset() override;
        void reset_timer() override{
            BaseAlgo::reset_timer();
            detail::reset_stats_timers_impl(_linear_solver, 0);
            timer_pre_proc_  = 0.;
            timer_mismatch_  = 0.;
            timer_ptdf_ = 0.;
            timer_lodf_ = 0.;
        }

        TimerJac get_timers_jacobian() const override
        {
            const LinearSolverStats lsstats = detail::get_stats_impl(_linear_solver, 0);
            TimerJac res;
            res.timer_Fx_         = timer_Fx_;
            res.timer_solve_      = lsstats.timer_solve_;
            res.timer_factor_     = lsstats.timer_factor_;
            res.timer_refactor_   = lsstats.timer_refactor_;
            res.timer_initialize_ = lsstats.timer_initialize_;
            res.timer_check_      = timer_check_;
            res.timer_total_nr_   = timer_total_nr_;
            res.timer_pre_proc_   = timer_pre_proc_;
            res.timer_mismatch_   = timer_mismatch_;
            return res;
        }

        // Per-call counters and timings for the underlying linear solver -- see
        // NRAlgo::get_linear_solver_stats for the full description.
        LinearSolverStats get_linear_solver_stats() const override {
            return detail::get_stats_impl(_linear_solver, 0);
        }

        TimerPTDFLODFType get_timers_ptdf_lodf() const override
        {
            TimerPTDFLODFType res = {
                timer_ptdf_,  
                timer_lodf_ - timer_ptdf_,
                -1.,  // not available yet so I put -1
            };
            return res;
        }

        // Native real-valued DC power flow: solves `Bbus . theta = Pbus`.
        // (the DC solver does not implement the complex `compute_pf`: it inherits the throwing
        //  default from BaseAlgo, since every DC code path goes through `compute_pf_dc`)
        bool compute_pf_dc(
            const EigenRefConstRealSpMat     & Bbus,
            const Eigen::Ref<const CplxVect> & V,
            const Eigen::Ref<const RealVect> & Pbus,
            const Eigen::Ref<const IntVect>  & slack_ids,
            const Eigen::Ref<const RealVect> & slack_weights,
            const Eigen::Ref<const IntVect>  & pv,
            const Eigen::Ref<const IntVect>  & pq
        ) override;

        // TOOD speed optim: return refs instead of plain structure
        RealMat get_ptdf() override;
        RealMat get_lodf(
            const Eigen::Ref<const IntVect> & from_bus,
            const Eigen::Ref<const IntVect> & to_bus
        ) override;

        Eigen::SparseMatrix<real_type> get_bsdf() override;  // TODO BSDF

        void update_internal_Ybus(const Coeff & coeff, bool add) override{
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

        // ----- bus masking ("handle disconnected grid" mode) -----------------------
        // ContingencyAnalysis can ask the DC solver to "mask" the buses of a
        // disconnected island: their reduced-system row becomes identity (theta = 0)
        // and their injection is dropped, so the largest connected component still
        // solves while the masked buses are reported as 0. The masking is applied to
        // a working copy of dcYbus_noslack_, so the (incrementally maintained)
        // persistent matrix and the symbolic factorization are left untouched (only a
        // numeric refactorize is needed). See compute_pf_dc.
        bool supports_bus_masking() const override { return true; }
        void set_masked_buses(const std::vector<int> & solver_bus_ids) override{
            masked_buses_ = solver_bus_ids;
        }

        // see BaseAlgo::set_lazy_v / lazy_v
        void set_lazy_v(bool value) override { _lazy_v_ = value; }
        bool lazy_v() const override { return _lazy_v_; }

    private:
        // no copy allowed
        BaseDCAlgo(const BaseDCAlgo&) = delete;
        BaseDCAlgo(BaseDCAlgo&&) = delete;
        BaseDCAlgo & operator=(BaseDCAlgo&&) = delete;
        BaseDCAlgo & operator=(const BaseDCAlgo&) = delete;

    protected:
        void fill_mat_bus_id(int nb_bus_solver);
        void fill_dcYbus_noslack(int nb_bus_solver, const Eigen::Ref<const Eigen::SparseMatrix<real_type>> & ref_mat);

        // hvdc angle-droop ("AC emulation") support: in dc, the linear-mode
        // droop lines contribute `p = p0 + k * (theta1 - theta2)`: the k term
        // goes into the dc matrix (like a branch of susceptance k), the p0
        // constant into the rhs. Saturated lines arrive through Sbus as fixed
        // injections (stamped by the HvdcLineContainer).
        // NB the PTDF / LODF then include the (lossless) droop sensitivity.
        bool update_hvdc_droop_data();  // pulls the data from the grid, true if it changed
        void add_droop_to_dcYbus();
        void add_droop_to_dcSbus();

        // remove_slack_buses: res_mat is resized and its compressed arrays are
        // rewritten from scratch in this function, so it needs a real reference,
        // not Eigen::Ref.
        template<typename ref_mat_type>  // ref_mat_type should be `real_type` or `cplx_type`
        void remove_slack_buses(int nb_bus_solver, const Eigen::Ref<const Eigen::SparseMatrix<ref_mat_type>> & ref_mat, Eigen::SparseMatrix<real_type> & res_mat);

    protected:
        LinearSolver  _linear_solver;
        bool need_factorize_;
        bool need_refactor_;

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

        // connected angle-droop hvdc lines (solver bus ids, pu), refreshed at
        // every compute_pf; a change forces a rebuild of dcYbus / dcSbus
        HvdcDroopSolverData hvdc_droop_data_;

        // solver bus ids (with-slack indexing) masked by the "handle disconnected
        // grid" mode (empty by default => no masking). See set_masked_buses.
        std::vector<int> masked_buses_;

        // batch DC theta-only fast path: see BaseAlgo::set_lazy_v. Deliberately not
        // touched by reset() -- it is a caller-controlled mode, set once per batch
        // compute() and expected to survive the internal reset()s a topology change
        // (eg a per-row contingency) may trigger mid-sweep.
        bool _lazy_v_;
};

#include "BaseDCAlgo.tpp"


} // namespace ls2g

#endif // BASE_DC_ALGO_H