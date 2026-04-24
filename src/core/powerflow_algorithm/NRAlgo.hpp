// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef NR_ALGO_H
#define NR_ALGO_H

#include "BaseAlgo.hpp"
#include "NRLayout.hpp"
#include "SlackPolicies.hpp"
#include "ScalingPolicies.hpp"

namespace ls2g {

/**
 * Unified Newton-Raphson solver parameterised by LinearSolver and SlackPolicy.
 *
 * SlackPolicy is either MultiSlackPolicy (lag=1, distributed slack) or
 * SingleSlackPolicy (lag=0, traditional single slack).  The policy is granted
 * friend access so its static template methods can reach protected members directly.
 */
template<class LinearSolver, class SlackPolicy, class ScalingPolicy = NoScalingPolicy>
class NRAlgo : public BaseAlgo
{
    friend SlackPolicy;
    friend ScalingPolicy;

public:
    NRAlgo() noexcept :
        BaseAlgo(true),
        need_factorize_(true),
        timer_refactor_(0.),
        timer_initialize_(0.),
        timer_dSbus_(0.),
        timer_fillJ_(0.),
        timer_Va_Vm_(0.),
        timer_pre_proc_(0.) {}

    virtual ~NRAlgo() noexcept = default;

    virtual
    Eigen::Ref<const Eigen::SparseMatrix<real_type>> get_J() const override {
        return J_;
    }

    virtual
    Eigen::SparseMatrix<real_type> get_J_python() const {
        Eigen::SparseMatrix<real_type> res = get_J();
        return res;
    }

    virtual
    TimerJacType get_timers_jacobian() const override
    {
        auto res = TimerJacType(timer_Fx_,
                                timer_solve_,
                                timer_refactor_,
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
    bool compute_pf(const Eigen::SparseMatrix<cplx_type>& Ybus,
                    CplxVect& V,
                    const CplxVect& Sbus,
                    Eigen::Ref<const IntVect> slack_ids,
                    const RealVect& slack_weights,
                    Eigen::Ref<const IntVect> pv,
                    Eigen::Ref<const IntVect> pq,
                    int max_iter,
                    real_type tol
                    ) override;

    virtual void reset() override;

    ScalingPolicy& get_scaling_policy() { return _scaling_policy_; }
    const ScalingPolicy& get_scaling_policy() const { return _scaling_policy_; }

    Eigen::SparseMatrix<real_type>
    create_jacobian_matrix_test(const Eigen::SparseMatrix<cplx_type>& Ybus,
                                const CplxVect& V,
                                const RealVect& slack_weights,
                                const Eigen::VectorXi& pq,
                                const Eigen::VectorXi& pvpq)
    {
        // DO NOT USE, FOR DEBUG ONLY (especially for multiple slacks)
        const auto& n_pvpq = pvpq.size();
        const auto& n_pq = pq.size();
        std::vector<int> pvpq_inv(V.size(), -1);
        for(int inv_id=0; inv_id < (int)n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
        std::vector<int> pq_inv(V.size(), -1);
        for(int inv_id=0; inv_id < (int)n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;
        const auto slack_ids_vec = extract_slack_bus_id(pvpq, pq, static_cast<unsigned int>(V.size()));
        const size_t slack_bus_id = static_cast<size_t>(slack_ids_vec(0));
        SlackPolicy::fill_jacobian_matrix(*this, Ybus, V, slack_bus_id, slack_weights, pq, pvpq, pq_inv, pvpq_inv);
        return J_;
    }

protected:
    virtual void reset_timer() override {
        BaseAlgo::reset_timer();
        timer_refactor_ = 0.;
        timer_dSbus_ = 0.;
        timer_fillJ_ = 0.;
        timer_Va_Vm_ = 0.;
        timer_pre_proc_ = 0.;
        timer_initialize_ = 0.;
    }

    void _dSbus_dV(const Eigen::Ref<const Eigen::SparseMatrix<cplx_type>>& Ybus,
                   const Eigen::Ref<const CplxVect>& V);

    void _get_values_J(int& nb_obj_this_col,
                       std::vector<Eigen::Index>& inner_index,
                       std::vector<real_type>& values,
                       const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
                       const std::vector<int>& index_row_inv,
                       const Eigen::VectorXi& index_col,
                       size_t col_id,
                       size_t row_lag,
                       size_t col_lag);

    void _get_values_J(int& nb_obj_this_col,
                       std::vector<Eigen::Index>& inner_index,
                       std::vector<real_type>& values,
                       const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
                       const std::vector<int>& index_row_inv,
                       size_t col_id_mat,
                       size_t row_lag,
                       size_t col_lag);

    void _do_voltage_update(const RealVect& dx,
                            const Eigen::VectorXi& my_pv,
                            Eigen::Ref<const IntVect> pq);

    void reset_if_needed() {
        if(err_ != ErrorType::NoError ||
           _solver_control.need_reset_solver() ||
           _solver_control.has_dimension_changed() ||
           _solver_control.ybus_change_sparsity_pattern() ||
           _solver_control.has_ybus_some_coeffs_zero() ||
           _solver_control.has_slack_participate_changed() ||
           _solver_control.has_pv_changed() ||
           _solver_control.has_pq_changed())
        {
            reset();
        }
    }

protected:
    // used linear solver
    LinearSolver _linear_solver;

    // Jacobian and power-flow derivatives
    Eigen::SparseMatrix<real_type> J_;
    Eigen::SparseMatrix<cplx_type> dS_dVm_;
    Eigen::SparseMatrix<cplx_type> dS_dVa_;
    bool need_factorize_;

    // pointer map from J_ entries back into dS_dVm_ / dS_dVa_
    std::vector<cplx_type*> value_map_;

    // block-index arithmetic for the NR state vector and Jacobian
    NRLayout _layout;

    // scaling policy instance (configurable at runtime for MaxVoltageChangePolicy)
    ScalingPolicy _scaling_policy_;

    // extra timers (beyond BaseAlgo's four)
    double timer_refactor_;
    double timer_initialize_;
    double timer_dSbus_;
    double timer_fillJ_;
    double timer_Va_Vm_;
    double timer_pre_proc_;

private:
    NRAlgo(const NRAlgo&) = delete;
    NRAlgo(NRAlgo&&) = delete;
    NRAlgo& operator=(NRAlgo&&) = delete;
    NRAlgo& operator=(const NRAlgo&) = delete;

    void print_J(int min_col=-1, int max_col=-1) const {
        auto size_J = J_.cols();
        if(max_col == -1) max_col = static_cast<int>(size_J);
        if(min_col == -1) min_col = 0;
        for(int col_id=min_col; col_id < max_col; ++col_id){
            for(Eigen::SparseMatrix<real_type>::InnerIterator it(J_, col_id); it; ++it)
                std::cout << it.row() << ", " << it.col() << ": " << it.value() << std::endl;
            std::cout << std::endl;
        }
    }

    Eigen::SparseMatrix<cplx_type>
    _make_diagonal_matrix(const Eigen::Ref<const CplxVect>& diag_val) {
        auto n = diag_val.size();
        Eigen::SparseMatrix<cplx_type> res(n, n);
        res.setIdentity();
        res.diagonal() = diag_val;
        return res;
    }
};

#include "NRAlgo.tpp"

} // namespace ls2g

#endif // NR_ALGO_H
