// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef MULTI_SLACK_NR_SYSTEM_H
#define MULTI_SLACK_NR_SYSTEM_H

#include "Utils.hpp"
#include "BaseConstants.hpp"
#include "NRLayout.hpp"
#include "CustTimer.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include <vector>

namespace ls2g {

/**
 * NRSystem implementation for distributed (multi) slack Newton-Raphson.
 *
 * Jacobian layout: J is (1 + nb_pvpq + nb_pq) x (1 + nb_pvpq + nb_pq).
 * Row/col 0 is the distributed-slack equation; the rest are the standard
 * dP/dQ blocks.  lag == 1.
 *
 * Sign convention: mismatch() returns Sbus - V*conj(Ybus*V) - correction
 * (NEGATED from the classic power-flow residual), so apply_step() uses +=.
 */
class MultiSlackNRSystem {
public:
    static constexpr int lag = 1;

    MultiSlackNRSystem() noexcept
        : Ybus_ptr_(nullptr), Sbus_ptr_(nullptr),
          slack_bus_id_(0), slack_absorbed_(0.),
          need_full_rebuild_(true),
          timer_dSbus_(0.), timer_fillJ_(0.) {}

    // ----- setup / state -------------------------------------------------------

    void setup(const Eigen::SparseMatrix<cplx_type>& Ybus,
               const CplxVect& V_init,
               const CplxVect& Sbus,
               Eigen::Ref<const IntVect> slack_ids,
               const RealVect& slack_weights,
               Eigen::Ref<const IntVect> pv,
               Eigen::Ref<const IntVect> pq);

    void clear_jacobian() {
        J_         = Eigen::SparseMatrix<real_type>();
        dS_dVm_    = Eigen::SparseMatrix<cplx_type>();
        dS_dVa_    = Eigen::SparseMatrix<cplx_type>();
        value_map_.clear();
        need_full_rebuild_ = true;
    }

    // ----- voltage state accessors (current iterate) ---------------------------

    const CplxVect& V()  const { return V_; }
    const RealVect& Va() const { return Va_; }
    const RealVect& Vm() const { return Vm_; }

    // ----- Jacobian ------------------------------------------------------------

    const Eigen::SparseMatrix<real_type>& J() const { return J_; }
    bool pattern_changed() const { return need_full_rebuild_; }
    void clear_pattern_changed() { need_full_rebuild_ = false; }

    void assemble_jacobian();

    // ----- NR iteration primitives ---------------------------------------------

    /** Negated residual: Sbus - V*conj(Ybus*V) - slack_correction */
    RealVect mismatch() const;

    /** Apply Newton step (+=  convention): Va += theta(dx), Vm += vm(dx) */
    void apply_step(const RealVect& dx);

    /**
     * Non-destructive probe: evaluate ||mismatch(x + dx)||^2 without
     * committing the update.  Used by step-scaling policies.
     */
    real_type mismatch_sq_norm_at(const RealVect& dx) const;

    /** NRLayout accessor (needed by NRAlgo for step-norm computation). */
    const NRLayout& layout() const { return layout_; }

    // ----- timers (accumulated across iterations) ------------------------------

    double timer_dSbus() const { return timer_dSbus_; }
    double timer_fillJ()  const { return timer_fillJ_; }
    void reset_timers() { timer_dSbus_ = 0.; timer_fillJ_ = 0.; }

private:
    // ---- stored problem parameters (set once per compute_pf) ------------------
    const Eigen::SparseMatrix<cplx_type>* Ybus_ptr_;
    const CplxVect*                        Sbus_ptr_;
    size_t                                 slack_bus_id_;
    RealVect                               slack_weights_;
    Eigen::VectorXi                        my_pv_;
    Eigen::VectorXi                        pq_;
    Eigen::VectorXi                        pvpq_;
    std::vector<int>                       pvpq_inv_;
    std::vector<int>                       pq_inv_;
    NRLayout                               layout_;

    // ---- iteration state ------------------------------------------------------
    RealVect  Va_;
    RealVect  Vm_;
    CplxVect  V_;
    real_type slack_absorbed_;

    // ---- Jacobian internals ---------------------------------------------------
    Eigen::SparseMatrix<real_type>  J_;
    Eigen::SparseMatrix<cplx_type>  dS_dVm_;
    Eigen::SparseMatrix<cplx_type>  dS_dVa_;
    std::vector<cplx_type*>         value_map_;
    bool                            need_full_rebuild_;

    // ---- timers ---------------------------------------------------------------
    double timer_dSbus_;
    double timer_fillJ_;

    // ---- private helpers ------------------------------------------------------
    void _dSbus_dV(const Eigen::SparseMatrix<cplx_type>& Ybus, const CplxVect& V);

    void _fill_jac_unknown_sparsity();
    void _fill_jac_known_sparsity();
    void _fill_value_map(bool reset_J);

    void _get_values_J(int& nb_obj_this_col,
                       std::vector<Eigen::Index>& inner_index,
                       std::vector<real_type>& values,
                       const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
                       const std::vector<int>& index_row_inv,
                       const Eigen::VectorXi& index_col,
                       size_t col_id,
                       size_t row_lag,
                       size_t col_lag) const;

    void _get_values_J(int& nb_obj_this_col,
                       std::vector<Eigen::Index>& inner_index,
                       std::vector<real_type>& values,
                       const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
                       const std::vector<int>& index_row_inv,
                       size_t col_id_mat,
                       size_t row_lag,
                       size_t col_lag) const;

    // Negated residual at an arbitrary voltage vector (used by mismatch_sq_norm_at)
    RealVect _mismatch_at(const CplxVect& V_trial, real_type slack_absorbed_trial) const;

    static CplxVect _reconstruct_V(const RealVect& Va, const RealVect& Vm);
};

#include "MultiSlackNRSystem.tpp"

} // namespace ls2g

#endif // MULTI_SLACK_NR_SYSTEM_H
