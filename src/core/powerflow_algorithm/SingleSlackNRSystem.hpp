// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SINGLE_SLACK_NR_SYSTEM_H
#define SINGLE_SLACK_NR_SYSTEM_H

#include "Utils.hpp"
#include "BaseConstants.hpp"
#include "NRLayout.hpp"
#include "CustTimer.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include <vector>

namespace ls2g {

/**
 * NRSystem implementation for single-slack Newton-Raphson.
 *
 * Jacobian layout: J is (nb_pvpq + nb_pq) x (nb_pvpq + nb_pq).
 * No distributed-slack row/col.  lag == 0.
 *
 * Sign convention: mismatch() returns Sbus - V*conj(Ybus*V)
 * (NEGATED from the classic power-flow residual), so apply_step() uses +=.
 */
class SingleSlackNRSystem {
public:
    static constexpr int lag = 0;

    SingleSlackNRSystem() noexcept
        : Ybus_ptr_(nullptr), Sbus_ptr_(nullptr),
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

    // ----- voltage state accessors ---------------------------------------------

    const CplxVect& V()  const { return V_; }
    const RealVect& Va() const { return Va_; }
    const RealVect& Vm() const { return Vm_; }

    // ----- Jacobian ------------------------------------------------------------

    const Eigen::SparseMatrix<real_type>& J() const { return J_; }
    bool pattern_changed() const { return need_full_rebuild_; }
    void clear_pattern_changed() { need_full_rebuild_ = false; }

    void assemble_jacobian();

    // ----- NR iteration primitives ---------------------------------------------

    RealVect mismatch() const;
    void apply_step(const RealVect& dx);
    real_type mismatch_sq_norm_at(const RealVect& dx) const;

    const NRLayout& layout() const { return layout_; }

    // ----- timers --------------------------------------------------------------

    double timer_dSbus() const { return timer_dSbus_; }
    double timer_fillJ()  const { return timer_fillJ_; }
    void reset_timers() { timer_dSbus_ = 0.; timer_fillJ_ = 0.; }

private:
    const Eigen::SparseMatrix<cplx_type>* Ybus_ptr_;
    const CplxVect*                        Sbus_ptr_;
    Eigen::VectorXi                        pv_;
    Eigen::VectorXi                        pq_;
    Eigen::VectorXi                        pvpq_;
    std::vector<int>                       pvpq_inv_;
    std::vector<int>                       pq_inv_;
    NRLayout                               layout_;

    RealVect Va_;
    RealVect Vm_;
    CplxVect V_;

    Eigen::SparseMatrix<real_type>  J_;
    Eigen::SparseMatrix<cplx_type>  dS_dVm_;
    Eigen::SparseMatrix<cplx_type>  dS_dVa_;
    std::vector<cplx_type*>         value_map_;
    bool                            need_full_rebuild_;

    double timer_dSbus_;
    double timer_fillJ_;

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

    RealVect _mismatch_at(const CplxVect& V_trial) const;
    static CplxVect _reconstruct_V(const RealVect& Va, const RealVect& Vm);
};

#include "SingleSlackNRSystem.tpp"

} // namespace ls2g

#endif // SINGLE_SLACK_NR_SYSTEM_H
