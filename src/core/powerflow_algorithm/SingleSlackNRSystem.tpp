// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// ---- setup ------------------------------------------------------------------

inline void SingleSlackNRSystem::setup(
    const Eigen::SparseMatrix<cplx_type>& Ybus,
    const CplxVect& V_init,
    const CplxVect& Sbus,
    Eigen::Ref<const IntVect> /*slack_ids*/,
    const RealVect& /*slack_weights*/,
    Eigen::Ref<const IntVect> pv,
    Eigen::Ref<const IntVect> pq)
{
    Ybus_ptr_ = &Ybus;
    Sbus_ptr_ = &Sbus;

    pv_ = IntVect(pv);
    pq_ = IntVect(pq);

    const int n_pv   = static_cast<int>(pv_.size());
    const int n_pq   = static_cast<int>(pq_.size());
    pvpq_.resize(n_pv + n_pq);
    pvpq_ << pv_, pq_;

    const int n_pvpq = static_cast<int>(pvpq_.size());
    pvpq_inv_.assign(V_init.size(), -1);
    for (int i = 0; i < n_pvpq; ++i) pvpq_inv_[pvpq_(i)] = i;
    pq_inv_.assign(V_init.size(), -1);
    for (int i = 0; i < n_pq;   ++i) pq_inv_[pq_(i)] = i;

    layout_ = NRLayout(n_pv, n_pq, lag);

    Va_ = V_init.array().arg();
    Vm_ = V_init.array().abs();
    V_  = V_init;
}

// ---- assemble_jacobian ------------------------------------------------------

inline void SingleSlackNRSystem::assemble_jacobian()
{
    _dSbus_dV(*Ybus_ptr_, V_);

    auto timer = CustTimer();
    const int n_pvpq = static_cast<int>(pvpq_.size());
    const int n_pq   = static_cast<int>(pq_.size());
    const int size_j = n_pvpq + n_pq;  // lag=0

    if (J_.cols() != size_j) {
        _fill_jac_unknown_sparsity();
        _fill_value_map(false);
        need_full_rebuild_ = true;
    } else {
        if (value_map_.empty()) _fill_value_map(true);
        _fill_jac_known_sparsity();
        need_full_rebuild_ = false;
    }
    timer_fillJ_ += timer.duration();
}

// ---- mismatch ---------------------------------------------------------------

inline RealVect SingleSlackNRSystem::mismatch() const
{
    return _mismatch_at(V_);
}

// ---- apply_step -------------------------------------------------------------

inline void SingleSlackNRSystem::apply_step(const RealVect& dx)
{
    if (layout_.nb_pv() > 0)
        Va_(pv_) += layout_.theta(dx).segment(0, layout_.nb_pv());
    if (layout_.nb_pq() > 0) {
        Va_(pq_) += layout_.theta(dx).segment(layout_.nb_pv(), layout_.nb_pq());
        Vm_(pq_) += layout_.vm(dx);
    }
    V_ = _reconstruct_V(Va_, Vm_);
    if (Vm_.minCoeff() < static_cast<real_type>(0.)) {
        Vm_ = V_.array().abs();
        Va_ = V_.array().arg();
    }
}

// ---- mismatch_sq_norm_at ----------------------------------------------------

inline real_type SingleSlackNRSystem::mismatch_sq_norm_at(const RealVect& dx) const
{
    RealVect Va_trial = Va_;
    RealVect Vm_trial = Vm_;
    if (layout_.nb_pv() > 0)
        Va_trial(pv_) += layout_.theta(dx).segment(0, layout_.nb_pv());
    if (layout_.nb_pq() > 0) {
        Va_trial(pq_) += layout_.theta(dx).segment(layout_.nb_pv(), layout_.nb_pq());
        Vm_trial(pq_) += layout_.vm(dx);
    }
    CplxVect V_trial = _reconstruct_V(Va_trial, Vm_trial);
    return _mismatch_at(V_trial).squaredNorm();
}

// ---- private helpers --------------------------------------------------------

inline CplxVect SingleSlackNRSystem::_reconstruct_V(const RealVect& Va, const RealVect& Vm)
{
    const cplx_type m_i = BaseConstants::my_i;
    return Vm.array() * (Va.array().cos().template cast<cplx_type>()
                         + m_i * Va.array().sin().template cast<cplx_type>());
}

inline RealVect SingleSlackNRSystem::_mismatch_at(const CplxVect& V_trial) const
{
    const int n_pv = static_cast<int>(pv_.size());
    const int n_pq = static_cast<int>(pq_.size());

    auto mis = V_trial.array() * (*Ybus_ptr_ * V_trial).array().conjugate()
               - Sbus_ptr_->array();
    const RealVect real_ = mis.real();
    const RealVect imag_ = mis.imag();

    // Result layout: [dP_pv, dP_pq, dQ_pq]  (negated)
    RealVect res(n_pv + 2 * n_pq);
    res.segment(0, n_pv)          = -real_(pv_);
    res.segment(n_pv, n_pq)       = -real_(pq_);
    res.segment(n_pv + n_pq, n_pq) = -imag_(pq_);
    return res;
}

inline void SingleSlackNRSystem::_dSbus_dV(
    const Eigen::SparseMatrix<cplx_type>& Ybus, const CplxVect& V)
{
    auto timer = CustTimer();
    const auto size_dS = V.size();
    const CplxVect Vnorm = V.array() / V.array().abs();
    const CplxVect Ibus = Ybus * V;
    const CplxVect conjIbus_Vnorm = Ibus.array().conjugate() * Vnorm.array();

    if (dS_dVm_.cols() != Ybus.cols()) {
        dS_dVm_ = Ybus;
        dS_dVa_ = Ybus;
    }

    cplx_type* ds_dvm_x = dS_dVm_.valuePtr();
    cplx_type* ds_dva_x = dS_dVa_.valuePtr();

    unsigned int pos = 0;
    for (int col_id = 0; col_id < size_dS; ++col_id) {
        for (Eigen::SparseMatrix<cplx_type>::InnerIterator it(Ybus, col_id); it; ++it) {
            const int row_id = static_cast<int>(it.row());
            const cplx_type el = it.value();

            cplx_type& dvm = ds_dvm_x[pos];
            cplx_type& dva = ds_dva_x[pos];

            dvm = el * Vnorm(col_id);
            dvm = std::conj(dvm) * V(row_id);

            dva = el * V(col_id);
            if (col_id == row_id) {
                dvm += conjIbus_Vnorm(row_id);
                ds_dva_x[pos] -= Ibus(row_id);
            }
            const cplx_type tmp = BaseConstants::my_i * V(row_id);
            dva = std::conj(-dva) * tmp;
            ++pos;
        }
    }
    timer_dSbus_ += timer.duration();
}

inline void SingleSlackNRSystem::_fill_jac_unknown_sparsity()
{
    const int n_pvpq = static_cast<int>(pvpq_.size());
    const int n_pq   = static_cast<int>(pq_.size());
    const int size_j = n_pvpq + n_pq;

    const Eigen::SparseMatrix<real_type> dS_dVa_r = dS_dVa_.real();
    const Eigen::SparseMatrix<real_type> dS_dVa_i = dS_dVa_.imag();
    const Eigen::SparseMatrix<real_type> dS_dVm_r = dS_dVm_.real();
    const Eigen::SparseMatrix<real_type> dS_dVm_i = dS_dVm_.imag();

    using StorageIndex = Eigen::SparseMatrix<cplx_type>::StorageIndex;

    J_ = Eigen::SparseMatrix<real_type>(size_j, size_j);

    std::vector<Eigen::Triplet<double>> coeffs;
    coeffs.reserve(2 * (dS_dVa_.nonZeros() + dS_dVm_.nonZeros()));

    int nb_obj = 0;
    std::vector<Eigen::Index> inner_index;
    std::vector<real_type> values;

    for (int col_id = 0; col_id < n_pvpq; ++col_id) {
        nb_obj = 0; inner_index.clear(); values.clear();
        _get_values_J(nb_obj, inner_index, values, dS_dVa_r, pvpq_inv_, pvpq_,
                      col_id, static_cast<size_t>(layout_.J_dP_row()), 0);
        _get_values_J(nb_obj, inner_index, values, dS_dVa_i, pq_inv_, pvpq_,
                      col_id, static_cast<size_t>(layout_.J_dQ_row()), 0);
        for (int k = 0; k < nb_obj; ++k)
            coeffs.push_back(Eigen::Triplet<double>(
                static_cast<int>(inner_index[k]), col_id, values[k]));
    }

    for (int col_id = 0; col_id < n_pq; ++col_id) {
        nb_obj = 0; inner_index.clear(); values.clear();
        _get_values_J(nb_obj, inner_index, values, dS_dVm_r, pvpq_inv_, pq_,
                      col_id, static_cast<size_t>(layout_.J_dP_row()), 0);
        _get_values_J(nb_obj, inner_index, values, dS_dVm_i, pq_inv_, pq_,
                      col_id, static_cast<size_t>(layout_.J_dQ_row()), 0);
        for (int k = 0; k < nb_obj; ++k)
            coeffs.push_back(Eigen::Triplet<double>(
                static_cast<int>(inner_index[k]), col_id + n_pvpq, values[k]));
    }

    J_.setFromTriplets(coeffs.begin(), coeffs.end());
    J_.makeCompressed();
}

inline void SingleSlackNRSystem::_fill_value_map(bool reset_J)
{
    const int n_pvpq = static_cast<int>(pvpq_.size());
    value_map_.clear();
    value_map_.reserve(J_.nonZeros());

    const int n_col = static_cast<int>(J_.cols());
    for (int col_ = 0; col_ < n_col; ++col_) {
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(J_, col_); it; ++it) {
            const int row_id = static_cast<int>(it.row());
            const int col_id = static_cast<int>(it.col());
            if (reset_J) it.valueRef() = static_cast<real_type>(0.);

            if (col_id < n_pvpq && row_id < n_pvpq) {
                value_map_.push_back(&dS_dVa_.coeffRef(pvpq_[row_id], pvpq_[col_id]));
            } else if (col_id < n_pvpq && row_id >= n_pvpq) {
                value_map_.push_back(&dS_dVa_.coeffRef(pq_[row_id - n_pvpq], pvpq_[col_id]));
            } else if (col_id >= n_pvpq && row_id < n_pvpq) {
                value_map_.push_back(&dS_dVm_.coeffRef(pvpq_[row_id], pq_[col_id - n_pvpq]));
            } else {
                value_map_.push_back(&dS_dVm_.coeffRef(pq_[row_id - n_pvpq], pq_[col_id - n_pvpq]));
            }
        }
    }
    dS_dVa_.makeCompressed();
    dS_dVm_.makeCompressed();
}

inline void SingleSlackNRSystem::_fill_jac_known_sparsity()
{
    const int n_pvpq_threshold = layout_.J_dQ_row();  // = theta_size for lag=0
    const int n_cols = static_cast<int>(J_.cols());
    unsigned int pos = 0;
    for (int col_id = 0; col_id < n_cols; ++col_id) {
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(J_, col_id); it; ++it) {
            it.valueRef() = it.row() < n_pvpq_threshold
                            ? std::real(*value_map_[pos])
                            : std::imag(*value_map_[pos]);
            ++pos;
        }
    }
}

inline void SingleSlackNRSystem::_get_values_J(
    int& nb_obj_this_col,
    std::vector<Eigen::Index>& inner_index,
    std::vector<real_type>& values,
    const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
    const std::vector<int>& index_row_inv,
    const Eigen::VectorXi& index_col,
    size_t col_id,
    size_t row_lag,
    size_t col_lag) const
{
    const int col_id_mat = index_col(static_cast<int>(col_id + col_lag));
    _get_values_J(nb_obj_this_col, inner_index, values, mat, index_row_inv,
                  col_id_mat, row_lag, col_lag);
}

inline void SingleSlackNRSystem::_get_values_J(
    int& nb_obj_this_col,
    std::vector<Eigen::Index>& inner_index,
    std::vector<real_type>& values,
    const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
    const std::vector<int>& index_row_inv,
    size_t col_id_mat,
    size_t row_lag,
    size_t /*col_lag*/) const
{
    const int start_id = mat.outerIndexPtr()[col_id_mat];
    const int end_id   = mat.outerIndexPtr()[col_id_mat + 1];
    const real_type* val_ptr = mat.valuePtr();
    for (int obj_id = start_id; obj_id < end_id; ++obj_id) {
        const int row_id_dS = mat.innerIndexPtr()[obj_id];
        const int row_id    = index_row_inv[row_id_dS];
        if (row_id >= 0) {
            inner_index.push_back(static_cast<Eigen::Index>(row_id) + row_lag);
            values.push_back(val_ptr[obj_id]);
            ++nb_obj_this_col;
        }
    }
}
