// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

namespace ls2g {

// ---- Phase 1: topology init --------------------------------------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::init_topology(
    Eigen::Ref<const IntVect>              slack_ids,
    Eigen::Ref<const RealVect>             slack_weights,
    Eigen::Ref<const IntVect>              pv,
    Eigen::Ref<const IntVect>              pq)
{
    // init the sparsity pattern of the dS matrices (values do not matter yet)
    dS_dVm_ = *Ybus_ptr_;
    dS_dVa_ = *Ybus_ptr_;
    map_dsdva_r_.clear();
    map_dsdva_i_.clear();
    map_dsdvm_r_.clear();
    map_dsdvm_i_.clear();

    // components compute their own index sets
    base_.init_topology(slack_ids, slack_weights, pv, pq);
    _init_topology_extensions(slack_ids, slack_weights, pv, pq, std::make_index_sequence<sizeof...(Rest)>{});

    // then claim their equations / unknowns in the ledger
    // (allocation order defines the J layout)
    const int n_bus = static_cast<int>(Ybus_ptr_->rows());
    ledger_.reset(n_bus);
    base_.register_in(ledger_);
    _register_in_extensions(ledger_, std::make_index_sequence<sizeof...(Rest)>{});
}

// ---- Phase 1.5: per-compute_pf state update ----------------------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::update_state(
    const LSGrid *                         lsgrid_ptr,
    const Eigen::SparseMatrix<cplx_type>&  Ybus,
    const CplxVect&                        V_init,
    Eigen::Ref<const CplxVect>              Sbus,
    Eigen::Ref<const RealVect>             slack_weights)
{
    lsgrid_ptr_ = lsgrid_ptr;
    Ybus_ptr_ = &Ybus;
    Sbus_data_ptr_ = Sbus.data();
    Sbus_size_ = Sbus.size();

    Va_ = V_init.array().arg();
    Vm_ = V_init.array().abs();
    V_  = V_init;

    // now inform the components
    base_.update_state(lsgrid_ptr, Ybus, Sbus, slack_weights);
    _update_state_extensions(lsgrid_ptr, Ybus, Sbus, slack_weights, std::make_index_sequence<sizeof...(Rest)>{});
}

// ---- Phase 2: build J sparsity + value maps -----------------------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::build_J_sparsity()
{
    J_ = Eigen::SparseMatrix<real_type, Eigen::ColMajor>();

    const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>& Ybus = *Ybus_ptr_;
    const int nnz_Y = static_cast<int>(Ybus.nonZeros());

    // generic dS pass: ONE loop over the Ybus nonzeros generates every
    // dS-derived entry, for all components at once. k is the running nnz index.
    std::vector<Contrib> cntrb;
    cntrb.reserve(4 * nnz_Y);  // pessimistic upper bound
    int k = 0;
    for (int outer = 0; outer < Ybus.outerSize(); ++outer) {
        for (Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>::InnerIterator
            it(Ybus, outer); it; ++it, ++k)
        {
            const int i = static_cast<int>(it.row()), j = static_cast<int>(it.col());
            const int pr = ledger_.p_row(i),     qr = ledger_.q_row(i);
            const int tc = ledger_.theta_col(j), vc = ledger_.vm_col(j);
            if (pr >= 0 && tc >= 0) cntrb.push_back({pr, tc, k, dSdVa_r});
            if (pr >= 0 && vc >= 0) cntrb.push_back({pr, vc, k, dSdVm_r});
            if (qr >= 0 && tc >= 0) cntrb.push_back({qr, tc, k, dSdVa_i});
            if (qr >= 0 && vc >= 0) cntrb.push_back({qr, vc, k, dSdVm_i});
        }
    }

    // feature (non dS-derived) entries declared by the components
    sink_.clear();
    base_.declare_feature_entries(sink_);
    _declare_feature_entries_extensions(sink_, std::make_index_sequence<sizeof...(Rest)>{});

    // build the matrix (duplicated positions are summed, all values are 0 here)
    std::vector<Eigen::Triplet<real_type> > triplets;
    triplets.reserve(cntrb.size() + sink_.size());
    for (const auto& c : cntrb) triplets.push_back({c.jrow(), c.jcol(), 0.});
    for (int h = 0; h < sink_.size(); ++h) triplets.push_back({sink_.row(h), sink_.col(h), 0.});

    const int dim_J = ledger_.size();
    J_.resize(dim_J, dim_J);
    J_.setFromTriplets(triplets.begin(), triplets.end());
    J_.makeCompressed();

    // resolve the dS value maps (Ybus nnz position -> J_.valuePtr() position)
    map_dsdva_r_.assign(nnz_Y, -1);
    map_dsdva_i_.assign(nnz_Y, -1);
    map_dsdvm_r_.assign(nnz_Y, -1);
    map_dsdvm_i_.assign(nnz_Y, -1);
    for (const auto& c : cntrb) {
        const int pos = Base::find_J_pos(J_, c.jrow(), c.jcol());
        switch (c.whichmat()) {
            case dSdVa_r: map_dsdva_r_[c.ybus_k()] = pos; break;
            case dSdVa_i: map_dsdva_i_[c.ybus_k()] = pos; break;
            case dSdVm_r: map_dsdvm_r_[c.ybus_k()] = pos; break;
            case dSdVm_i: map_dsdvm_i_[c.ybus_k()] = pos; break;
        }
    }

    // resolve the feature positions
    feature_pos_.resize(sink_.size());
    for (int h = 0; h < sink_.size(); ++h)
        feature_pos_[h] = Base::find_J_pos(J_, sink_.row(h), sink_.col(h));
}

// ---- Phase 3: fill J values (fast, called every factorisation) ---------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::fill_J()
{
    auto timer = CustTimer();

    // zero J first: every write below accumulates (+=), since several
    // contributions may resolve to the same position
    real_type* J_values = J_.valuePtr();
    std::fill(J_values, J_values + J_.nonZeros(), static_cast<real_type>(0.));

    const cplx_type* ds_dvm = dS_dVm_.valuePtr();
    const cplx_type* ds_dva = dS_dVa_.valuePtr();
    size_t i = 0;
    for (auto& c : map_dsdva_r_) {
        if (c != -1) J_values[c] += std::real(ds_dva[i]);
        ++i;
    }
    i = 0;
    for (auto& c : map_dsdva_i_) {
        if (c != -1) J_values[c] += std::imag(ds_dva[i]);
        ++i;
    }
    i = 0;
    for (auto& c : map_dsdvm_r_) {
        if (c != -1) J_values[c] += std::real(ds_dvm[i]);
        ++i;
    }
    i = 0;
    for (auto& c : map_dsdvm_i_) {
        if (c != -1) J_values[c] += std::imag(ds_dvm[i]);
        ++i;
    }

    // feature (non dS-derived) values
    FeatureWriter writer(J_values, feature_pos_);
    base_.fill_feature_values(writer, Va_);
    _fill_feature_values_extensions(writer, std::make_index_sequence<sizeof...(Rest)>{});

    // bus masking: overwrite the masked buses' rows with the identity (zero the
    // whole row, then set the diagonal to 1 -- ones last, since the diagonal is
    // itself part of a masked row). Pure value-level edit, J sparsity unchanged.
    if (!masked_buses_.empty()) {
        if (masked_dirty_) _recompute_mask_positions();
        for (int p : masked_zero_pos_) J_values[p] = static_cast<real_type>(0.);
        for (int p : masked_one_pos_)  J_values[p] = static_cast<real_type>(1.);
    }

    timer_fillJ_ += timer.duration();
}

template <typename... Rest>
inline void NRSystem<Base, Rest...>::fill_internal_variables()
{
    auto timer = CustTimer();
    const Eigen::SparseMatrix<cplx_type>& Ybus = *Ybus_ptr_;
    const auto size_dS = V_.size();
    const CplxVect Vnorm = V_.array() / V_.array().abs();
    const CplxVect Ibus  = Ybus * V_;
    const CplxVect conjIbus_Vnorm = Ibus.array().conjugate() * Vnorm.array();

    cplx_type * ds_dvm_val_ptr = dS_dVm_.valuePtr();
    cplx_type * ds_dva_val_ptr = dS_dVa_.valuePtr();

    size_t pos = 0;
    for (size_t col_id = 0; col_id < static_cast<size_t>(size_dS); ++col_id) {
        for (Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>::InnerIterator it(Ybus, col_id); it; ++it) {
            const size_t row_id = static_cast<size_t>(it.row());
            const cplx_type el_ybus = it.value();

            cplx_type& dvm = ds_dvm_val_ptr[pos];
            cplx_type& dva = ds_dva_val_ptr[pos];

            // use formula derived from pandapower
            dvm = el_ybus * Vnorm(col_id);
            dvm = std::conj(dvm) * V_(row_id);

            dva = el_ybus * V_(col_id);
            if (col_id == row_id) {
                dvm += conjIbus_Vnorm(row_id);
                dva -= Ibus(row_id);
            }
            const cplx_type tmp = BaseConstants::my_i * V_(row_id);
            dva = std::conj(-dva) * tmp;
            ++pos;
        }
    }
    timer_dSbus_ += timer.duration();
}

// ---- NR primitives -----------------------------------------------------------
template <typename... Rest>
inline RealVect NRSystem<Base, Rest...>::mismatch() const
{
    // current state: no step (dx == 0), residual evaluated at V_
    const auto n = static_cast<Eigen::Index>(total_state_variables());
    if(dx_zero_cache_.size() != n) dx_zero_cache_ = RealVect::Zero(n);
    return _residual(V_, dx_zero_cache_);
}

template <typename... Rest>
inline void NRSystem<Base, Rest...>::apply_step(const RealVect& dx)
{
    // generic voltage updates, driven by the ledger's (bus, col) pair lists
    const std::vector<int>& theta_buses = ledger_.theta_buses();
    const std::vector<int>& theta_cols  = ledger_.theta_cols();
    for (size_t k = 0; k < theta_buses.size(); ++k) Va_(theta_buses[k]) += dx(theta_cols[k]);
    const std::vector<int>& vm_buses = ledger_.vm_buses();
    const std::vector<int>& vm_cols  = ledger_.vm_cols();
    for (size_t k = 0; k < vm_buses.size(); ++k) Vm_(vm_buses[k]) += dx(vm_cols[k]);

    // component-owned non-voltage state (e.g. slack absorbed)
    base_.apply_step(dx);
    _apply_step_extensions(dx, std::make_index_sequence<sizeof...(Rest)>{});

    V_ = _reconstruct_V(Va_, Vm_);
    if (Vm_.minCoeff() < static_cast<real_type>(0.)) {
        Vm_ = V_.array().abs();
        Va_ = V_.array().arg();
    }
}

template <typename... Rest>
inline real_type NRSystem<Base, Rest...>::mismatch_sq_norm_at(const RealVect& dx) const
{
    return _residual(_compute_trial_V(dx), dx).squaredNorm();
}

template <typename... Rest>
inline CplxVect NRSystem<Base, Rest...>::_reconstruct_V(const RealVect& Va, const RealVect& Vm)
{
    const cplx_type m_i = BaseConstants::my_i;
    return Vm.array() * (Va.array().cos().template cast<cplx_type>()
                         + m_i * Va.array().sin().template cast<cplx_type>());
}

template <typename... Rest>
inline CplxVect NRSystem<Base, Rest...>::_compute_trial_V(const RealVect& dx) const
{
    // same voltage loops as apply_step, on local copies; the components' extra
    // state is handled through the dx argument of adjust_mismatch / fill_custom_rows
    RealVect Va_t = Va_;
    RealVect Vm_t = Vm_;
    const std::vector<int>& theta_buses = ledger_.theta_buses();
    const std::vector<int>& theta_cols  = ledger_.theta_cols();
    for (size_t k = 0; k < theta_buses.size(); ++k) Va_t(theta_buses[k]) += dx(theta_cols[k]);
    const std::vector<int>& vm_buses = ledger_.vm_buses();
    const std::vector<int>& vm_cols  = ledger_.vm_cols();
    for (size_t k = 0; k < vm_buses.size(); ++k) Vm_t(vm_buses[k]) += dx(vm_cols[k]);
    return _reconstruct_V(Va_t, Vm_t);
}

template <typename... Rest>
inline RealVect NRSystem<Base, Rest...>::_residual(const CplxVect& V_t, const RealVect& dx) const
{
    // per-bus complex power mismatch: V .* conj(Ybus V) - Sbus
    CplxVect mis = V_t.array() * (*Ybus_ptr_ * V_t).array().conjugate()
                   - _Sbus_view().array();
    // components adjust the complex injection (e.g. + slack_absorbed * slack_weights,
    // + the theta-dependent hvdc droop flows)
    base_.adjust_mismatch(V_t, dx, mis);
    _adjust_mismatch_extensions(V_t, dx, mis, std::make_index_sequence<sizeof...(Rest)>{});

    // generic residual rows, driven by the ledger's (bus, row) pair lists
    // (accumulate: duplicated rows must sum)
    RealVect res = RealVect::Zero(static_cast<Eigen::Index>(total_state_variables()));
    const std::vector<int>& p_buses = ledger_.p_buses();
    const std::vector<int>& p_rows  = ledger_.p_rows();
    for (size_t k = 0; k < p_buses.size(); ++k) res(p_rows[k]) -= std::real(mis(p_buses[k]));
    const std::vector<int>& q_buses = ledger_.q_buses();
    const std::vector<int>& q_rows  = ledger_.q_rows();
    for (size_t k = 0; k < q_buses.size(); ++k) res(q_rows[k]) -= std::imag(mis(q_buses[k]));

    // component-owned custom rows (none for Base / MultiSlack)
    base_.fill_custom_rows(res, Va_, Vm_, dx);
    _fill_custom_rows_extensions(res, Va_, Vm_, dx, std::make_index_sequence<sizeof...(Rest)>{});

    // bus masking: the masked buses' P/Q residuals are forced to 0 so the trivial
    // identity rows of J (see fill_J) yield dx == 0 on those buses.
    if (!masked_buses_.empty()) {
        for (int b : masked_buses_) {
            const int pr = ledger_.p_row(b);
            if (pr >= 0) res(pr) = static_cast<real_type>(0.);
            const int qr = ledger_.q_row(b);
            if (qr >= 0) res(qr) = static_cast<real_type>(0.);
        }
    }
    return res;
}

} // namespace ls2g
