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
    const Eigen::Ref<const IntVect> &              slack_ids,
    const Eigen::Ref<const RealVect> &             slack_weights,
    const Eigen::Ref<const IntVect> &              pv,
    const Eigen::Ref<const IntVect> &              pq)
{
    // init the sparsity pattern of the dS matrices (values do not matter yet)
    // TODO speed: copy only the sparsity pattern and not the values
    dS_dVm_ = Ybus_ref_;  // Ybus_ref_ is a reference, dS_dVm_ is a plain struct, this forces a copy
    dS_dVa_ = Ybus_ref_;
    map_dsdva_r_.clear();
    map_dsdva_i_.clear();
    map_dsdvm_r_.clear();
    map_dsdvm_i_.clear();

    // components compute their own index sets
    base_.init_topology(slack_ids, slack_weights, pv, pq);
    _init_topology_extensions(slack_ids, slack_weights, pv, pq, std::make_index_sequence<sizeof...(Rest)>{});

    // On the python-facing entry point only (see BaseAlgo::deep_validation_),
    // check every bus id the components are about to hand the ledger, because
    // register_in is the first thing to use them as raw indices into its
    // n_bus-sized maps (`theta_col_of_bus_[bus] = col` is an out-of-bounds WRITE
    // for a bus id nobody validated). This is the earliest point at which the
    // full picture exists: update_state has set the extension data, and the
    // component init_topology calls just above the pv / pq / slack index sets.
    //
    // It is deliberately NOT run on the internal path (LSGrid::ac_pf,
    // BaseBatchSolverSynch), where it would be per-element work inside a tight
    // loop and would duplicate checks already made: pv / pq / slack_ids are
    // validated by BaseAlgo::check_pf_inputs, and the extension bus ids by
    // LSGrid::fill_hvdc_droop_solver_data / fill_voltage_control_solver_data,
    // which reject a disconnected bus before ever emitting a solver id. What is
    // NOT covered by either -- the caches going stale against refreshed data --
    // is checked on every solve by validate_registration(), in O(1).
    const int n_bus = static_cast<int>(Ybus_ref_.rows());
    if (deep_validation_) _validate_data(n_bus);

    // then claim their equations / unknowns in the ledger
    // (allocation order defines the J layout)
    ledger_.reset(n_bus);
    base_.register_in(ledger_);
    _register_in_extensions(ledger_, std::make_index_sequence<sizeof...(Rest)>{});
}

// O(1): four size comparisons, no allocation, no per-bus loop
template <typename... Rest>
inline void NRSystem<Base, Rest...>::_validate_shared_sizes(int n_bus) const
{
    if (Va_.size() != n_bus || Vm_.size() != n_bus || V_.size() != n_bus)
        nr_state_error("NRSystem", "the voltage vectors do not have one entry per Ybus bus");
    if (Sbus_size_ != n_bus)
        nr_state_error("NRSystem", "Sbus does not have one entry per Ybus bus");
}

// ---- pre-register_in validation (topology rebuilds only) ---------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::_validate_data(int n_bus) const
{
    _validate_shared_sizes(n_bus);
    base_.validate_data(n_bus);
    _validate_data_extensions(n_bus, std::make_index_sequence<sizeof...(Rest)>{});
}

// ---- Phase 1.5: per-compute_pf state update ----------------------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::update_state(
    const LSGrid                     * lsgrid_ptr,
    const EigenRefConstCplxSpMat     & Ybus,
    const Eigen::Ref<const CplxVect> & V_init,
    const Eigen::Ref<const CplxVect> & Sbus,
    const Eigen::Ref<const RealVect> & slack_weights,
    const AlgoControl                & solver_control)
{
    lsgrid_ptr_ = lsgrid_ptr;
    Ybus_ref_ = Ybus;
    Sbus_data_ptr_ = Sbus.data();
    Sbus_size_ = Sbus.size();

    Va_ = V_init.array().arg();
    Vm_ = V_init.array().abs();
    V_  = V_init;

    // now inform the components
    base_.update_state(lsgrid_ptr, Ybus, Sbus, slack_weights, solver_control);
    _update_state_extensions(lsgrid_ptr, Ybus, Sbus, slack_weights, solver_control, std::make_index_sequence<sizeof...(Rest)>{});
    // NB: no validation here -- MultiSlack / Base only learn their index sets in
    // init_topology, which runs AFTER this. The data pulled above is checked
    // there (pre-register_in) on a rebuild, and by validate_registration()
    // otherwise; both run before anything indexes it.
}

// ---- Phase 1.75: per-solve validation of the registered layout ----------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::validate_registration() const
{
    // This is the ONLY validation on the internal hot path, so it is kept to
    // comparisons of quantities the components already hold: no per-element
    // sweep, no allocation, nothing proportional to the number of buses. The
    // per-element range checks live in _validate_data(), which runs only on the
    // python-facing entry point (see init_topology).
    //
    // What is checked here is the one thing nothing else can catch: update_state
    // re-pulls the extension data from the grid on EVERY solve, while
    // register_in re-derives the row / column caches that index it only on a
    // topology rebuild. If a grid change altered the number of hvdc droop lines
    // or voltage controllers without signalling that rebuild, the two disagree
    // and the per-iteration loops run off the end of the cached arrays --
    // through FeatureWriter, an out-of-bounds write into J. Comparing the sizes
    // catches exactly that, in a fixed handful of integer comparisons.
    _validate_shared_sizes(static_cast<int>(Ybus_ref_.rows()));

    // `ledger_.size()` collapses the row and column counts into one number and
    // only compares them through an assert, which -DNDEBUG removes. Do it here
    // for real: a J built non-square would be resized to (n_rows, n_rows) and
    // then written at a column index beyond it.
    if (ledger_.n_rows() != ledger_.n_cols()) {
        std::ostringstream oss;
        oss << "the components registered " << ledger_.n_rows() << " equations for "
            << ledger_.n_cols() << " unknowns; the Jacobian must be square";
        nr_state_error("NRSystem", oss.str());
    }
    const int dim_J = ledger_.n_rows();
    base_.validate_registration(dim_J);
    _validate_registration_extensions(dim_J, std::make_index_sequence<sizeof...(Rest)>{});
}

// ---- Phase 2: build J sparsity + value maps -----------------------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::build_J_sparsity()
{
    J_ = Eigen::SparseMatrix<real_type, Eigen::ColMajor>();

    // const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>& Ybus = *Ybus_ptr_;
    const int nnz_Y = static_cast<int>(Ybus_ref_.nonZeros());

    // generic dS pass: ONE loop over the Ybus nonzeros generates every
    // dS-derived entry, for all components at once. k is the running nnz index.
    std::vector<Contrib> cntrb;
    cntrb.reserve(4 * nnz_Y);  // pessimistic upper bound
    int k = 0;
    for (int outer = 0; outer < Ybus_ref_.outerSize(); ++outer) {
        for (EigenRefConstCplxSpMat::InnerIterator
            it(Ybus_ref_, outer); it; ++it, ++k)
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
inline void NRSystem<Base, Rest...>::update_trailing_feature_values(int count, const Eigen::Ref<const RealVect>& deltas)
{
    assert(deltas.size() == count);
    const int first_h = static_cast<int>(sink_.size()) - count;
    real_type* J_values = J_.valuePtr();
    for (int i = 0; i < count; ++i)
        J_values[feature_pos_[first_h + i]] += deltas(i);
}

template <typename... Rest>
inline void NRSystem<Base, Rest...>::fill_internal_variables()
{
    auto timer = CustTimer();
    // const Eigen::SparseMatrix<cplx_type>& Ybus = *Ybus_ptr_;
    const auto size_dS = V_.size();
    const CplxVect Vnorm = V_.array() / V_.array().abs();
    const CplxVect Ibus  = Ybus_ref_ * V_;
    const CplxVect conjIbus_Vnorm = Ibus.array().conjugate() * Vnorm.array();

    cplx_type * ds_dvm_val_ptr = dS_dVm_.valuePtr();
    cplx_type * ds_dva_val_ptr = dS_dVa_.valuePtr();

    size_t pos = 0;
    for (size_t col_id = 0; col_id < static_cast<size_t>(size_dS); ++col_id) {
        for (EigenRefConstCplxSpMat::InnerIterator it(Ybus_ref_, col_id); it; ++it) {
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
inline void NRSystem<Base, Rest...>::apply_step(const Eigen::Ref<const RealVect>& dx)
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
inline real_type NRSystem<Base, Rest...>::mismatch_sq_norm_at(const Eigen::Ref<const RealVect>& dx) const
{
    return _residual(_compute_trial_V(dx), dx).squaredNorm();
}

template <typename... Rest>
inline CplxVect NRSystem<Base, Rest...>::_reconstruct_V(const Eigen::Ref<const RealVect>& Va, const Eigen::Ref<const RealVect>& Vm)
{
    const cplx_type m_i = BaseConstants::my_i;
    return Vm.array() * (Va.array().cos().template cast<cplx_type>()
                         + m_i * Va.array().sin().template cast<cplx_type>());
}

template <typename... Rest>
inline CplxVect NRSystem<Base, Rest...>::_compute_trial_V(const Eigen::Ref<const RealVect>& dx) const
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
inline RealVect NRSystem<Base, Rest...>::_residual(const Eigen::Ref<const CplxVect>& V_t, const Eigen::Ref<const RealVect>& dx) const
{
    // per-bus complex power mismatch: V .* conj(Ybus V) - Sbus
    CplxVect mis = V_t.array() * (Ybus_ref_ * V_t).array().conjugate()
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
