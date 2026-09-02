// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include <cassert>

namespace ls2g {

// ---- Phase 1: topology init --------------------------------------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::init_topology(
    const Eigen::Ref<const IntVect> &              slack_ids,
    const Eigen::Ref<const RealVect> &             slack_weights,
    const Eigen::Ref<const IntVect> &              pv,
    const Eigen::Ref<const IntVect> &              pq)
{
    // size the dS value arrays: one coefficient per Ybus nonzero, in Ybus'
    // own storage order (they are pure value arrays, see their declaration --
    // no structure to copy, so nothing of Ybus is duplicated here)
    const Eigen::Index nnz_ybus = Ybus_ref_.nonZeros();
    if(dS_dVm_vals_.size() != nnz_ybus) dS_dVm_vals_.resize(nnz_ybus);
    if(dS_dVa_vals_.size() != nnz_ybus) dS_dVa_vals_.resize(nnz_ybus);
    map_dsdva_r_.clear();
    map_dsdva_i_.clear();
    map_dsdvm_r_.clear();
    map_dsdvm_i_.clear();

    // components compute their own index sets
    base_.init_topology(slack_ids, slack_weights, pv, pq);
    _init_topology_extensions(slack_ids, slack_weights, pv, pq, std::make_index_sequence<sizeof...(Rest)>{});

    // then claim their equations / unknowns in the ledger
    // (allocation order defines the J layout)
    const int n_bus = static_cast<int>(Ybus_ref_.rows());
    ledger_.reset(n_bus);
    base_.register_in(ledger_);
    _register_in_extensions(ledger_, std::make_index_sequence<sizeof...(Rest)>{});
}

// ---- Phase 1.5: per-compute_pf state update ----------------------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::update_state(
    const LSGrid                     * lsgrid_ptr,
    const EigenRefConstCplxSpMat     & Ybus,
    const Eigen::Ref<const CplxVect> & V_init,
    const Eigen::Ref<const CplxVect> & Sbus,
    const Eigen::Ref<const RealVect> & slack_weights)
{
    lsgrid_ptr_ = lsgrid_ptr;
    Ybus_ref_ = Ybus;
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

    // build the matrix (duplicated positions are summed, all values are 0 here).
    // The feature entries are appended to cntrb so that one range describes the
    // whole pattern and Contrib can be fed to setFromTriplets as is; only the
    // first n_ds entries carry a dS value.
    const std::size_t n_ds = cntrb.size();
    cntrb.reserve(n_ds + static_cast<std::size_t>(sink_.size()));
    for (int h = 0; h < sink_.size(); ++h)
        cntrb.push_back(Contrib::structural(sink_.row(h), sink_.col(h)));

    const int dim_J = ledger_.size();
    J_.resize(dim_J, dim_J);
    J_.setFromTriplets(cntrb.begin(), cntrb.end());
    J_.makeCompressed();

    // resolve the dS value maps (Ybus nnz position -> J_.valuePtr() position).
    // J_'s index arrays are read once here, not once per coefficient: this loop
    // runs four times the Ybus nonzeros.
    const int* J_outer = J_.outerIndexPtr();
    const int* J_inner = J_.innerIndexPtr();
    map_dsdva_r_.assign(nnz_Y, -1);
    map_dsdva_i_.assign(nnz_Y, -1);
    map_dsdvm_r_.assign(nnz_Y, -1);
    map_dsdvm_i_.assign(nnz_Y, -1);
    for (std::size_t idx = 0; idx < n_ds; ++idx) {
        const Contrib& c = cntrb[idx];
        const int pos = Base::find_J_pos(J_outer, J_inner, c.jrow(), c.jcol());
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
        feature_pos_[h] = Base::find_J_pos(J_outer, J_inner, sink_.row(h), sink_.col(h));

#ifndef NDEBUG
    // fill_J assigns the dS values and only ACCUMULATES the feature ones, which
    // is correct as long as no Jacobian coefficient is claimed by two dS
    // entries (see _assign_ds) and every stored coefficient has a writer at
    // all. Both are properties of the ledger, so they are checked here, where
    // the layout is decided, rather than assumed in the hot loop.
    {
        std::vector<int> writers(static_cast<std::size_t>(J_.nonZeros()), 0);
        for (std::size_t idx = 0; idx < n_ds; ++idx) {
            const Contrib& c = cntrb[idx];
            const int pos = Base::find_J_pos(J_outer, J_inner, c.jrow(), c.jcol());
            if (pos >= 0) ++writers[static_cast<std::size_t>(pos)];
        }
        for (int w : writers) assert(w <= 1 && "two dS entries share one J coefficient");
        for (int pos : feature_pos_)
            if (pos >= 0) ++writers[static_cast<std::size_t>(pos)];
        for (int w : writers) assert(w >= 1 && "a J coefficient no one writes would keep a stale value");
    }
#endif
}

// ---- Phase 3: fill J values (fast, called every factorisation) ---------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::fill_J()
{
    auto timer = CustTimer();

    real_type* J_values = J_.valuePtr();

    // Only the FEATURE coefficients are zeroed, not the whole value array.
    // This used to be a std::fill over every nonzero of J -- a megabyte on a
    // 9241-bus grid, rewritten at every factorisation -- because every write
    // below accumulated. It no longer has to: the dS passes assign (no two of
    // them ever claim the same coefficient, see _assign_ds), so they bring
    // their own coefficients up to date whatever was there before. What still
    // accumulates is the feature entries, because a component may legitimately
    // add to a dS coefficient -- an hvdc droop slope adds to the dP/dtheta of
    // its end buses -- and those need a zero to start from.
    //
    // Hence the order: zero the feature positions, let the dS pass overwrite
    // the ones it shares, then add the feature values on top. Values are
    // unchanged: `0 + x` is `x`.
    for (int pos : feature_pos_)
        if (pos >= 0) J_values[pos] = static_cast<real_type>(0.);

    const cplx_type* ds_dvm = dS_dVm_vals_.data();
    const cplx_type* ds_dva = dS_dVa_vals_.data();
    _assign_ds<true>(J_values, map_dsdva_r_, ds_dva);
    _assign_ds<false>(J_values, map_dsdva_i_, ds_dva);
    _assign_ds<true>(J_values, map_dsdvm_r_, ds_dvm);
    _assign_ds<false>(J_values, map_dsdvm_i_, ds_dvm);

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
    const Eigen::Index size_dS = V_.size();

    // ---- per-bus quantities ------------------------------------------------
    // Two vectors, both persistent scratch (plain assignment resizes only when
    // the dimension actually changed): Ybus . V, and the reciprocal of |V|,
    // which the pass below multiplies by once per nonzero. The unit phasors
    // V / |V| and the conj(Ibus) . * products the diagonal wants are NOT
    // materialised -- each is used once per bus, on the diagonal, so building
    // them there costs the same arithmetic without more nb_bus vectors sitting
    // in cache for the whole solve.
    Ibus_cache_.noalias() = Ybus_ref_ * V_;      // noalias: no product temporary
    inv_vm_cache_ = V_.array().abs().inverse();  // 1 / |V|

    cplx_type * ds_dvm_val_ptr = dS_dVm_vals_.data();
    cplx_type * ds_dva_val_ptr = dS_dVa_vals_.data();
    // read through raw pointers: these are members, and a compiler that cannot
    // prove the writes below do not reach them would reload each base address
    // on every nonzero
    const cplx_type * const v_ptr     = V_.data();
    const cplx_type * const ibus_ptr  = Ibus_cache_.data();
    const real_type * const invvm_ptr = inv_vm_cache_.data();

    // ---- the dS pass -------------------------------------------------------
    //
    // Both derivatives of the power injected at bus i with respect to the state
    // of bus j come out of ONE product. Writing Y for Ybus(i, j), the formulas
    // derived from pandapower are
    //
    //     dS_dVm(i, j) = conj(Y . Vnorm_j) . V_i    [+ conj(Ibus_i) . Vnorm_i on the diagonal]
    //     dS_dVa(i, j) = -conj(Y . V_j) . i . V_i   [Ibus_i subtracted inside the conjugate there]
    //
    // and Vnorm_j is V_j / |V_j|, so conj(Y . Vnorm_j) . V_i is nothing but
    // conj(Y . V_j) . V_i divided by the real |V_j|. Both entries are therefore
    //
    //     B = conj(Y . V_j) . V_i     (two complex products)
    //     dS_dVm = B / |V_j|          (a real scaling)
    //     dS_dVa = -i . B             (a swap and a sign flip)
    //
    // -- four complex products per nonzero become two. The diagonal's two extra
    // terms share a product as well: subtracting Ibus_i inside the conjugate is
    // conj(Y . V_i - Ibus_i) . V_i == B - conj(Ibus_i) . V_i, and its dS_dVm
    // term conj(Ibus_i) . Vnorm_i is that same conj(Ibus_i) . V_i over |V_i|.
    // So one complex product per bus covers both.
    //
    // The arithmetic is spelled out on real and imaginary parts rather than left
    // to std::complex, whose operator* carries a NaN-recovery branch after every
    // product. Measured together on a 9241-bus grid, this loop goes from 364 to
    // 264 us. Values move by less than one ulp relative (1.9e-16 on that grid):
    // the rounding that differs is dividing by |V_j| after the product rather
    // than before it.
    Eigen::Index pos = 0;
    for (Eigen::Index col_id = 0; col_id < size_dS; ++col_id) {
        // loop-invariant over the inner (row) loop: read once per column
        const real_type V_col_r = std::real(v_ptr[col_id]);
        const real_type V_col_i = std::imag(v_ptr[col_id]);
        const real_type inv_vm_col = invvm_ptr[col_id];
        for (EigenRefConstCplxSpMat::InnerIterator it(Ybus_ref_, col_id); it; ++it) {
            const Eigen::Index row_id = it.row();
            const real_type y_r = std::real(it.value()), y_i = std::imag(it.value());
            const real_type V_row_r = std::real(v_ptr[row_id]), V_row_i = std::imag(v_ptr[row_id]);

            // A = Y . V_col
            const real_type a_r = y_r * V_col_r - y_i * V_col_i;
            const real_type a_i = y_r * V_col_i + y_i * V_col_r;
            // B = conj(A) . V_row
            const real_type b_r = a_r * V_row_r + a_i * V_row_i;
            const real_type b_i = a_r * V_row_i - a_i * V_row_r;

            real_type dvm_r = b_r * inv_vm_col, dvm_i = b_i * inv_vm_col;
            real_type k_r = b_r, k_i = b_i;
            if (col_id == row_id) {
                // conj(Ibus_i) . V_i, then that same value over |V_i|
                const real_type ib_r =  std::real(ibus_ptr[row_id]);
                const real_type ib_i = -std::imag(ibus_ptr[row_id]);   // conj
                const real_type c_r = ib_r * V_row_r - ib_i * V_row_i;
                const real_type c_i = ib_r * V_row_i + ib_i * V_row_r;
                dvm_r += c_r * inv_vm_col;   // col_id == row_id, so this is 1 / |V_i|
                dvm_i += c_i * inv_vm_col;
                k_r -= c_r;
                k_i -= c_i;
            }
            ds_dvm_val_ptr[pos] = cplx_type(dvm_r, dvm_i);
            ds_dva_val_ptr[pos] = cplx_type(k_i, -k_r);   // -i . K
            ++pos;
        }
    }
    timer_dSbus_ += timer.duration();
}

// ---- NR primitives -----------------------------------------------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::mismatch_into(Eigen::Ref<RealVect> res) const
{
    // current state: no step (dx == 0), residual evaluated at V_
    const auto n = static_cast<Eigen::Index>(total_state_variables());
    assert(res.size() == n);
    if(dx_zero_cache_.size() != n) dx_zero_cache_ = RealVect::Zero(n);
    _residual_into(res, V_, dx_zero_cache_);
}

template <typename... Rest>
inline RealVect NRSystem<Base, Rest...>::mismatch() const
{
    RealVect res(static_cast<Eigen::Index>(total_state_variables()));
    mismatch_into(res);
    return res;
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

    _reconstruct_V_into(V_, Va_, Vm_);
    if (Vm_.minCoeff() < static_cast<real_type>(0.)) {
        Vm_ = V_.array().abs();
        Va_ = V_.array().arg();
    }
}

template <typename... Rest>
inline real_type NRSystem<Base, Rest...>::mismatch_sq_norm_at(const Eigen::Ref<const RealVect>& dx) const
{
    const auto n = static_cast<Eigen::Index>(total_state_variables());
    if(res_cache_.size() != n) res_cache_.resize(n);
    _compute_trial_V_into(V_trial_cache_, dx);
    _residual_into(res_cache_, V_trial_cache_, dx);
    return res_cache_.squaredNorm();
}

template <typename... Rest>
inline real_type NRSystem<Base, Rest...>::mismatch_sq_norm_at_current() const
{
    // deliberately the zero-step call, not a shortcut through V_: the trial
    // voltages are rebuilt from (Va_, Vm_) exactly as for any other step, so
    // this is bit-for-bit what mismatch_sq_norm_at(zero vector) returned.
    const auto n = static_cast<Eigen::Index>(total_state_variables());
    if(dx_zero_cache_.size() != n) dx_zero_cache_ = RealVect::Zero(n);
    return mismatch_sq_norm_at(dx_zero_cache_);
}

template <typename... Rest>
inline void NRSystem<Base, Rest...>::_reconstruct_V_into(
    CplxVect& V_out,
    const Eigen::Ref<const RealVect>& Va,
    const Eigen::Ref<const RealVect>& Vm)
{
    // V = Vm * (cos(Va) + i.sin(Va)), straight into the caller's vector: the
    // same expression as before, but assigned instead of returned, so the
    // per-call nb_bus complex temporary is gone (Eigen resizes V_out only when
    // the dimension actually changed, so nothing is allocated after the first
    // call of a topology).
    const cplx_type m_i = BaseConstants::my_i;
    V_out = Vm.array() * (Va.array().cos().template cast<cplx_type>()
                          + m_i * Va.array().sin().template cast<cplx_type>());
}

template <typename... Rest>
inline void NRSystem<Base, Rest...>::_compute_trial_V_into(CplxVect& V_out, const Eigen::Ref<const RealVect>& dx) const
{
    // same voltage loops as apply_step, on (persistent) copies; the components'
    // extra state is handled through the dx argument of adjust_mismatch /
    // fill_custom_rows
    Va_trial_cache_ = Va_;
    Vm_trial_cache_ = Vm_;
    const std::vector<int>& theta_buses = ledger_.theta_buses();
    const std::vector<int>& theta_cols  = ledger_.theta_cols();
    for (size_t k = 0; k < theta_buses.size(); ++k) Va_trial_cache_(theta_buses[k]) += dx(theta_cols[k]);
    const std::vector<int>& vm_buses = ledger_.vm_buses();
    const std::vector<int>& vm_cols  = ledger_.vm_cols();
    for (size_t k = 0; k < vm_buses.size(); ++k) Vm_trial_cache_(vm_buses[k]) += dx(vm_cols[k]);
    _reconstruct_V_into(V_out, Va_trial_cache_, Vm_trial_cache_);
}

template <typename... Rest>
inline void NRSystem<Base, Rest...>::_residual_into(
    Eigen::Ref<RealVect> res,
    const Eigen::Ref<const CplxVect>& V_t,
    const Eigen::Ref<const RealVect>& dx) const
{
    // per-bus complex power mismatch: V .* conj(Ybus V) - Sbus.
    // Ybus * V_t goes into its own persistent buffer first: left inside the
    // coefficient-wise expression it is a sparse-times-dense product, which
    // Eigen can only evaluate into a heap temporary -- one allocation per call,
    // and this runs at least twice per NR iteration (once per backtracking
    // trial of a line search).
    assert(res.size() == static_cast<Eigen::Index>(total_state_variables()));
    ybus_v_cache_.noalias() = Ybus_ref_ * V_t;
    mis_cache_ = V_t.array() * ybus_v_cache_.array().conjugate()
                 - _Sbus_view().array();
    // components adjust the complex injection (e.g. + slack_absorbed * slack_weights,
    // + the theta-dependent hvdc droop flows)
    base_.adjust_mismatch(V_t, dx, mis_cache_);
    _adjust_mismatch_extensions(V_t, dx, mis_cache_, std::make_index_sequence<sizeof...(Rest)>{});

    // generic residual rows, driven by the ledger's (bus, row) pair lists
    // (accumulate: duplicated rows must sum)
    res.setZero();
    const std::vector<int>& p_buses = ledger_.p_buses();
    const std::vector<int>& p_rows  = ledger_.p_rows();
    for (size_t k = 0; k < p_buses.size(); ++k) res(p_rows[k]) -= std::real(mis_cache_(p_buses[k]));
    const std::vector<int>& q_buses = ledger_.q_buses();
    const std::vector<int>& q_rows  = ledger_.q_rows();
    for (size_t k = 0; k < q_buses.size(); ++k) res(q_rows[k]) -= std::imag(mis_cache_(q_buses[k]));

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
}

// ---- value-returning wrappers (kept for out-of-tree derived algorithms) ------
template <typename... Rest>
inline CplxVect NRSystem<Base, Rest...>::_reconstruct_V(const Eigen::Ref<const RealVect>& Va, const Eigen::Ref<const RealVect>& Vm)
{
    CplxVect res;
    _reconstruct_V_into(res, Va, Vm);
    return res;
}

template <typename... Rest>
inline CplxVect NRSystem<Base, Rest...>::_compute_trial_V(const Eigen::Ref<const RealVect>& dx) const
{
    CplxVect res;
    _compute_trial_V_into(res, dx);
    return res;
}

template <typename... Rest>
inline RealVect NRSystem<Base, Rest...>::_residual(const Eigen::Ref<const CplxVect>& V_t, const Eigen::Ref<const RealVect>& dx) const
{
    RealVect res(static_cast<Eigen::Index>(total_state_variables()));
    _residual_into(res, V_t, dx);
    return res;
}

} // namespace ls2g
