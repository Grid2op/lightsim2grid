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
    const RealVect&                        slack_weights,
    Eigen::Ref<const IntVect>              pv,
    Eigen::Ref<const IntVect>              pq)
{

    // init the sparsity pattern
    // I think we don't really care about the
    // values
    dS_dVm_ = *Ybus_ptr_;
    dS_dVa_ = *Ybus_ptr_;
    map_dsdva_r_.clear();
    map_dsdva_i_.clear();
    map_dsdvm_r_.clear();
    map_dsdvm_i_.clear();
    
    // now init the extra features
    base_.init_topology(slack_ids, slack_weights, pv, pq);
    _init_topology_extensions(slack_ids, slack_weights, pv, pq, std::make_index_sequence<sizeof...(Rest)>{});
    
    // now compute the jacobian size
    _update_total_state_variables(std::make_index_sequence<sizeof...(Rest)>{});

    // build the bus_id -> Jacobian column converters (base block + extensions)
    const int n_bus = static_cast<int>(Ybus_ptr_->rows());
    theta_to_J_col_.assign(n_bus, -1);
    vm_to_J_col_.assign(n_bus, -1);
    q_to_J_col_.assign(n_bus, -1);
    base_.fill_col_converters(theta_to_J_col_, vm_to_J_col_, q_to_J_col_);
    _fill_col_converters_extensions(theta_to_J_col_, vm_to_J_col_, q_to_J_col_,
                                    std::make_index_sequence<sizeof...(Rest)>{});

    need_full_rebuild_ = true;
}

// ---- Phase 1.5: per-compute_pf state update ----------------------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::update_state(
    const LSGrid *                         lsgrid_ptr,
    const Eigen::SparseMatrix<cplx_type>&  Ybus,
    const CplxVect&                        V_init,
    const CplxVect&                        Sbus,
    Eigen::Ref<const RealVect>             slack_weights)
{
    lsgrid_ptr_ = lsgrid_ptr;
    Ybus_ptr_ = &Ybus;
    Sbus_ptr_ = &Sbus;

    Va_ = V_init.array().arg();
    Vm_ = V_init.array().abs();
    V_  = V_init;

    // now inform the extensions
    base_.update_state(lsgrid_ptr, Ybus, Sbus, slack_weights);
    _update_state_extensions(lsgrid_ptr, Ybus, Sbus, slack_weights, std::make_index_sequence<sizeof...(Rest)>{});
}

// ---- Phase 2: build J sparsity + value_map -----------------------------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::build_J_sparsity()
{
    // reset J
    J_         = Eigen::SparseMatrix<real_type>();

    // compute its sparsity pattern
    const Eigen::SparseMatrix<cplx_type> & Ybus = *Ybus_ptr_;
    const int n_bus  = Ybus.rows();
    const size_t dim_J  =  total_state_variables(); // nb_pvpq_ + nb_pq_;
    // base_.build_J_sparsity();

    // const Eigen::SparseMatrix<cplx_type> & Ybus = *Ybus_ptr_;
    // const int nnz_Y  = Ybus.nonZeros();

    // // get the triplets (will be in a virtual function later)
    // std::vector< std::vector<Contrib> > cijs(4); // stores in order c11, c12, c21 and c22
    // int c11 = 0, c12 = 1, c21 = 2, c22 = 3;  // TODO NR refacto: move this as static attrs

    // // TODO this might be overly pessimistic (only pvpq comp are kept for c11 and c21)
    // // TODO this might be overly pessimistic (only pq comp are kept for c12 and c22) 
    // for(auto & el: cijs) el.reserve(nnz_Y);
    // int k = 0;
    // for (int outer = 0; outer < Ybus.outerSize(); ++outer) {
    //     for (Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>::InnerIterator
    //          it(Ybus, outer); it; ++it, ++k)
    //     {
    //         int i = (int)it.row(), j = (int)it.col();
    //         int ri = pvpq_inv_[i], rq = pq_inv_[i];
    //         int ci = pvpq_inv_[j], cq = pq_inv_[j];
    //         if (ri >= 0 && ci >= 0) cijs[c11].push_back({ri,          ci,          k});
    //         if (ri >= 0 && cq >= 0) cijs[c12].push_back({ri,          nb_pvpq_ + cq, k});
    //         if (rq >= 0 && ci >= 0) cijs[c21].push_back({nb_pvpq_ + rq, ci,          k});
    //         if (rq >= 0 && cq >= 0) cijs[c22].push_back({nb_pvpq_ + rq, nb_pvpq_ + cq, k});
    //     }
    // }
    std::vector< std::vector<Contrib> > contribs;
    contribs.reserve(sizeof...(Rest) + 1);
    contribs.push_back(base_.build_J_contrib());
    _build_J_contrib_extensions(contribs, std::make_index_sequence<sizeof...(Rest)>{});

    size_t expected_size = 0;
    for(const auto & el: contribs) expected_size += el.size();
    // std::cout << "\texpected_size " << expected_size << std::endl;
    // reserve enough space
    std::vector<Eigen::Triplet<real_type> > triplets;
    triplets.reserve(expected_size);

    // now fill the triplets
    for(const auto& cij : contribs)
    {
        for (auto& c : cij)
        {
            triplets.push_back({c.jrow(), c.jcol(), 0.});
            // if((c.jrow() >= 23) || (c.jcol() >= 23)){
            //     std::cout << "\t error in NRSystem.tpp: " << c.jrow() << " " << c.jcol() << std::endl;
            // }
        }
    }

    // and build the matrix
    // std::cout << "and build the matrix, dim_J: " << dim_J << std::endl;
    J_.resize(dim_J, dim_J);
    J_.setFromTriplets(triplets.begin(), triplets.end());
    J_.makeCompressed();

    // and finally build the value maps
    // std::cout << "_build_value_map " << std::endl;
    _build_value_map(contribs);  // will call build_value_map_extensions
}

// ---- Phase 3: fill J values (fast, called every factorisation) ---------------
template <typename... Rest>
inline void NRSystem<Base, Rest...>::fill_J()
{
    auto timer = CustTimer();

    const cplx_type* ds_dvm = dS_dVm_.valuePtr();
    const cplx_type* ds_dva = dS_dVa_.valuePtr();
    size_t i = 0;
    for(auto & c : map_dsdva_r_){
        if(c == -1){
            // coeff of J11 not used in J
            i++;
            continue;
        }
        J_.valuePtr()[c] = std::real(ds_dva[i]);
        i++;
    }

    i = 0;
    for(auto & c : map_dsdva_i_){
        if(c == -1){
            // coeff of J21 not used in J
            i++;
            continue;
        }
        J_.valuePtr()[c] = std::imag(ds_dva[i]);
        i++;
    }

    i = 0;
    for(auto & c : map_dsdvm_r_){
        if(c == -1){
            // coeff of J12 not used in J
            i++;
            continue;
        }
        J_.valuePtr()[c] = std::real(ds_dvm[i]);
        i++;
    }

    i = 0;
    for(auto & c : map_dsdvm_i_){
        if(c == -1){
            // coeff of J22 not used in J
            i++;
            continue;
        }
        J_.valuePtr()[c] = std::imag(ds_dvm[i]);
        i++;
    }
    _fill_J_extensions(std::make_index_sequence<sizeof...(Rest)>{});
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
    for (size_t col_id = 0; col_id < size_dS; ++col_id) {
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
    return _residual(V_, RealVect::Zero(total_state_variables()));
}

template <typename... Rest>
inline void NRSystem<Base, Rest...>::apply_step(const RealVect& dx)
{
    // base block: theta at pv/pq, vm at pq
    if (base_.nb_pv() > 0)
        Va_(base_.pv()) += theta_base(dx).segment(0, base_.nb_pv());
    if (base_.nb_pq() > 0) {
        Va_(base_.pq()) += theta_base(dx).segment(base_.nb_pv(), base_.nb_pq());
        Vm_(base_.pq()) += vm_base(dx);
    }
    // extension blocks: e.g. pv_slack angles + slack absorbed
    _apply_step_extensions(dx, Va_, Vm_, std::make_index_sequence<sizeof...(Rest)>{});

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
    // TODO check that dx has the proper size
    RealVect Va_t = Va_;
    RealVect Vm_t = Vm_;
    if (base_.nb_pv() > 0) Va_t(base_.pv()) += theta_base(dx).segment(0, base_.nb_pv());
    if (base_.nb_pq() > 0) {
        Va_t(base_.pq()) += theta_base(dx).segment(base_.nb_pv(), base_.nb_pq());
        Vm_t(base_.pq()) += vm_base(dx);
    }
    // extension blocks may also carry angle state (e.g. pv_slack)
    _apply_trial_voltages_extensions(dx, Va_t, Vm_t, std::make_index_sequence<sizeof...(Rest)>{});
    return _reconstruct_V(Va_t, Vm_t);
}

template <typename... Rest>
inline RealVect NRSystem<Base, Rest...>::_residual(const CplxVect& V_t, const RealVect& dx) const
{
    // per-bus complex power mismatch: V .* conj(Ybus V) - Sbus
    CplxVect mis = V_t.array() * (*Ybus_ptr_ * V_t).array().conjugate()
                   - Sbus_ptr_->array();
    // extensions adjust the complex injection (e.g. + slack_absorbed * slack_weights)
    _adjust_mismatch_extensions(dx, mis, std::make_index_sequence<sizeof...(Rest)>{});

    const RealVect real_ = mis.real();
    const RealVect imag_ = mis.imag();

    RealVect res(total_state_variables());
    base_.fill_mismatch(res.segment(0, base_.get_size()), real_, imag_);
    _fill_mismatch_extensions(res, real_, imag_, std::make_index_sequence<sizeof...(Rest)>{});
    return res;
}

} // namespace ls2g
