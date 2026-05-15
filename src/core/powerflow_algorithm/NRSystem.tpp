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
    return _mismatch_core(V_);
}

template <typename... Rest>
inline void NRSystem<Base, Rest...>::apply_step(const RealVect& dx)
{
    if (base_.nb_pv() > 0)
        Va_(base_.pv()) += theta(dx).segment(0, base_.nb_pv());
    if (base_.nb_pq() > 0) {
        Va_(base_.pq()) += theta(dx).segment(base_.nb_pv(), base_.nb_pq());
        Vm_(base_.pq()) += vm(dx);
    }
    V_ = _reconstruct_V(Va_, Vm_);
    if (Vm_.minCoeff() < static_cast<real_type>(0.)) {
        Vm_ = V_.array().abs();
        Va_ = V_.array().arg();
    }
}

template <typename... Rest>
inline real_type NRSystem<Base, Rest...>::mismatch_sq_norm_at(const RealVect& dx) const
{
    return _mismatch_core(_compute_trial_V(dx)).squaredNorm();
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
    RealVect Va_t = Va_;
    RealVect Vm_t = Vm_;
    if (base_.nb_pv() > 0) Va_t(base_.pv()) += theta(dx).segment(0, base_.nb_pv());
    if (base_.nb_pq() > 0) {
        Va_t(base_.pq()) += theta(dx).segment(base_.nb_pv(), base_.nb_pq());
        Vm_t(base_.pq()) += vm(dx);
    }
    return _reconstruct_V(Va_t, Vm_t);
}

template <typename... Rest>
inline RealVect NRSystem<Base, Rest...>::_mismatch_core(const CplxVect& V_trial) const
{
    auto mis = V_trial.array() * (*Ybus_ptr_ * V_trial).array().conjugate()
               - Sbus_ptr_->array();
    const RealVect real_ = mis.real();
    const RealVect imag_ = mis.imag();

    RealVect res(base_.get_size());
    res.segment(0,               base_.nb_pv()) = -real_(base_.pv());
    res.segment(base_.nb_pv(),          base_.nb_pq()) = -real_(base_.pq());
    res.segment(base_.nb_pv() + base_.nb_pq(), base_.nb_pq()) = -imag_(base_.pq());
    return res;
}

// // =============================================================================
// //  NRSystem<MultiSlack, Rest...> — distributed-slack extension
// // =============================================================================

// // ---- Phase 1: topology init --------------------------------------------------

// template <typename... Rest>
// void NRSystem<MultiSlack, Rest...>::init_topology(
//     const Eigen::SparseMatrix<cplx_type>& Ybus,
//     const CplxVect&                        Sbus,
//     Eigen::Ref<const IntVect>              slack_ids,
//     const RealVect&                        slack_weights,
//     Eigen::Ref<const IntVect>              pv,
//     Eigen::Ref<const IntVect>              pq)
// {
//     // std::cout << "NRSystem<MultiSlack>::init_topology \n";
//     slack_bus_id_  = static_cast<size_t>(slack_ids(0));
//     slack_weights_ = slack_weights;

//     // Build my_pv = (extra slack buses beyond the primary) ++ pv
//     const int nb_slack_added = static_cast<int>(slack_ids.size()) - 1;
//     IntVect my_pv;
//     if (nb_slack_added > 0) {
//         my_pv.resize(static_cast<int>(pv.size()) + nb_slack_added);
//         for (int i = 0; i < nb_slack_added; ++i) my_pv(i) = slack_ids(i + 1);
//         for (int i = 0; i < static_cast<int>(pv.size()); ++i)
//             my_pv(i + nb_slack_added) = pv(i);
//     } else {
//         my_pv = IntVect(pv);
//     }

//     Base::init_topology(Ybus, Sbus, slack_ids, slack_weights, my_pv, pq);
//     // this->lag_ += 1;   // MultiSlack always contributes exactly one extra row/col
// }

// // ---- Phase 1.5: per-compute_pf state update ----------------------------------

// template <typename... Rest>
// void NRSystem<MultiSlack, Rest...>::update_state(
//     const Eigen::SparseMatrix<cplx_type>& Ybus,
//     const CplxVect&                        V_init,
//     const CplxVect&                        Sbus)
// {
//     Base::update_state(Ybus, V_init, Sbus);
//     slack_absorbed_ = std::real(Sbus.sum());
// }

// // ---- NR primitives -----------------------------------------------------------

// template <typename... Rest>
// RealVect NRSystem<MultiSlack, Rest...>::mismatch() const
// {
//     return _mismatch_with_slack(Base::V(), slack_absorbed_);
// }

// template <typename... Rest>
// void NRSystem<MultiSlack, Rest...>::apply_step(const RealVect& dx)
// {
//     Base::apply_step(dx);
//     slack_absorbed_ += dx(_J_slack_row());
// }

// template <typename... Rest>
// real_type NRSystem<MultiSlack, Rest...>::mismatch_sq_norm_at(const RealVect& dx) const
// {
//     const real_type sa_trial = slack_absorbed_ + dx(_J_slack_row());
//     return _mismatch_with_slack(this->_compute_trial_V(dx), sa_trial).squaredNorm();
// }

// // ---- Virtual hooks -----------------------------------------------------------


// // ---- Private helpers ---------------------------------------------------------

// template <typename... Rest>
// void NRSystem<MultiSlack, Rest...>::_append_slack_triplets(
//     std::vector<Eigen::Triplet<double>>& coeffs) const
// {
//     const int slack_row = _J_slack_row();
//     const int slack_col = this->total() - 1;     // last column
//     const int n_pvpq    = this->theta_size();

//     const Eigen::SparseMatrix<real_type> dS_dVa_r = Base::dS_dVa().real();
//     const Eigen::SparseMatrix<real_type> dS_dVm_r = Base::dS_dVm().real();
//     const int nb_bus = dS_dVa_r.cols();

//     // Slack bus row: dP_slack / dTheta  (dTheta columns 0..n_pvpq-1)
//     for (int col_id = 0; col_id < nb_bus; ++col_id) {
//         const int J_col = this->pvpq_inv()[col_id];
//         if (J_col < 0) continue;
//         for (Eigen::SparseMatrix<real_type>::InnerIterator it(dS_dVa_r, col_id); it; ++it) {
//             if (it.row() != static_cast<Eigen::Index>(slack_bus_id_)) continue;
//             coeffs.push_back(Eigen::Triplet<double>(slack_row, J_col, it.value()));
//         }
//     }
//     // Slack bus row: dP_slack / dVm  (dVm columns n_pvpq..n_pvpq+n_pq-1)
//     for (int col_id = 0; col_id < nb_bus; ++col_id) {
//         const int J_col = Base::pq_inv()[col_id];
//         if (J_col < 0) continue;
//         for (Eigen::SparseMatrix<real_type>::InnerIterator it(dS_dVm_r, col_id); it; ++it) {
//             if (it.row() != static_cast<Eigen::Index>(slack_bus_id_)) continue;
//             coeffs.push_back(Eigen::Triplet<double>(slack_row, J_col + n_pvpq, it.value()));
//         }
//     }

//     // Slack column: diagonal entry (slack_row, slack_col)
//     coeffs.push_back(Eigen::Triplet<double>(
//         slack_row, slack_col,
//         slack_weights_[static_cast<int>(slack_bus_id_)]));

//     // Slack column: coupling to pvpq rows (rows 0..n_pvpq-1)
//     for (int i = 0; i < static_cast<int>(Base::pvpq().size()); ++i) {
//         const real_type sl_w = slack_weights_(Base::pvpq()(i));
//         if (std::abs(sl_w) > BaseConstants::_tol_equal_float)
//             coeffs.push_back(Eigen::Triplet<double>(i, slack_col, sl_w));
//     }
// }

// template <typename... Rest>
// RealVect NRSystem<MultiSlack, Rest...>::_mismatch_with_slack(
//     const CplxVect& V_trial, real_type sa) const
// {
//     const int n_pv = Base::nb_pv();
//     const int n_pq = Base::nb_pq();

//     auto mis = V_trial.array() * (*this->Ybus_ptr_ * V_trial).array().conjugate()
//                - this->Sbus_ptr_->array()
//                + sa * slack_weights_.array();
//     const RealVect real_ = mis.real();
//     const RealVect imag_ = mis.imag();

//     // Layout: [dP_pv, dP_pq, dQ_pq, slack_eq]  (negated)
//     RealVect res(n_pv + 2 * n_pq + 1);
//     res.segment(0,               n_pv) = -real_(Base::pv());
//     res.segment(n_pv,            n_pq) = -real_(Base::pq());
//     res.segment(n_pv + n_pq,     n_pq) = -imag_(Base::pq());
//     res(n_pv + 2 * n_pq) = -real_(static_cast<int>(slack_bus_id_));
//     return res;
// }

} // namespace ls2g
