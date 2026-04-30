// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

namespace ls2g {

// =============================================================================
//  NRSystem<> — base specialisation (single-slack, lag = 0)
// =============================================================================

// ---- Phase 1: topology init --------------------------------------------------

inline void NRSystem<>::init_topology(  // TODO rename
    const Eigen::SparseMatrix<cplx_type>& Ybus,
    const CplxVect&                        Sbus,
    Eigen::Ref<const IntVect>              /*slack_ids*/,
    const RealVect&                        /*slack_weights*/,
    Eigen::Ref<const IntVect>              pv,
    Eigen::Ref<const IntVect>              pq)
{
    Ybus_ptr_ = &Ybus;
    Sbus_ptr_ = &Sbus;

    pv_ = IntVect(pv);
    pq_ = IntVect(pq);

    nb_pv_ = static_cast<int>(pv_.size());
    nb_pq_ = static_cast<int>(pq_.size());
    lag_   = 0;

    pvpq_.resize(nb_pv_ + nb_pq_);
    pvpq_ << pv_, pq_;

    nb_pvpq_ = static_cast<int>(pvpq_.size());
    const int n_bus  = static_cast<int>(Ybus.rows());

    pvpq_inv_.assign(n_bus, -1);
    for (int i = 0; i < nb_pvpq_; ++i) pvpq_inv_[pvpq_(i)] = i;
    pq_inv_.assign(n_bus, -1);
    for (int i = 0; i < nb_pq_; ++i) pq_inv_[pq_(i)] = i;

    // init the sparsity pattern
    // I think we don't really care about the
    // values
    dS_dVm_ = Ybus;
    dS_dVa_ = Ybus;
    
    need_full_rebuild_ = true;
}

// ---- Phase 1.5: per-compute_pf state update ----------------------------------

inline void NRSystem<>::update_state(
    const Eigen::SparseMatrix<cplx_type>& Ybus,
    const CplxVect&                        V_init,
    const CplxVect&                        Sbus)
{
    Ybus_ptr_ = &Ybus;
    Sbus_ptr_ = &Sbus;

    Va_ = V_init.array().arg();
    Vm_ = V_init.array().abs();
    V_  = V_init;
}

// ---- Phase 2: build J sparsity + value_map -----------------------------------

inline void NRSystem<>::build_J_sparsity()
{
    // reset J
    J_         = Eigen::SparseMatrix<real_type>();

    // compute its sparsity pattern
    const Eigen::SparseMatrix<cplx_type> & Ybus = *Ybus_ptr_;
    const int n_bus  = Ybus.rows();
    const int dim_J  = nb_pvpq_ + nb_pq_;
    const int nnz_Y  = Ybus.nonZeros();

    std::vector<Contrib> c11, c12, c21, c22;
    // c11.reserve()  // TODO !
    int k = 0;
    for (int outer = 0; outer < Ybus.outerSize(); ++outer) {
        for (Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>::InnerIterator
             it(Ybus, outer); it; ++it, ++k)
        {
            int i = (int)it.row(), j = (int)it.col();
            int ri = pvpq_inv_[i], rq = pq_inv_[i];
            int ci = pvpq_inv_[j], cq = pq_inv_[j];
            if (ri >= 0 && ci >= 0) c11.push_back({ri,          ci,          k});
            if (ri >= 0 && cq >= 0) c12.push_back({ri,          nb_pvpq_ + cq, k});
            if (rq >= 0 && ci >= 0) c21.push_back({nb_pvpq_ + rq, ci,          k});
            if (rq >= 0 && cq >= 0) c22.push_back({nb_pvpq_ + rq, nb_pvpq_ + cq, k});
        }
    }

    std::vector<Eigen::Triplet<real_type> > triplets;
    triplets.reserve(c11.size() + c12.size() + c21.size() + c22.size());
    for (auto& c : c11) triplets.push_back({c.jrow, c.jcol, 0.});
    for (auto& c : c12) triplets.push_back({c.jrow, c.jcol, 0.});
    for (auto& c : c21) triplets.push_back({c.jrow, c.jcol, 0.});
    for (auto& c : c22) triplets.push_back({c.jrow, c.jcol, 0.});

    J_.resize(dim_J, dim_J);
    J_.setFromTriplets(triplets.begin(), triplets.end());
    J_.makeCompressed();
    // std::cout << "J.shape " << J_.rows() << " ," << J_.cols() << "\n";
    // std::cout << "J.nnz " << J_.nonZeros() << "\n";
    _build_value_map(c11, c21, c12, c22);
}

inline void NRSystem<>::_build_value_map(
    const std::vector<Contrib> & c11,
    const std::vector<Contrib> & c21,
    const std::vector<Contrib> & c12,
    const std::vector<Contrib> & c22
)
{
    const int nnz_Y  = Ybus_ptr_->nonZeros();

    // now synch the map_j{1,2}{1,2} with the new J_ matrix
    map_j11_.assign(nnz_Y, -1);
    map_j12_.assign(nnz_Y, -1);
    map_j21_.assign(nnz_Y, -1);
    map_j22_.assign(nnz_Y, -1);

    for (auto& c : c11) map_j11_[c.ybus_k] = find_J_pos(c.jrow, c.jcol);
    for (auto& c : c12) map_j12_[c.ybus_k] = find_J_pos(c.jrow, c.jcol);
    for (auto& c : c21) map_j21_[c.ybus_k] = find_J_pos(c.jrow, c.jcol);
    for (auto& c : c22) map_j22_[c.ybus_k] = find_J_pos(c.jrow, c.jcol);
    // std::cout << "J11: \n";
    // int i = 0;
    // for(auto & c : map_j11_)
    // {
    //     std::cout << "\t map_j11_[" << i << "] =" << c << std::endl;
    //     i++;
    // }
    // std::cout << "-------------------------------\n";
    need_full_rebuild_ = false;
}

// ---- Phase 3: fill J values (fast, called every factorisation) ---------------

inline void NRSystem<>::fill_J()
{
    auto timer = CustTimer();

    std::tuple<int, int> ref_ = {-1, -1};
    std::vector<std::tuple<int, int> > map_pos(J_.nonZeros(), ref_);
    const cplx_type* ds_dvm = dS_dVm_.valuePtr();
    const cplx_type* ds_dva = dS_dVa_.valuePtr();
    size_t i = 0;
    for(auto & c : map_j11_){
        if(c == -1){
            // coeff of J11 not used in J
            i++;
            continue;
        }
        // if(map_pos[c] != ref_) throw std::runtime_error("Error in NRSystem::fill_J (J11)");
        J_.valuePtr()[c] = std::real(ds_dva[i]);
        // map_pos[c] = {0, i};
        i++;
    }

    i = 0;
    for(auto & c : map_j21_){
        if(c == -1){
            // coeff of J21 not used in J
            i++;
            continue;
        }
        // if(map_pos[c] != ref_) throw std::runtime_error("Error in NRSystem::fill_J (J21)");
        J_.valuePtr()[c] = std::imag(ds_dva[i]);
        // map_pos[c] = {1, i};
        i++;
    }

    i = 0;
    for(auto & c : map_j12_){
        if(c == -1){
            // coeff of J12 not used in J
            i++;
            continue;
        }
        // if(map_pos[c] != ref_) throw std::runtime_error("Error in NRSystem::fill_J (J12)");
        J_.valuePtr()[c] = std::real(ds_dvm[i]);
        // map_pos[c] = {2, i};
        i++;
    }

    i = 0;
    for(auto & c : map_j22_){
        if(c == -1){
            // coeff of J22 not used in J
            i++;
            continue;
        }
        // if(map_pos[c] != ref_) throw std::runtime_error("Error in NRSystem::fill_J (J22)");
        J_.valuePtr()[c] = std::imag(ds_dvm[i]);
        // map_pos[c] = {3, i};
        i++;
    }
    // std::cout << "map_pos: \n";
    // i = 0;
    // for(auto & c : map_pos)
    // {
    //     std::cout << "\t map_pos[" << i << "] = " << std::get<0>(c) << " " << std::get<1>(c) << std::endl;
    //     i++;
    // }
    // std::cout << "-------------------------------\n";
    timer_fillJ_ += timer.duration();
}


inline void NRSystem<>::fill_internal_variables()
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

inline RealVect NRSystem<>::mismatch() const
{
    return _mismatch_core(V_);
}

inline void NRSystem<>::apply_step(const RealVect& dx)
{
    if (nb_pv_ > 0)
        Va_(pv_) += theta(dx).segment(0, nb_pv_);
    if (nb_pq_ > 0) {
        Va_(pq_) += theta(dx).segment(nb_pv_, nb_pq_);
        Vm_(pq_) += vm(dx);
    }
    V_ = _reconstruct_V(Va_, Vm_);
    if (Vm_.minCoeff() < static_cast<real_type>(0.)) {
        Vm_ = V_.array().abs();
        Va_ = V_.array().arg();
    }
}

inline real_type NRSystem<>::mismatch_sq_norm_at(const RealVect& dx) const
{
    return _mismatch_core(_compute_trial_V(dx)).squaredNorm();
}

// ---- Virtual hooks -----------------------------------------------------------

// inline void NRSystem<>::_collect_J_triplets(
//     std::vector<Eigen::Triplet<double>>& coeffs) const
// {
//     const int n_pvpq = theta_size();

//     const Eigen::SparseMatrix<real_type> dS_dVa_r = dS_dVa_.real();
//     const Eigen::SparseMatrix<real_type> dS_dVa_i = dS_dVa_.imag();
//     const Eigen::SparseMatrix<real_type> dS_dVm_r = dS_dVm_.real();
//     const Eigen::SparseMatrix<real_type> dS_dVm_i = dS_dVm_.imag();

//     int nb_obj = 0;
//     std::vector<Eigen::Index> inner_index;
//     std::vector<real_type> values;

//     // dTheta columns (0..n_pvpq-1)
//     for (int col_id = 0; col_id < n_pvpq; ++col_id) {
//         nb_obj = 0; inner_index.clear(); values.clear();
//         _get_values_J(nb_obj, inner_index, values, dS_dVa_r, pvpq_inv_, pvpq_,
//                       static_cast<size_t>(col_id), 0, 0);
//         _get_values_J(nb_obj, inner_index, values, dS_dVa_i, pq_inv_, pvpq_,
//                       static_cast<size_t>(col_id), static_cast<size_t>(n_pvpq), 0);
//         for (int k = 0; k < nb_obj; ++k)
//             coeffs.push_back(Eigen::Triplet<double>(
//                 static_cast<int>(inner_index[k]), col_id, values[k]));
//     }

//     // dVm columns (n_pvpq..n_pvpq+n_pq-1)
//     for (int col_id = 0; col_id < nb_pq_; ++col_id) {
//         nb_obj = 0; inner_index.clear(); values.clear();
//         _get_values_J(nb_obj, inner_index, values, dS_dVm_r, pvpq_inv_, pq_,
//                       static_cast<size_t>(col_id), 0, 0);
//         _get_values_J(nb_obj, inner_index, values, dS_dVm_i, pq_inv_, pq_,
//                       static_cast<size_t>(col_id), static_cast<size_t>(n_pvpq), 0);
//         for (int k = 0; k < nb_obj; ++k)
//             coeffs.push_back(Eigen::Triplet<double>(
//                 static_cast<int>(inner_index[k]), col_id + n_pvpq, values[k]));
//     }
// }

// inline cplx_type* NRSystem<>::_get_entry_ptr(int row, int col)
// {
//     const int n_pvpq = theta_size();
//     if (col < n_pvpq) {
//         if (row < n_pvpq)
//             return &dS_dVa_.coeffRef(pvpq_[row],            pvpq_[col]);
//         else
//             return &dS_dVa_.coeffRef(pq_[row - n_pvpq],     pvpq_[col]);
//     } else {
//         const int pq_col = col - n_pvpq;
//         if (row < n_pvpq)
//             return &dS_dVm_.coeffRef(pvpq_[row],            pq_[pq_col]);
//         else
//             return &dS_dVm_.coeffRef(pq_[row - n_pvpq],     pq_[pq_col]);
//     }
// }

// ---- Shared helpers (non-virtual) --------------------------------------------

// inline void NRSystem<>::_collect_value_map()
// {
//     const int end_col = J_.cols() - lag_;
//     value_map_.clear();
//     value_map_.reserve(static_cast<size_t>(J_.nonZeros()));

//     for (int col = 0; col < end_col; ++col) {
//         for (Eigen::SparseMatrix<real_type>::InnerIterator it(J_, col); it; ++it) {
//             value_map_.push_back(_get_entry_ptr(static_cast<int>(it.row()), col));
//         }
//     }
//     dS_dVa_.makeCompressed();
//     dS_dVm_.makeCompressed();
// }

// inline void NRSystem<>::_dSbus_dV(
//     const Eigen::SparseMatrix<cplx_type>& Ybus, const CplxVect& V)
// {
//     auto timer = CustTimer();
//     const auto size_dS = V.size();
//     const CplxVect Vnorm = V.array() / V.array().abs();
//     const CplxVect Ibus  = Ybus * V;
//     const CplxVect conjIbus_Vnorm = Ibus.array().conjugate() * Vnorm.array();

//     if(need_full_rebuild_){
//         // I had to rebuild the system
//         dS_dVm_ = Ybus;
//         dS_dVa_ = Ybus;
//         // TODO init from Ybus sparsity pattern with all 0 instead
//     }

//     cplx_type* ds_dvm_x = dS_dVm_.valuePtr();
//     cplx_type* ds_dva_x = dS_dVa_.valuePtr();

//     unsigned int pos = 0;
//     for (int col_id = 0; col_id < size_dS; ++col_id) {
//         for (Eigen::SparseMatrix<cplx_type>::InnerIterator it(Ybus, col_id); it; ++it) {
//             const int row_id = static_cast<int>(it.row());
//             const cplx_type el = it.value();

//             cplx_type& dvm = ds_dvm_x[pos];
//             cplx_type& dva = ds_dva_x[pos];

//             dvm = el * Vnorm(col_id);
//             dvm = std::conj(dvm) * V(row_id);

//             dva = el * V(col_id);
//             if (col_id == row_id) {
//                 dvm += conjIbus_Vnorm(row_id);
//                 ds_dva_x[pos] -= Ibus(row_id);
//             }
//             const cplx_type tmp = BaseConstants::my_i * V(row_id);
//             dva = std::conj(-dva) * tmp;
//             ++pos;
//         }
//     }
//     timer_dSbus_ += timer.duration();
// }

inline CplxVect NRSystem<>::_reconstruct_V(const RealVect& Va, const RealVect& Vm)
{
    const cplx_type m_i = BaseConstants::my_i;
    return Vm.array() * (Va.array().cos().template cast<cplx_type>()
                         + m_i * Va.array().sin().template cast<cplx_type>());
}

inline CplxVect NRSystem<>::_compute_trial_V(const RealVect& dx) const
{
    RealVect Va_t = Va_;
    RealVect Vm_t = Vm_;
    if (nb_pv_ > 0) Va_t(pv_) += theta(dx).segment(0, nb_pv_);
    if (nb_pq_ > 0) {
        Va_t(pq_) += theta(dx).segment(nb_pv_, nb_pq_);
        Vm_t(pq_) += vm(dx);
    }
    return _reconstruct_V(Va_t, Vm_t);
}

inline RealVect NRSystem<>::_mismatch_core(const CplxVect& V_trial) const
{
    auto mis = V_trial.array() * (*Ybus_ptr_ * V_trial).array().conjugate()
               - Sbus_ptr_->array();
    const RealVect real_ = mis.real();
    const RealVect imag_ = mis.imag();

    RealVect res(nb_pv_ + 2 * nb_pq_);
    res.segment(0,               nb_pv_) = -real_(pv_);
    res.segment(nb_pv_,          nb_pq_) = -real_(pq_);
    res.segment(nb_pv_ + nb_pq_, nb_pq_) = -imag_(pq_);
    return res;
}

// inline void NRSystem<>::_get_values_J(
//     int& nb_obj_this_col,
//     std::vector<Eigen::Index>& inner_index,
//     std::vector<real_type>& values,
//     const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
//     const std::vector<int>& index_row_inv,
//     const Eigen::VectorXi& index_col,
//     size_t col_id,
//     size_t row_lag,
//     size_t col_lag) const
// {
//     const int col_id_mat = index_col(static_cast<int>(col_id + col_lag));
//     _get_values_J(nb_obj_this_col, inner_index, values, mat, index_row_inv,
//                   static_cast<size_t>(col_id_mat), row_lag, col_lag);
// }

// inline void NRSystem<>::_get_values_J(
//     int& nb_obj_this_col,
//     std::vector<Eigen::Index>& inner_index,
//     std::vector<real_type>& values,
//     const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
//     const std::vector<int>& index_row_inv,
//     size_t col_id_mat,
//     size_t row_lag,
//     size_t /*col_lag*/) const
// {
//     const int start_id = mat.outerIndexPtr()[col_id_mat];
//     const int end_id   = mat.outerIndexPtr()[col_id_mat + 1];
//     const real_type* val_ptr = mat.valuePtr();
//     for (int obj_id = start_id; obj_id < end_id; ++obj_id) {
//         const int row_id_dS = mat.innerIndexPtr()[obj_id];
//         const int row_id    = index_row_inv[row_id_dS];
//         if (row_id >= 0) {
//             inner_index.push_back(static_cast<Eigen::Index>(row_id) + row_lag);
//             values.push_back(val_ptr[obj_id]);
//             ++nb_obj_this_col;
//         }
//     }
// }

// =============================================================================
//  NRSystem<MultiSlack, Rest...> — distributed-slack extension
// =============================================================================

// ---- Phase 1: topology init --------------------------------------------------

template <typename... Rest>
void NRSystem<MultiSlack, Rest...>::init_topology(
    const Eigen::SparseMatrix<cplx_type>& Ybus,
    const CplxVect&                        Sbus,
    Eigen::Ref<const IntVect>              slack_ids,
    const RealVect&                        slack_weights,
    Eigen::Ref<const IntVect>              pv,
    Eigen::Ref<const IntVect>              pq)
{
    // std::cout << "NRSystem<MultiSlack>::init_topology \n";
    slack_bus_id_  = static_cast<size_t>(slack_ids(0));
    slack_weights_ = slack_weights;

    // Build my_pv = (extra slack buses beyond the primary) ++ pv
    const int nb_slack_added = static_cast<int>(slack_ids.size()) - 1;
    IntVect my_pv;
    if (nb_slack_added > 0) {
        my_pv.resize(static_cast<int>(pv.size()) + nb_slack_added);
        for (int i = 0; i < nb_slack_added; ++i) my_pv(i) = slack_ids(i + 1);
        for (int i = 0; i < static_cast<int>(pv.size()); ++i)
            my_pv(i + nb_slack_added) = pv(i);
    } else {
        my_pv = IntVect(pv);
    }

    Base::init_topology(Ybus, Sbus, slack_ids, slack_weights, my_pv, pq);
    this->lag_ += 1;   // MultiSlack always contributes exactly one extra row/col
}

// ---- Phase 1.5: per-compute_pf state update ----------------------------------

template <typename... Rest>
void NRSystem<MultiSlack, Rest...>::update_state(
    const Eigen::SparseMatrix<cplx_type>& Ybus,
    const CplxVect&                        V_init,
    const CplxVect&                        Sbus)
{
    Base::update_state(Ybus, V_init, Sbus);
    slack_absorbed_ = std::real(Sbus.sum());
}

// ---- NR primitives -----------------------------------------------------------

template <typename... Rest>
RealVect NRSystem<MultiSlack, Rest...>::mismatch() const
{
    return _mismatch_with_slack(this->V_, slack_absorbed_);
}

template <typename... Rest>
void NRSystem<MultiSlack, Rest...>::apply_step(const RealVect& dx)
{
    Base::apply_step(dx);
    slack_absorbed_ += dx(_J_slack_row());
}

template <typename... Rest>
real_type NRSystem<MultiSlack, Rest...>::mismatch_sq_norm_at(const RealVect& dx) const
{
    const real_type sa_trial = slack_absorbed_ + dx(_J_slack_row());
    return _mismatch_with_slack(this->_compute_trial_V(dx), sa_trial).squaredNorm();
}

// ---- Virtual hooks -----------------------------------------------------------

// template <typename... Rest>
// void NRSystem<MultiSlack, Rest...>::_collect_J_triplets(
//     std::vector<Eigen::Triplet<double>>& coeffs) const
// {
//     Base::_collect_J_triplets(coeffs);
//     _append_slack_triplets(coeffs);
// }

// template <typename... Rest>
// cplx_type* NRSystem<MultiSlack, Rest...>::_get_entry_ptr(int row, int col)
// {
//     if (row == _J_slack_row()) return nullptr;   // slack row is frozen
//     return Base::_get_entry_ptr(row, col);
// }

// ---- Private helpers ---------------------------------------------------------

template <typename... Rest>
void NRSystem<MultiSlack, Rest...>::_append_slack_triplets(
    std::vector<Eigen::Triplet<double>>& coeffs) const
{
    const int slack_row = _J_slack_row();
    const int slack_col = this->total() - 1;     // last column
    const int n_pvpq    = this->theta_size();

    const Eigen::SparseMatrix<real_type> dS_dVa_r = this->dS_dVa_.real();
    const Eigen::SparseMatrix<real_type> dS_dVm_r = this->dS_dVm_.real();
    const int nb_bus = dS_dVa_r.cols();

    // Slack bus row: dP_slack / dTheta  (dTheta columns 0..n_pvpq-1)
    for (int col_id = 0; col_id < nb_bus; ++col_id) {
        const int J_col = this->pvpq_inv_[col_id];
        if (J_col < 0) continue;
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(dS_dVa_r, col_id); it; ++it) {
            if (it.row() != static_cast<Eigen::Index>(slack_bus_id_)) continue;
            coeffs.push_back(Eigen::Triplet<double>(slack_row, J_col, it.value()));
        }
    }
    // Slack bus row: dP_slack / dVm  (dVm columns n_pvpq..n_pvpq+n_pq-1)
    for (int col_id = 0; col_id < nb_bus; ++col_id) {
        const int J_col = this->pq_inv_[col_id];
        if (J_col < 0) continue;
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(dS_dVm_r, col_id); it; ++it) {
            if (it.row() != static_cast<Eigen::Index>(slack_bus_id_)) continue;
            coeffs.push_back(Eigen::Triplet<double>(slack_row, J_col + n_pvpq, it.value()));
        }
    }

    // Slack column: diagonal entry (slack_row, slack_col)
    coeffs.push_back(Eigen::Triplet<double>(
        slack_row, slack_col,
        slack_weights_[static_cast<int>(slack_bus_id_)]));

    // Slack column: coupling to pvpq rows (rows 0..n_pvpq-1)
    for (int i = 0; i < static_cast<int>(this->pvpq_.size()); ++i) {
        const real_type sl_w = slack_weights_(this->pvpq_(i));
        if (std::abs(sl_w) > BaseConstants::_tol_equal_float)
            coeffs.push_back(Eigen::Triplet<double>(i, slack_col, sl_w));
    }
}

template <typename... Rest>
RealVect NRSystem<MultiSlack, Rest...>::_mismatch_with_slack(
    const CplxVect& V_trial, real_type sa) const
{
    const int n_pv = this->nb_pv_;
    const int n_pq = this->nb_pq_;

    auto mis = V_trial.array() * (*this->Ybus_ptr_ * V_trial).array().conjugate()
               - this->Sbus_ptr_->array()
               + sa * slack_weights_.array();
    const RealVect real_ = mis.real();
    const RealVect imag_ = mis.imag();

    // Layout: [dP_pv, dP_pq, dQ_pq, slack_eq]  (negated)
    RealVect res(n_pv + 2 * n_pq + 1);
    res.segment(0,               n_pv) = -real_(this->pv_);
    res.segment(n_pv,            n_pq) = -real_(this->pq_);
    res.segment(n_pv + n_pq,     n_pq) = -imag_(this->pq_);
    res(n_pv + 2 * n_pq) = -real_(static_cast<int>(slack_bus_id_));
    return res;
}

} // namespace ls2g
