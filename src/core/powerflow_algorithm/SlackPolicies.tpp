// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// ============================================================
//  MultiSlackPolicy
// ============================================================

template<class NRAlgoT>
Eigen::VectorXi MultiSlackPolicy::get_my_pv(NRAlgoT& algo,
                                              Eigen::Ref<const IntVect> slack_ids,
                                              Eigen::Ref<const IntVect> pv)
{
    return algo.retrieve_pv_with_slack(slack_ids, pv);
}

inline real_type MultiSlackPolicy::initial_slack_absorbed(const CplxVect& Sbus)
{
    return std::real(Sbus.sum());
}

template<class NRAlgoT>
RealVect MultiSlackPolicy::evaluate_Fx(NRAlgoT& algo,
                                        const Eigen::SparseMatrix<cplx_type>& Ybus,
                                        const CplxVect& V,
                                        const CplxVect& Sbus,
                                        size_t slack_bus_id,
                                        real_type slack_absorbed,
                                        const RealVect& slack_weights,
                                        const Eigen::VectorXi& my_pv,
                                        Eigen::Ref<const IntVect> pq)
{
    return algo._evaluate_Fx(Ybus, V, Sbus, slack_bus_id, slack_absorbed, slack_weights, my_pv, pq);
}

template<class NRAlgoT>
void MultiSlackPolicy::update_slack_absorbed(NRAlgoT& algo,
                                              const RealVect& dx,
                                              real_type& slack_absorbed)
{
    slack_absorbed -= dx(algo._layout.J_slack_row());
}

template<class NRAlgoT>
void MultiSlackPolicy::fill_jacobian_matrix(NRAlgoT& algo,
                                             const Eigen::SparseMatrix<cplx_type>& Ybus,
                                             const CplxVect& V,
                                             size_t slack_bus_id,
                                             const RealVect& slack_weights,
                                             const Eigen::VectorXi& pq,
                                             const Eigen::VectorXi& pvpq,
                                             const std::vector<int>& pq_inv,
                                             const std::vector<int>& pvpq_inv)
{
    algo._dSbus_dV(Ybus, V);

    auto timer = CustTimer();
    const auto n_pvpq = pvpq.size();
    const auto n_pq = pq.size();
    const auto size_j = n_pvpq + n_pq + lag;  // lag=1 for multi-slack

    if(algo.J_.cols() != (Eigen::Index)size_j)
    {
        fill_jac_unknown_sparsity(algo, Ybus, V, slack_bus_id, slack_weights, pq, pvpq, pq_inv, pvpq_inv);
        fill_value_map_impl(algo, slack_bus_id, pq, pvpq, false);
    }
    else
    {
        if(algo.value_map_.size() == 0) fill_value_map_impl(algo, slack_bus_id, pq, pvpq, true);
        fill_jac_known_sparsity(algo, slack_bus_id, pq, pvpq);
    }
    algo.timer_fillJ_ += timer.duration();
}

template<class NRAlgoT>
void MultiSlackPolicy::fill_jac_unknown_sparsity(NRAlgoT& algo,
                                                   const Eigen::SparseMatrix<cplx_type>& Ybus,
                                                   const CplxVect& V,
                                                   size_t slack_bus_id,
                                                   const RealVect& slack_weights,
                                                   const Eigen::VectorXi& pq,
                                                   const Eigen::VectorXi& pvpq,
                                                   const std::vector<int>& pq_inv,
                                                   const std::vector<int>& pvpq_inv)
{
    using StorageIndex = Eigen::SparseMatrix<cplx_type>::StorageIndex;

    const size_t n_pvpq = pvpq.size();
    const size_t n_pq = pq.size();
    const auto size_j = n_pvpq + n_pq + lag;  // lag=1

    const Eigen::SparseMatrix<real_type> dS_dVa_r = algo.dS_dVa_.real();
    const Eigen::SparseMatrix<real_type> dS_dVa_i = algo.dS_dVa_.imag();
    const Eigen::SparseMatrix<real_type> dS_dVm_r = algo.dS_dVm_.real();
    const Eigen::SparseMatrix<real_type> dS_dVm_i = algo.dS_dVm_.imag();

    if(algo.J_.cols() != (Eigen::Index)size_j)
        algo.J_ = Eigen::SparseMatrix<real_type>(size_j, size_j);

    std::vector<Eigen::Triplet<double>> coeffs;
    coeffs.reserve(2*(algo.dS_dVa_.nonZeros() + algo.dS_dVm_.nonZeros()) + slack_weights.size());

    int nb_obj_this_col = 0;
    std::vector<Eigen::Index> inner_index;
    std::vector<real_type> values;

    for(Eigen::Index col_id=0; col_id < (Eigen::Index)n_pvpq; ++col_id){
        nb_obj_this_col = 0;
        inner_index.clear();
        values.clear();

        algo._get_values_J(nb_obj_this_col, inner_index, values,
                           dS_dVa_r,
                           pvpq_inv, pvpq,
                           col_id,
                           static_cast<size_t>(algo._layout.J_dP_row()),
                           0);
        algo._get_values_J(nb_obj_this_col, inner_index, values,
                           dS_dVa_i,
                           pq_inv, pvpq,
                           col_id,
                           static_cast<size_t>(algo._layout.J_dQ_row()),
                           0);

        for(Eigen::Index in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
            StorageIndex row_id = static_cast<StorageIndex>(inner_index[in_ind]);
            coeffs.push_back(Eigen::Triplet<double>(row_id,
                             static_cast<StorageIndex>(col_id) + static_cast<StorageIndex>(algo._layout.lag()),
                             values[in_ind]));
        }
    }

    for(Eigen::Index col_id=0; col_id < (Eigen::Index)n_pq; ++col_id){
        nb_obj_this_col = 0;
        inner_index.clear();
        values.clear();

        algo._get_values_J(nb_obj_this_col, inner_index, values,
                           dS_dVm_r,
                           pvpq_inv, pq,
                           col_id,
                           static_cast<size_t>(algo._layout.J_dP_row()),
                           0);
        algo._get_values_J(nb_obj_this_col, inner_index, values,
                           dS_dVm_i,
                           pq_inv, pq,
                           col_id,
                           static_cast<size_t>(algo._layout.J_dQ_row()),
                           0);

        for(Eigen::Index in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
            auto row_id = static_cast<StorageIndex>(inner_index[in_ind]);
            coeffs.push_back(Eigen::Triplet<double>(row_id,
                             static_cast<StorageIndex>(col_id + n_pvpq) + static_cast<StorageIndex>(algo._layout.lag()),
                             values[in_ind]));
        }
    }

    // slack bus row: equation for the reference slack bus
    const int nb_bus = dS_dVa_r.cols();
    for(Eigen::Index col_id=0; col_id < nb_bus; ++col_id){
        const auto J_col = pvpq_inv[col_id];
        if(J_col < 0) continue;
        for(Eigen::SparseMatrix<real_type>::InnerIterator it(dS_dVa_r, col_id); it; ++it)
        {
            if(it.row() != (Eigen::Index)slack_bus_id) continue;
            coeffs.push_back(Eigen::Triplet<double>(algo._layout.J_slack_row(),
                             J_col + algo._layout.lag(),
                             it.value()));
        }
    }
    for(Eigen::Index col_id=0; col_id < nb_bus; ++col_id){
        const auto J_col = pq_inv[col_id];
        if(J_col < 0) continue;
        for(Eigen::SparseMatrix<real_type>::InnerIterator it(dS_dVm_r, col_id); it; ++it)
        {
            if(it.row() != (Eigen::Index)slack_bus_id) continue;
            coeffs.push_back(Eigen::Triplet<double>(algo._layout.J_slack_row(),
                             static_cast<StorageIndex>(J_col + n_pvpq) + static_cast<StorageIndex>(algo._layout.lag()),
                             it.value()));
        }
    }

    // slack column: distributed-slack equation
    const StorageIndex last_col = static_cast<StorageIndex>(algo._layout.J_slack_col());
    coeffs.push_back(Eigen::Triplet<double>(algo._layout.J_slack_row(), last_col, slack_weights[slack_bus_id]));
    auto row_j = algo._layout.lag();
    for(auto ind : pvpq){
        auto sl_w = slack_weights(ind);
        if(std::abs(sl_w) > BaseConstants::_tol_equal_float)
            coeffs.push_back(Eigen::Triplet<double>(row_j, last_col, sl_w));
        ++row_j;
    }
    algo.J_.setFromTriplets(coeffs.begin(), coeffs.end());
    algo.J_.makeCompressed();
}

template<class NRAlgoT>
void MultiSlackPolicy::fill_value_map_impl(NRAlgoT& algo,
                                            size_t slack_bus_id,
                                            const Eigen::VectorXi& pq,
                                            const Eigen::VectorXi& pvpq,
                                            bool reset_J)
{
    const int n_pvpq = static_cast<int>(pvpq.size());
    algo.value_map_.clear();
    algo.value_map_.reserve(algo.J_.nonZeros());

    const auto n_row = algo.J_.cols();
    for(Eigen::Index col_=static_cast<Eigen::Index>(algo._layout.lag()); col_ < n_row; ++col_){
        for(Eigen::SparseMatrix<real_type>::InnerIterator it(algo.J_, col_); it; ++it)
        {
            auto row_id = it.row();
            const auto col_id = it.col() - static_cast<Eigen::Index>(algo._layout.lag());
            if(reset_J) it.valueRef() = 0.;

            if(row_id == static_cast<Eigen::Index>(algo._layout.J_slack_row())){
                // slack bus row
                const size_t row_id_dS_dVx_r = slack_bus_id;
                if(col_id < n_pvpq){
                    const int col_id_dS_dVa_r = pvpq[col_id];
                    algo.value_map_.push_back(&algo.dS_dVa_.coeffRef(row_id_dS_dVx_r, col_id_dS_dVa_r));
                } else {
                    const int col_id_dS_dVm_r = pq[col_id - n_pvpq];
                    algo.value_map_.push_back(&algo.dS_dVm_.coeffRef(row_id_dS_dVx_r, col_id_dS_dVm_r));
                }
            } else {
                row_id -= static_cast<Eigen::Index>(algo._layout.lag());
                if((col_id < n_pvpq) && (row_id < n_pvpq)){
                    // J11 (dS_dVa_r)
                    const int row_id_dS = pvpq[row_id];
                    const int col_id_dS = pvpq[col_id];
                    algo.value_map_.push_back(&algo.dS_dVa_.coeffRef(row_id_dS, col_id_dS));
                } else if((col_id < n_pvpq) && (row_id >= n_pvpq)){
                    // J21 (dS_dVa_i)
                    const int row_id_dS = pq[row_id - n_pvpq];
                    const int col_id_dS = pvpq[col_id];
                    algo.value_map_.push_back(&algo.dS_dVa_.coeffRef(row_id_dS, col_id_dS));
                } else if((col_id >= n_pvpq) && (row_id < n_pvpq)){
                    // J12 (dS_dVm_r)
                    const int row_id_dS = pvpq[row_id];
                    const int col_id_dS = pq[col_id - n_pvpq];
                    algo.value_map_.push_back(&algo.dS_dVm_.coeffRef(row_id_dS, col_id_dS));
                } else {
                    // J22 (dS_dVm_i)
                    const int row_id_dS = pq[row_id - n_pvpq];
                    const int col_id_dS = pq[col_id - n_pvpq];
                    algo.value_map_.push_back(&algo.dS_dVm_.coeffRef(row_id_dS, col_id_dS));
                }
            }
        }
    }
    algo.dS_dVa_.makeCompressed();
    algo.dS_dVm_.makeCompressed();
}

template<class NRAlgoT>
void MultiSlackPolicy::fill_jac_known_sparsity(NRAlgoT& algo,
                                                size_t slack_bus_id,
                                                const Eigen::VectorXi& pq,
                                                const Eigen::VectorXi& pvpq)
{
    const Eigen::Index n_pvpq_1 = static_cast<Eigen::Index>(algo._layout.J_dQ_row());  // lag + theta_size
    const auto n_cols = algo.J_.cols();
    unsigned int pos_el = 0;
    for(Eigen::Index col_id=static_cast<Eigen::Index>(algo._layout.lag()); col_id < n_cols; ++col_id){
        for(Eigen::SparseMatrix<real_type>::InnerIterator it(algo.J_, col_id); it; ++it)
        {
            const auto row_id = it.row();
            it.valueRef() = row_id < n_pvpq_1
                            ? std::real(*algo.value_map_[pos_el])
                            : std::imag(*algo.value_map_[pos_el]);
            ++pos_el;
        }
    }
}


// ============================================================
//  SingleSlackPolicy
// ============================================================

template<class NRAlgoT>
Eigen::VectorXi SingleSlackPolicy::get_my_pv(NRAlgoT& algo,
                                               Eigen::Ref<const IntVect> slack_ids,
                                               Eigen::Ref<const IntVect> pv)
{
    return IntVect(pv);
}

inline real_type SingleSlackPolicy::initial_slack_absorbed(const CplxVect& Sbus)
{
    return static_cast<real_type>(0.);
}

template<class NRAlgoT>
RealVect SingleSlackPolicy::evaluate_Fx(NRAlgoT& algo,
                                         const Eigen::SparseMatrix<cplx_type>& Ybus,
                                         const CplxVect& V,
                                         const CplxVect& Sbus,
                                         size_t slack_bus_id,
                                         real_type slack_absorbed,
                                         const RealVect& slack_weights,
                                         const Eigen::VectorXi& my_pv,
                                         Eigen::Ref<const IntVect> pq)
{
    return algo._evaluate_Fx(Ybus, V, Sbus, my_pv, pq);
}

template<class NRAlgoT>
void SingleSlackPolicy::update_slack_absorbed(NRAlgoT& algo,
                                               const RealVect& dx,
                                               real_type& slack_absorbed)
{
    // single-slack: no distributed-slack variable in the state vector
}

template<class NRAlgoT>
void SingleSlackPolicy::fill_jacobian_matrix(NRAlgoT& algo,
                                              const Eigen::SparseMatrix<cplx_type>& Ybus,
                                              const CplxVect& V,
                                              size_t slack_bus_id,
                                              const RealVect& slack_weights,
                                              const Eigen::VectorXi& pq,
                                              const Eigen::VectorXi& pvpq,
                                              const std::vector<int>& pq_inv,
                                              const std::vector<int>& pvpq_inv)
{
    algo._dSbus_dV(Ybus, V);

    auto timer = CustTimer();
    const int n_pvpq = static_cast<int>(pvpq.size());
    const int n_pq = static_cast<int>(pq.size());
    const int size_j = n_pvpq + n_pq;  // lag=0 for single-slack

    if(algo.J_.cols() != size_j)
    {
        fill_jac_unknown_sparsity(algo, Ybus, V, pq, pvpq, pq_inv, pvpq_inv);
        fill_value_map_impl(algo, pq, pvpq, false);
    }
    else
    {
        if(algo.value_map_.size() == 0) fill_value_map_impl(algo, pq, pvpq, true);
        fill_jac_known_sparsity(algo, pq, pvpq);
    }
    algo.timer_fillJ_ += timer.duration();
}

template<class NRAlgoT>
void SingleSlackPolicy::fill_jac_unknown_sparsity(NRAlgoT& algo,
                                                    const Eigen::SparseMatrix<cplx_type>& Ybus,
                                                    const CplxVect& V,
                                                    const Eigen::VectorXi& pq,
                                                    const Eigen::VectorXi& pvpq,
                                                    const std::vector<int>& pq_inv,
                                                    const std::vector<int>& pvpq_inv)
{
    const int n_pvpq = static_cast<int>(pvpq.size());
    const int n_pq = static_cast<int>(pq.size());
    const int size_j = n_pvpq + n_pq;

    const Eigen::SparseMatrix<real_type> dS_dVa_r = algo.dS_dVa_.real();
    const Eigen::SparseMatrix<real_type> dS_dVa_i = algo.dS_dVa_.imag();
    const Eigen::SparseMatrix<real_type> dS_dVm_r = algo.dS_dVm_.real();
    const Eigen::SparseMatrix<real_type> dS_dVm_i = algo.dS_dVm_.imag();

    using StorageIndex = Eigen::SparseMatrix<cplx_type>::StorageIndex;

    if(algo.J_.cols() != size_j)
    {
        algo.J_ = Eigen::SparseMatrix<real_type>(size_j, size_j);
        algo.J_.reserve(2*(algo.dS_dVa_.nonZeros() + algo.dS_dVm_.nonZeros()));
    }

    std::vector<Eigen::Triplet<double>> coeffs;
    coeffs.reserve(2*(algo.dS_dVa_.nonZeros() + algo.dS_dVm_.nonZeros()));

    int nb_obj_this_col = 0;
    std::vector<Eigen::Index> inner_index;
    std::vector<real_type> values;

    for(int col_id=0; col_id < n_pvpq; ++col_id){
        nb_obj_this_col = 0;
        inner_index.clear();
        values.clear();

        algo._get_values_J(nb_obj_this_col, inner_index, values,
                           dS_dVa_r,
                           pvpq_inv, pvpq,
                           col_id,
                           static_cast<size_t>(algo._layout.J_dP_row()),
                           0);
        algo._get_values_J(nb_obj_this_col, inner_index, values,
                           dS_dVa_i,
                           pq_inv, pvpq,
                           col_id,
                           static_cast<size_t>(algo._layout.J_dQ_row()),
                           0);

        for(int in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
            int row_id = inner_index[in_ind];
            coeffs.push_back(Eigen::Triplet<double>(row_id, col_id, values[in_ind]));
        }
    }

    for(int col_id=0; col_id < n_pq; ++col_id){
        nb_obj_this_col = 0;
        inner_index.clear();
        values.clear();

        algo._get_values_J(nb_obj_this_col, inner_index, values,
                           dS_dVm_r,
                           pvpq_inv, pq,
                           col_id,
                           static_cast<size_t>(algo._layout.J_dP_row()),
                           0);
        algo._get_values_J(nb_obj_this_col, inner_index, values,
                           dS_dVm_i,
                           pq_inv, pq,
                           col_id,
                           static_cast<size_t>(algo._layout.J_dQ_row()),
                           0);

        for(int in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
            int row_id = inner_index[in_ind];
            coeffs.push_back(Eigen::Triplet<double>(row_id, col_id + n_pvpq, values[in_ind]));
        }
    }
    algo.J_.setFromTriplets(coeffs.begin(), coeffs.end());
    algo.J_.makeCompressed();
}

template<class NRAlgoT>
void SingleSlackPolicy::fill_value_map_impl(NRAlgoT& algo,
                                             const Eigen::VectorXi& pq,
                                             const Eigen::VectorXi& pvpq,
                                             bool reset_J)
{
    const int n_pvpq = static_cast<int>(pvpq.size());
    algo.value_map_.clear();
    algo.value_map_.reserve(algo.J_.nonZeros());

    const int n_col = static_cast<int>(algo.J_.cols());
    const int col_start = algo._layout.lag();  // 0 for single-slack
    for(int col_=col_start; col_ < n_col; ++col_){
        for(Eigen::SparseMatrix<real_type>::InnerIterator it(algo.J_, col_); it; ++it)
        {
            const int row_id = static_cast<int>(it.row());
            const int col_id = static_cast<int>(it.col());
            if(reset_J) it.valueRef() = 0.;

            if((col_id < n_pvpq) && (row_id < n_pvpq)){
                // J11 (dS_dVa_r)
                const int row_id_dS = pvpq[row_id];
                const int col_id_dS = pvpq[col_id];
                algo.value_map_.push_back(&algo.dS_dVa_.coeffRef(row_id_dS, col_id_dS));
            } else if((col_id < n_pvpq) && (row_id >= n_pvpq)){
                // J21 (dS_dVa_i)
                const int row_id_dS = pq[row_id - n_pvpq];
                const int col_id_dS = pvpq[col_id];
                algo.value_map_.push_back(&algo.dS_dVa_.coeffRef(row_id_dS, col_id_dS));
            } else if((col_id >= n_pvpq) && (row_id < n_pvpq)){
                // J12 (dS_dVm_r)
                const int row_id_dS = pvpq[row_id];
                const int col_id_dS = pq[col_id - n_pvpq];
                algo.value_map_.push_back(&algo.dS_dVm_.coeffRef(row_id_dS, col_id_dS));
            } else {
                // J22 (dS_dVm_i)
                const int row_id_dS = pq[row_id - n_pvpq];
                const int col_id_dS = pq[col_id - n_pvpq];
                algo.value_map_.push_back(&algo.dS_dVm_.coeffRef(row_id_dS, col_id_dS));
            }
        }
    }
}

template<class NRAlgoT>
void SingleSlackPolicy::fill_jac_known_sparsity(NRAlgoT& algo,
                                                  const Eigen::VectorXi& pq,
                                                  const Eigen::VectorXi& pvpq)
{
    const int n_pvpq_threshold = algo._layout.J_dQ_row();  // lag + theta_size = n_pvpq for lag=0
    const int n_cols = static_cast<int>(algo.J_.cols());
    unsigned int pos_el = 0;
    for(int col_id=algo._layout.lag(); col_id < n_cols; ++col_id){
        for(Eigen::SparseMatrix<real_type>::InnerIterator it(algo.J_, col_id); it; ++it)
        {
            const auto row_id = it.row();
            it.valueRef() = row_id < n_pvpq_threshold
                            ? std::real(*algo.value_map_[pos_el])
                            : std::imag(*algo.value_map_[pos_el]);
            ++pos_el;
        }
    }
}
