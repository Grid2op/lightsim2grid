// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "YbusPolicy.hpp"

namespace ls2g {

std::vector<Coeff> YbusPolicy::Contingency::_coeffs_for_branch_ids(
    const std::vector<int> & branch_ids,
    const LSGrid & grid_model,
    bool ac_solver_used,
    const SolverBusIdVect & id_me_to_solver,
    size_t n_line)
{
    const auto & powerlines = grid_model.get_powerlines_as_data();
    const auto & trafos = grid_model.get_trafos_as_data();
    std::vector<Coeff> res;
    res.reserve(branch_ids.size() * 4);  // usually there are 4 coeffs per powerlines / trafos
    int bus_1_id, bus_2_id;
    cplx_type y_ff, y_ft, y_tf, y_tt;
    bool status;
    for(auto br_id : branch_ids){
        int el_id;
        const BranchContainer *p_branch;
        if(static_cast<size_t>(br_id) < n_line)
        {
            // this is a powerline
            el_id = br_id;
            p_branch = & powerlines;
        }else{
            // this is a trafo
            el_id = br_id - static_cast<int>(n_line);
            p_branch = & trafos;
        }

        GlobalBusId glob_bus_1 = p_branch->get_bus_side_1(el_id);
        GlobalBusId glob_bus_2 = p_branch->get_bus_side_2(el_id);
        bus_1_id = glob_bus_1.cast_int() == GenericContainer::_deactivated_bus_id ?
            GenericContainer::_deactivated_bus_id :
            id_me_to_solver[glob_bus_1.cast_int()].cast_int();
        bus_2_id = glob_bus_2.cast_int() == GenericContainer::_deactivated_bus_id ?
            GenericContainer::_deactivated_bus_id :
            id_me_to_solver[glob_bus_2.cast_int()].cast_int();
        status = p_branch->get_status_global()[el_id];
        if(ac_solver_used){
            y_ff = p_branch->yac_eff_11()[el_id];
            y_ft = p_branch->yac_eff_12()[el_id];
            y_tf = p_branch->yac_eff_21()[el_id];
            y_tt = p_branch->yac_eff_22()[el_id];
        }else{
            y_ff = p_branch->ydc_11()[el_id];
            y_ft = p_branch->ydc_12()[el_id];
            y_tf = p_branch->ydc_21()[el_id];
            y_tt = p_branch->ydc_22()[el_id];
        }

        if(status)
        {
            // element is connected, update coeffs based on status of each powerlines
            if((bus_1_id != GenericContainer::_deactivated_bus_id)) res.push_back({bus_1_id, bus_1_id, y_ff});
            if((bus_2_id != GenericContainer::_deactivated_bus_id)) res.push_back({bus_2_id, bus_2_id, y_tt});
            if((bus_1_id != GenericContainer::_deactivated_bus_id) && (bus_2_id != GenericContainer::_deactivated_bus_id)){
                res.push_back({bus_1_id, bus_2_id, y_ft});
                res.push_back({bus_2_id, bus_1_id, y_tf});
            }
        }
    }
    return res;
}

void YbusPolicy::Contingency::init_li_coeffs(
    const LSGrid & grid_model,
    bool ac_solver_used,
    const SolverBusIdVect & id_me_to_solver,
    size_t n_line)
{
    li_coeffs.clear();
    li_coeffs.reserve(li_defaults.size());
    for(const auto & this_cont_id: li_defaults){
        const std::vector<int> branch_ids(this_cont_id.begin(), this_cont_id.end());
        li_coeffs.push_back(_coeffs_for_branch_ids(branch_ids, grid_model, ac_solver_used, id_me_to_solver, n_line));
    }
}

std::vector<int> YbusPolicy::Contingency::branch_ids_for_row(Eigen::Index row, size_t n_line) const
{
    std::vector<int> branch_ids;
    if(line_mask.rows() > 0){
        for(Eigen::Index col = 0; col < line_mask.cols(); ++col){
            if(line_mask(row, col)) branch_ids.push_back(static_cast<int>(col));
        }
    }
    if(trafo_mask.rows() > 0){
        for(Eigen::Index col = 0; col < trafo_mask.cols(); ++col){
            if(trafo_mask(row, col)) branch_ids.push_back(static_cast<int>(n_line + static_cast<size_t>(col)));
        }
    }
    return branch_ids;
}

void YbusPolicy::Contingency::init_li_coeffs_from_masks(
    const LSGrid & grid_model,
    bool ac_solver_used,
    const SolverBusIdVect & id_me_to_solver,
    size_t n_line,
    Eigen::Index nb_steps)
{
    li_coeffs.clear();
    li_coeffs.reserve(static_cast<size_t>(nb_steps));
    for(Eigen::Index row = 0; row < nb_steps; ++row){
        li_coeffs.push_back(_coeffs_for_branch_ids(branch_ids_for_row(row, n_line), grid_model, ac_solver_used, id_me_to_solver, n_line));
    }
}

void YbusPolicy::Contingency::remove_from_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                                               const std::vector<Coeff> & coeffs,
                                               bool ac_solver_used,
                                               AlgorithmSelector & algo)
{
    if(ac_solver_used)
    {
        for(const auto & coeff_to_remove: coeffs){
            Ybus.coeffRef(coeff_to_remove.row_id, coeff_to_remove.col_id) -= coeff_to_remove.value;
        }
    } else{
        // DC solver stores the ybus internally, I update it
        // instead of building it over and over
        for(const Coeff& coeff : coeffs){
            algo.update_internal_Ybus(coeff, false);  // false => remove the coeff (using -= )
        }
    }
}

void YbusPolicy::Contingency::readd_to_Ybus(
    Eigen::SparseMatrix<cplx_type> & Ybus,
    const std::vector<Coeff> & coeffs,
    bool ac_solver_used,
    AlgorithmSelector & algo)
{
    if(ac_solver_used){
        for(const Coeff & coeff_to_remove: coeffs){
            Ybus.coeffRef(coeff_to_remove.row_id, coeff_to_remove.col_id) += coeff_to_remove.value;
        }
    } else {
        // DC solver stores the ybus internally, I update it
        // instead of building it over and over
        for(const Coeff& coeff : coeffs){
            algo.update_internal_Ybus(coeff, true);  // true => add back the coeff (using += )
        }
    }
}

} // namespace ls2g
