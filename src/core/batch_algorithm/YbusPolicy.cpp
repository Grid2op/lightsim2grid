// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "YbusPolicy.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>

namespace ls2g {

std::vector<Coeff> YbusPolicy::Contingency::_coeffs_for_branch_ids(
    const std::vector<int> & branch_ids,
    const LSGrid & grid_model,
    bool ac_solver_used,
    const SolverBusIdVect & id_me_to_solver,
    size_t n_line,
    std::vector<std::pair<int, cplx_type> > & dc_shift_dp)
{
    const auto & powerlines = grid_model.get_powerlines_as_data();
    const auto & trafos = grid_model.get_trafos_as_data();
    Eigen::Ref<const RealVect> dc_x_tau_shift = trafos.dc_x_tau_shift();
    dc_shift_dp.clear();
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

        // DC: a branch open at one end is not in Bbus at all (see fillBdc: "disco on one
        // side == disco on both sides"), so disconnecting it has nothing to remove. AC
        // keeps its Kron-reduced side, which yac_eff_* already describes.
        if(!ac_solver_used &&
           (bus_1_id == GenericContainer::_deactivated_bus_id || bus_2_id == GenericContainer::_deactivated_bus_id)) continue;
        if(status)
        {
            // element is connected, update coeffs based on status of each powerlines
            if((bus_1_id != GenericContainer::_deactivated_bus_id)) res.push_back({bus_1_id, bus_1_id, y_ff});
            if((bus_2_id != GenericContainer::_deactivated_bus_id)) res.push_back({bus_2_id, bus_2_id, y_tt});
            if((bus_1_id != GenericContainer::_deactivated_bus_id) && (bus_2_id != GenericContainer::_deactivated_bus_id)){
                res.push_back({bus_1_id, bus_2_id, y_ft});
                res.push_back({bus_2_id, bus_1_id, y_tf});
                // a DC phase shifter: undo its injection pair, the exact opposite of
                // TrafoContainer::hack_Sbus_for_dc_phase_shifter (same test, same values)
                if(!ac_solver_used && p_branch == &trafos &&
                   std::abs(dc_x_tau_shift(el_id)) >= BaseConstants::_tol_equal_float){
                    dc_shift_dp.push_back({bus_1_id, cplx_type(dc_x_tau_shift(el_id), 0.)});
                    dc_shift_dp.push_back({bus_2_id, cplx_type(-dc_x_tau_shift(el_id), 0.)});
                }
            }
        }
    }
    return res;
}

namespace {
// the solver bus of a one-sided element that is in the base grid's matrix / injection,
// -1 for one off (status off, or on a disconnected bus)
template<class T>
int live_solver_bus(const T & container, int el_id, const SolverBusIdVect & id_me_to_solver)
{
    if(!container.get_status()[static_cast<size_t>(el_id)]) return -1;
    const int bus_me = container.get_bus_id()(el_id).cast_int();
    if(bus_me == BaseConstants::_deactivated_bus_id) return -1;
    const int bus_solver = id_me_to_solver[bus_me].cast_int();
    return bus_solver == BaseConstants::_deactivated_bus_id ? -1 : bus_solver;
}
// x / sn_mva the way the injection builders divide (only when sn_mva is not 1, so the
// values come out bit for bit the ones they stamped)
cplx_type per_unit(cplx_type x, real_type sn_mva)
{
    if(std::abs(sn_mva - 1.) > BaseConstants::_tol_equal_float) x /= sn_mva;
    return x;
}
}  // anonymous namespace

void YbusPolicy::Contingency::_append_shunt_trips(const std::vector<int> & shunt_ids,
                                                  const LSGrid & grid_model,
                                                  bool ac_solver_used,
                                                  const SolverBusIdVect & id_me_to_solver,
                                                  std::vector<Coeff> & coeffs,
                                                  std::vector<std::pair<int, cplx_type> > & sbus_dp)
{
    const ShuntContainer & shunts = grid_model.get_shunts();
    const real_type sn_mva = grid_model.get_sn_mva();
    const Eigen::Ref<const RealVect> target_p = shunts.get_target_p();
    const Eigen::Ref<const RealVect> target_q = shunts.get_target_q();
    for(int shunt_id : shunt_ids){
        const int bus = live_solver_bus(shunts, shunt_id, id_me_to_solver);
        if(bus < 0) continue;
        if(ac_solver_used){
            // the very coefficient ShuntContainer::_fillYbus stamps
            coeffs.push_back({bus, bus, per_unit(cplx_type(target_p(shunt_id), -target_q(shunt_id)), sn_mva)});
        } else {
            // DC: no matrix entry, its active part in Pbus (ShuntContainer::_fillSbus
            // subtracts it), given back
            sbus_dp.push_back({bus, per_unit(cplx_type(target_p(shunt_id), 0.), sn_mva)});
        }
    }
}

void YbusPolicy::Contingency::_append_load_trips(const std::vector<int> & load_ids,
                                                 const LSGrid & grid_model,
                                                 bool ac_solver_used,
                                                 const SolverBusIdVect & id_me_to_solver,
                                                 std::vector<std::pair<int, cplx_type> > & sbus_dp)
{
    const auto & loads = grid_model.get_loads_as_data();
    const real_type sn_mva = grid_model.get_sn_mva();
    const Eigen::Ref<const RealVect> target_p = loads.get_target_p();
    const Eigen::Ref<const RealVect> target_q = loads.get_target_q();
    for(int load_id : load_ids){
        const int bus = live_solver_bus(loads, load_id, id_me_to_solver);
        if(bus < 0) continue;
        // a load draws -(p + jq): given back
        const real_type q = ac_solver_used ? target_q(load_id) : 0.;
        sbus_dp.push_back({bus, per_unit(cplx_type(target_p(load_id), q), sn_mva)});
    }
}

void YbusPolicy::Contingency::_append_gen_trips(const std::vector<int> & gen_ids,
                                                const LSGrid & grid_model,
                                                bool ac_solver_used,
                                                const SolverBusIdVect & id_me_to_solver,
                                                std::vector<std::pair<int, cplx_type> > & sbus_dp)
{
    const auto & generators = grid_model.get_generators_as_data();
    const real_type sn_mva = grid_model.get_sn_mva();
    const Eigen::Ref<const RealVect> target_p = grid_model.get_gen_target_p();
    for(int gen_id : gen_ids){
        const int bus = live_solver_bus(generators, gen_id, id_me_to_solver);
        if(bus < 0) continue;
        // what SbusPolicy::Vary::fill_row takes back for a generator a row disconnects:
        // its active power, and the reactive setpoint of one that does not regulate
        real_type q = 0.;
        if(ac_solver_used && !generators.get_voltage_regulator_on(gen_id)) q = generators.get_target_q_mvar(gen_id);
        sbus_dp.push_back({bus, -per_unit(cplx_type(target_p(gen_id), q), sn_mva)});
    }
}

void YbusPolicy::Contingency::init_li_coeffs(
    const LSGrid & grid_model,
    bool ac_solver_used,
    const SolverBusIdVect & id_me_to_solver,
    size_t n_line)
{
    li_coeffs.clear();
    li_coeffs.reserve(li_defaults.size());
    li_sbus_dp.assign(li_defaults.size(), {});
    li_rows.assign(li_defaults.begin(), li_defaults.end());
    size_t cont_id = 0;
    for(const auto & this_cont: li_defaults){
        std::vector<std::pair<int, cplx_type> > & sbus_dp = li_sbus_dp[cont_id];
        const std::vector<int> branch_ids(this_cont.branches.begin(), this_cont.branches.end());
        std::vector<Coeff> coeffs = _coeffs_for_branch_ids(branch_ids, grid_model, ac_solver_used, id_me_to_solver, n_line,
                                                           sbus_dp);
        _append_shunt_trips(std::vector<int>(this_cont.shunts.begin(), this_cont.shunts.end()),
                            grid_model, ac_solver_used, id_me_to_solver, coeffs, sbus_dp);
        _append_load_trips(std::vector<int>(this_cont.loads.begin(), this_cont.loads.end()),
                           grid_model, ac_solver_used, id_me_to_solver, sbus_dp);
        _append_gen_trips(std::vector<int>(this_cont.gens.begin(), this_cont.gens.end()),
                          grid_model, ac_solver_used, id_me_to_solver, sbus_dp);
        li_coeffs.push_back(std::move(coeffs));
        ++cont_id;
    }
}

std::vector<int> YbusPolicy::Contingency::shunt_ids_for_row(Eigen::Index row) const
{
    std::vector<int> shunt_ids;
    if(shunt_mask.rows() == 0 || row >= shunt_mask.rows()) return shunt_ids;
    for(Eigen::Index col = 0; col < shunt_mask.cols(); ++col){
        if(shunt_mask(row, col)) shunt_ids.push_back(static_cast<int>(col));
    }
    return shunt_ids;
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
    if(static_cast<size_t>(row) < topo_branches_off.size() && !topo_branches_off[static_cast<size_t>(row)].empty()){
        // both lists are sorted and, by construction, disjoint (see
        // BaseBatchSweep::_maybe_resolve_topology)
        const std::vector<int> & topo = topo_branches_off[static_cast<size_t>(row)];
        std::vector<int> merged;
        merged.reserve(branch_ids.size() + topo.size());
        std::set_union(branch_ids.begin(), branch_ids.end(), topo.begin(), topo.end(), std::back_inserter(merged));
        branch_ids.swap(merged);
    }
    return branch_ids;
}

void YbusPolicy::Contingency::_append_placement_coeffs(const BranchPlacement & placement,
                                                        const LSGrid & grid_model,
                                                        bool ac_solver_used,
                                                        size_t n_line,
                                                        std::vector<Coeff> & out)
{
    const bool is_trafo = static_cast<size_t>(placement.branch_id) >= n_line;
    const int el_id = is_trafo ? placement.branch_id - static_cast<int>(n_line) : placement.branch_id;
    const BranchContainer & branch = is_trafo ? static_cast<const BranchContainer &>(grid_model.get_trafos_as_data())
                                              : static_cast<const BranchContainer &>(grid_model.get_powerlines_as_data());
    cplx_type y_ff, y_ft, y_tf, y_tt;
    if(ac_solver_used){
        y_ff = branch.yac_11()[el_id];
        y_ft = branch.yac_12()[el_id];
        y_tf = branch.yac_21()[el_id];
        y_tt = branch.yac_22()[el_id];
    }else{
        y_ff = branch.ydc_11()[el_id];
        y_ft = branch.ydc_12()[el_id];
        y_tf = branch.ydc_21()[el_id];
        y_tt = branch.ydc_22()[el_id];
    }
    const int b1 = placement.bus1_solver;
    const int b2 = placement.bus2_solver;
    out.push_back({b1, b1, -y_ff});
    out.push_back({b2, b2, -y_tt});
    out.push_back({b1, b2, -y_ft});
    out.push_back({b2, b1, -y_tf});
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
    li_sbus_dp.assign(static_cast<size_t>(nb_steps), {});
    li_rows.clear();
    for(Eigen::Index row = 0; row < nb_steps; ++row){
        std::vector<std::pair<int, cplx_type> > & row_shift_dp = li_sbus_dp[static_cast<size_t>(row)];
        std::vector<Coeff> coeffs = _coeffs_for_branch_ids(branch_ids_for_row(row, n_line), grid_model, ac_solver_used, id_me_to_solver, n_line,
                                                           row_shift_dp);
        // its shunts (the loads and generators it trips are SbusPolicy::Vary's)
        _append_shunt_trips(shunt_ids_for_row(row), grid_model, ac_solver_used, id_me_to_solver, coeffs, row_shift_dp);
        if(static_cast<size_t>(row) < topo_branches_moved.size()){
            for(const BranchPlacement & placement : topo_branches_moved[static_cast<size_t>(row)]){
                // the base contribution out (where there is one), the row's own in. A row
                // that moves a branch is AC only (a DC row with an action is refused), so
                // there is no phase-shift injection to carry over.
                if(placement.base_on){
                    std::vector<std::pair<int, cplx_type> > unused_shift_dp;
                    const std::vector<Coeff> base = _coeffs_for_branch_ids(
                        std::vector<int>(1, placement.branch_id), grid_model, ac_solver_used, id_me_to_solver, n_line,
                        unused_shift_dp);
                    coeffs.insert(coeffs.end(), base.begin(), base.end());
                }
                _append_placement_coeffs(placement, grid_model, ac_solver_used, n_line, coeffs);
            }
        }
        li_coeffs.push_back(coeffs);
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
