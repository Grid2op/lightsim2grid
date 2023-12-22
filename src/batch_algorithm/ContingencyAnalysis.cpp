// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ContingencyAnalysis.h"

#include <queue>
#include <math.h>       /* isfinite */

bool ContingencyAnalysis::check_invertible(const Eigen::SparseMatrix<cplx_type> & Ybus) const{
    std::vector<bool> visited(Ybus.cols(), false); 
    std::vector<bool> already_added(Ybus.cols(), false);
    std::queue<Eigen::Index> neighborhood;
    Eigen::Index col_id = 0;  // start by node 0, why not
    while (true)
    {
        visited[col_id] = true;
        for (Eigen::SparseMatrix<cplx_type>::InnerIterator it(Ybus, col_id); it; ++it)
        {
            // add in the queue all my neighbor (if the coefficient is big enough)
            if(!visited[it.row()] && !already_added[it.row()] && abs(it.value()) > 1e-8){
                neighborhood.push(it.row());
                already_added[it.row()] = true;
            }
        }
        if(neighborhood.empty()) break;
        col_id = neighborhood.front();
        neighborhood.pop();
    }
    
    bool ok = true;
    for(auto el: visited){
        if(!el)
        {
            // this node has not been visited, there is an error
            ok=false;
            break;
        }
    }
    return ok;
}

void ContingencyAnalysis::init_li_coeffs(bool ac_solver_used){
    _li_coeffs.clear();
    _li_coeffs.reserve(_li_defaults.size());
    const auto & powerlines = _grid_model.get_powerlines_as_data();
    const auto & trafos = _grid_model.get_trafos_as_data();
    const auto & id_me_to_solver = ac_solver_used ? _grid_model.id_me_to_ac_solver(): _grid_model.id_me_to_dc_solver();
    Eigen::Index bus_1_id, bus_2_id;
    cplx_type y_ff, y_ft, y_tf, y_tt;
    for(const auto & this_cont_id: _li_defaults){
        std::vector<Coeff> this_cont_coeffs;
        this_cont_coeffs.reserve(this_cont_id.size() * 4);  // usually there are 4 coeffs per powerlines / trafos
        for(auto line_id : this_cont_id){
            if(line_id < n_line_)
            {
                // this is a powerline
                bus_1_id = id_me_to_solver[powerlines.get_bus_from()[line_id]];
                bus_2_id = id_me_to_solver[powerlines.get_bus_to()[line_id]];
                if(ac_solver_used){
                    y_ff = powerlines.yac_ff()[line_id];
                    y_ft = powerlines.yac_ft()[line_id];
                    y_tf = powerlines.yac_tf()[line_id];
                    y_tt = powerlines.yac_tt()[line_id];
                }else{
                    y_ff = powerlines.ydc_ff()[line_id];
                    y_ft = powerlines.ydc_ft()[line_id];
                    y_tf = powerlines.ydc_tf()[line_id];
                    y_tt = powerlines.ydc_tt()[line_id];
                }
            }else{
                // this is a trafo
                const auto trafo_id = line_id - n_line_;
                bus_1_id = id_me_to_solver[trafos.get_bus_from()[trafo_id]];
                bus_2_id = id_me_to_solver[trafos.get_bus_to()[trafo_id]];
                if(ac_solver_used){
                    y_ff = trafos.yac_ff()[trafo_id];
                    y_ft = trafos.yac_ft()[trafo_id];
                    y_tf = trafos.yac_tf()[trafo_id];
                    y_tt = trafos.yac_tt()[trafo_id];
                }else{
                    y_ff = trafos.ydc_ff()[trafo_id];
                    y_ft = trafos.ydc_ft()[trafo_id];
                    y_tf = trafos.ydc_tf()[trafo_id];
                    y_tt = trafos.ydc_tt()[trafo_id];
                }
            }

            if(bus_1_id != GenericContainer::_deactivated_bus_id && bus_2_id != GenericContainer::_deactivated_bus_id)
            {
                // element is connected
                this_cont_coeffs.push_back({bus_1_id, bus_1_id, y_ff});
                this_cont_coeffs.push_back({bus_1_id, bus_2_id, y_ft});
                this_cont_coeffs.push_back({bus_2_id, bus_1_id, y_tf});
                this_cont_coeffs.push_back({bus_2_id, bus_2_id, y_tt});
            }
        }
        _li_coeffs.push_back(this_cont_coeffs);
    }
}


bool ContingencyAnalysis::remove_from_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                                        const std::vector<Coeff> & coeffs) const
{
    for(const auto & coeff_to_remove: coeffs){
        Ybus.coeffRef(coeff_to_remove.row_id, coeff_to_remove.col_id) -= coeff_to_remove.value;
    }
    return check_invertible(Ybus);
}
void ContingencyAnalysis::readd_to_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                                     const std::vector<Coeff> & coeffs) const
{
    for(const auto & coeff_to_remove: coeffs){
        Ybus.coeffRef(coeff_to_remove.row_id, coeff_to_remove.col_id) += coeff_to_remove.value;
    }
}

void ContingencyAnalysis::compute(const CplxVect & Vinit, int max_iter, real_type tol)
{
    auto timer = CustTimer();
    auto timer_preproc = CustTimer();

    _timer_modif_Ybus = 0.;
    _timer_pre_proc = 0.;
    _timer_total = 0.;
    _timer_solver = 0.;

    const Eigen::Index nb_total_bus = _grid_model.total_bus();
    if(Vinit.size() != nb_total_bus){
        std::ostringstream exc_;
        exc_ << "SecurityAnalysis::compute: Size of the Vinit should be the same as the total number of buses. Currently:  ";
        exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_total_bus << " buses.";
        exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
        throw std::runtime_error(exc_.str());
    }

    // read from the grid the usefull information
    const auto & sn_mva = _grid_model.get_sn_mva();
    const bool ac_solver_used = _solver.ac_solver_used();
    Eigen::SparseMatrix<cplx_type> Ybus = ac_solver_used ? _grid_model.get_Ybus() : _grid_model.get_dcYbus();
    const Eigen::Index nb_buses_solver = Ybus.cols();
    const auto & id_solver_to_me = ac_solver_used ? _grid_model.id_ac_solver_to_me() : _grid_model.id_dc_solver_to_me();
    const Eigen::VectorXi & bus_pv = _grid_model.get_pv();
    const Eigen::VectorXi & bus_pq = _grid_model.get_pq();
    const Eigen::VectorXi & slack_ids = ac_solver_used ? _grid_model.get_slack_ids(): _grid_model.get_slack_ids_dc();
    const RealVect & slack_weights = _grid_model.get_slack_weights();
    const auto & id_me_to_solver = ac_solver_used ? _grid_model.id_me_to_ac_solver() :  _grid_model.id_me_to_dc_solver();
    
    // get the proper Sbus vector
    CplxVect Sbus = CplxVect::Zero(nb_buses_solver);
    _grid_model.fillSbus_other(Sbus, ac_solver_used, id_me_to_solver); 

    // initialize properly the coefficients that I will need to remove
    init_li_coeffs(ac_solver_used);
    Eigen::Index nb_steps = _li_defaults.size();

    // init the results matrices
    _voltages = BaseBatchSolverSynch::CplxMat::Zero(nb_steps, nb_total_bus); 
    _amps_flows = RealMat::Zero(0, n_total_);

    // reset the solver
    _solver.reset();
    _solver_control.tell_none_changed();
    _solver_control.tell_recompute_ybus();
    // _solver_control.tell_ybus_some_coeffs_zero();
    // ybus does not change sparsity pattern here

    // compute the right Vinit to send to the solver
    CplxVect Vinit_solver = extract_Vsolver_from_Vinit(Vinit, nb_buses_solver, nb_total_bus, id_me_to_solver);

    // end of pre processing
    _timer_pre_proc = timer_preproc.duration();

    // now perform the security analysis
    Eigen::Index cont_id = 0;
    bool conv;
    CplxVect V;
    // int contingency = 0;
    for(const auto & coeffs_modif: _li_coeffs){
        auto timer_modif_Ybus = CustTimer();
        bool invertible = remove_from_Ybus(Ybus, coeffs_modif);
        _timer_modif_Ybus += timer_modif_Ybus.duration();
        conv = false;

        // I have absolutely no idea why, but if i add this "if"
        // which reduces the computation time, it somehow increase it by A LOT !
        // 5.2ms without it vs 81.9ms with it (for the iee 118)
        // So better make the computation, even if it's not used...

        if(invertible)
        {
            V = Vinit_solver; // Vinit is reused for each contingencies
            conv = compute_one_powerflow(Ybus, V, Sbus,
                                         slack_ids, slack_weights,
                                         bus_pv, bus_pq,
                                         max_iter,
                                         tol / sn_mva);
        }
        // std::string conv_str =  conv ? "has converged" : "has diverged";
        // std::cout << "contingency " << contingency << ": " << conv_str << std::endl;
        // if(!conv) std::cout << "\t error was: " << _solver.get_error() << std::endl;
        // ++contingency;

        timer_modif_Ybus = CustTimer();
        readd_to_Ybus(Ybus, coeffs_modif);
        _timer_modif_Ybus += timer_modif_Ybus.duration();
        if (conv && invertible) _voltages.row(cont_id)(id_solver_to_me) = V.array();
        ++cont_id;
    }
    _timer_total = timer.duration();
}

void ContingencyAnalysis::clean_flows(bool is_amps)
{
    auto timer = CustTimer();
    Eigen::Index cont_id = 0;
    for(const auto & l_id_this_cont: _li_defaults){
        for(auto l_id : l_id_this_cont){
            real_type & el = is_amps ? _amps_flows(cont_id, l_id): _active_power_flows(cont_id, l_id);
            if(isfinite(el)) el = 0.;
        }
        ++cont_id;
    }
    if (is_amps) _timer_compute_A += timer.duration();
    else _timer_compute_P += timer.duration();
}
