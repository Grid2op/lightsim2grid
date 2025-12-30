// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ContingencyAnalysis.hpp"

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

void ContingencyAnalysis::init_li_coeffs(
    bool ac_solver_used,
    const std::vector<SolverBusId> &id_me_to_solver)
{
    _li_coeffs.clear();
    _li_coeffs.reserve(_li_defaults.size());
    const auto & powerlines = _grid_model.get_powerlines_as_data();
    const auto & trafos = _grid_model.get_trafos_as_data();
    // const auto & id_me_to_solver = ac_solver_used ? _grid_model.id_me_to_ac_solver(): _grid_model.id_me_to_dc_solver();
    int bus_1_id, bus_2_id;
    cplx_type y_ff, y_ft, y_tf, y_tt;
    bool status;
    for(const auto & this_cont_id: _li_defaults){
        std::vector<Coeff> this_cont_coeffs;
        this_cont_coeffs.reserve(this_cont_id.size() * 4);  // usually there are 4 coeffs per powerlines / trafos
        for(auto br_id : this_cont_id){
            int el_id;
            const TwoSidesContainer_rxh_A<OneSideContainer_ForBranch> *p_branch;
            if(br_id < n_line_)
            {
                // this is a powerline
                el_id = br_id;
                p_branch = & powerlines;
            }else{
                // this is a trafo
                el_id = br_id - n_line_;
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
            // TODO disconnected one side !
            if(ac_solver_used){
                y_ff = p_branch->yac_11()[el_id];
                y_ft = p_branch->yac_12()[el_id];
                y_tf = p_branch->yac_21()[el_id];
                y_tt = p_branch->yac_22()[el_id];
            }else{
                y_ff = p_branch->ydc_11()[el_id];
                y_ft = p_branch->ydc_12()[el_id];
                y_tf = p_branch->ydc_21()[el_id];
                y_tt = p_branch->ydc_22()[el_id];
            }

            if(status)
            {
                // element is connected, update coeffs based on status of each powerlines
                if((bus_1_id != GenericContainer::_deactivated_bus_id)) this_cont_coeffs.push_back({bus_1_id, bus_1_id, y_ff});
                if((bus_2_id != GenericContainer::_deactivated_bus_id)) this_cont_coeffs.push_back({bus_2_id, bus_2_id, y_tt});
                if((bus_1_id != GenericContainer::_deactivated_bus_id) && (bus_2_id != GenericContainer::_deactivated_bus_id)){
                    this_cont_coeffs.push_back({bus_1_id, bus_2_id, y_ft});
                    this_cont_coeffs.push_back({bus_2_id, bus_1_id, y_tf});
                }
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

IntVect ContingencyAnalysis::is_grid_connected_after_contingency(){
    const bool ac_solver_used = _solver.ac_solver_used();
    Eigen::SparseMatrix<cplx_type> Ybus = ac_solver_used ? _grid_model.get_Ybus_solver() : _grid_model.get_dcYbus_solver();
    IntVect res = IntVect::Constant(_li_coeffs.size(), 0);
    int cont_id = 0;
    for(const auto & coeffs_modif: _li_coeffs){
        if(remove_from_Ybus(Ybus, coeffs_modif)) res(cont_id) = 1;
        else res(cont_id) = 0;
        readd_to_Ybus(Ybus, coeffs_modif);
        ++cont_id;
    }
    return res;
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
    // std::cout << "entering compute\n";
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

    // redo a powerflow to prepare the solver
    // if(ac_solver_used){
    //     if(_grid_model.get_solver().get_type() != _solver.get_type())
    //     {
    //         _grid_model.change_solver(_solver.get_type());
    //         _grid_model.ac_pf(Vinit, max_iter, tol);
    //         _grid_model.unset_changes();
    //     }
    // }else{
    //     if(_grid_model.get_dc_solver().get_type() != _solver.get_type())
    //     {
    //         _grid_model.change_solver(_solver.get_type());
    //         _grid_model.dc_pf(Vinit, max_iter, tol);
    //         _grid_model.unset_changes();
    //     }
    // }

    // prepare the gridmodel (compute Ybus, Sbus etc.)
    CplxVect Vinit_solver = prepare_solver_input_base(Vinit, ac_solver_used);

    // _solver_control.tell_all_changed();
    // CplxVect Sbus;
    // Eigen::SparseMatrix<cplx_type> Ybus;
    // std::vector<GlobalBusId> id_solver_to_me;
    // std::vector<SolverBusId> id_me_to_solver;
    // SolverBusIdVect slack_ids_solver;
    // GlobalBusIdVect slack_ids_me;
    // CplxVect Vinit_solver = _grid_model.pre_process_solver(
    //     Vinit, 
    //     Sbus,
    //     Ybus,
    //     id_me_to_solver,
    //     id_solver_to_me,
    //     slack_ids_me,
    //     slack_ids_solver,
    //     ac_solver_used,
    //     _solver_control);

    // std::cout << "end pre_process_solver\n";

    // Eigen::SparseMatrix<cplx_type> Ybus = ac_solver_used ? _grid_model.get_Ybus_solver() : _grid_model.get_dcYbus_solver();
    // const Eigen::Index nb_buses_solver = Ybus_.cols();
    // // const std::vector<GlobalBusId> & id_solver_to_me = ac_solver_used ? _grid_model.id_ac_solver_to_me() : _grid_model.id_dc_solver_to_me();
    // const SolverBusIdVect & gm_bus_pv = _grid_model.get_pv_solver();
    // // std::cout << "end get_pv_solver\n";
    // const SolverBusIdVect & gm_bus_pq = _grid_model.get_pq_solver();
    // // std::cout << "end get_pq_solver\n";
    // const IntVect & bus_pv = _to_intvect(gm_bus_pv);
    // // std::cout << "end bus_pv\n";
    // const IntVect & bus_pq = _to_intvect(gm_bus_pq);
    // // std::cout << "end gm_bus_pq\n";
    // const IntVect & slack_ids = _to_intvect(slack_ids_solver_);
    // // const IntVect & slack_ids = 
    // //     ac_solver_used ? 
    // //     _to_intvect(static_cast<const SolverBusIdVect &>(_grid_model.get_slack_ids_solver())): 
    // //     _to_intvect(static_cast<const SolverBusIdVect &>(_grid_model.get_slack_ids_dc_solver()));
    // const RealVect & slack_weights = _grid_model.get_slack_weights_solver();
    // const auto & id_me_to_solver = ac_solver_used ? _grid_model.id_me_to_ac_solver() :  _grid_model.id_me_to_dc_solver();
    
    // get the proper Sbus vector
    // CplxVect Sbus = CplxVect::Zero(nb_buses_solver);
    // _grid_model.fillSbus_other(Sbus, ac_solver_used, id_me_to_solver); 

    // initialize properly the coefficients that I will need to remove
    init_li_coeffs(ac_solver_used, id_me_to_solver_);
    Eigen::Index nb_steps = _li_defaults.size();

    // init the results matrices
    _voltages = BaseBatchSolverSynch::CplxMat::Zero(nb_steps, nb_total_bus); 
    _amps_flows = RealMat::Zero(0, n_total_);

    // reset the solver
    _solver.reset();

    // compute the right Vinit to send to the solver
    // CplxVect Vinit_solver = extract_Vsolver_from_Vinit(Vinit, nb_buses_solver, nb_total_bus, id_me_to_solver);

    // perform the initial powerflow / "powerflow in n"
    _solver_control.tell_all_changed();
    _solver.tell_solver_control(_solver_control);
    // std::cout << "Ybus.size() " << Ybus.size() << std::endl;
    // std::cout << "Vinit_solver.size() " << Vinit_solver.size() << std::endl;
    // std::cout << "Sbus.size() " << Sbus.size() << std::endl;
    // std::cout << "slack_ids.size() " << slack_ids.size() << std::endl;
    // std::cout << "slack_weights.size() " << slack_weights.size() << std::endl;
    // std::cout << "bus_pv.size() " << bus_pv.size() << std::endl;
    // std::cout << "bus_pq.size() " << bus_pq.size() << std::endl;
    // _grid_model.get_generators().set_vm(Vinit_solver, id_me_to_solver_);
    CplxVect Vinit_solver2 = Vinit_solver;
    bool conv = _solver.compute_pf(
        Ybus_,
        Vinit_solver2,
        Sbus_,
        slack_ids_,
        slack_weights_,
        bus_pv_,
        bus_pq_,
        max_iter,
        tol);

    // end of pre processing
    _timer_pre_proc = timer_preproc.duration();
    if(!conv) return;
    _solver_control.tell_none_changed();

    // now perform the security analysis
    Eigen::Index cont_id = 0;
    CplxVect V;
    // std::cout << "entering SA " << std::endl;
    for(const auto & coeffs_modif: _li_coeffs){
        auto timer_modif_Ybus = CustTimer();
        bool invertible = true;
        // no need to add to this Ybus as DC solver have an internal Ybus which is updated with _solver.update_internal_Ybus
        if (ac_solver_used) invertible = remove_from_Ybus(Ybus_, coeffs_modif);
        _timer_modif_Ybus += timer_modif_Ybus.duration();
        conv = false;

        if(invertible)
        {
            if(!ac_solver_used)
            {
                // DC solver stores the ybus internally, I update it
                // instead of building it over and over
                for(const Coeff& coeff : coeffs_modif){
                    _solver.update_internal_Ybus(coeff, false);  // false => remove the coeff (using -= )
                }
            }
            V = Vinit_solver; // Vinit is reused for each contingencies
            // _solver_control.tell_all_changed();
            // _solver_control.tell_solver_need_reset();
            // std::cout << "\t compute_one_powerflow\n";
            conv = compute_one_powerflow(
                Ybus_,
                V,
                Sbus_,
                slack_ids_,
                slack_weights_,
                bus_pv_,
                bus_pq_,
                max_iter,
                tol / sn_mva);
            if(!ac_solver_used)
            {
                // DC solver stores the ybus internally, I update it
                // instead of building it over and over
                for(const Coeff& coeff : coeffs_modif){
                    _solver.update_internal_Ybus(coeff, true);  // true => add back the coeff (using += )
                }
            }
        }
        // std::string conv_str =  conv ? "has converged" : "has diverged";
        // std::cout << "contingency " << cont_id << ": " << conv_str << std::endl;
        // if(!conv) std::cout << "\t error was: " << _solver.get_error() << std::endl;

        timer_modif_Ybus = CustTimer();
        // std::cout << "\t end compute_one_powerflow\n";
        // no need to add to this Ybus as DC solver have an internal Ybus which is updated with _solver.update_internal_Ybus
        if (ac_solver_used) readd_to_Ybus(Ybus_, coeffs_modif); 
        _timer_modif_Ybus += timer_modif_Ybus.duration();
        if (conv && invertible) _voltages.row(cont_id)(reinterpret_cast<const std::vector<int> & >(id_solver_to_me_)) = V.array();
        ++cont_id;
    }
    _timer_total = timer.duration();
    // std::cout << "end SA\n";
}

// by default the flows are not 0 when the powerline is connected in the original topology
// this function sorts this out
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
