// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "TimeSeries.h"
#include <iostream>
#include <sstream>

int TimeSeries::compute_Vs(Eigen::Ref<const RealMat> gen_p,
                           Eigen::Ref<const RealMat> sgen_p,
                           Eigen::Ref<const RealMat> load_p,
                           Eigen::Ref<const RealMat> load_q,
                           const CplxVect & Vinit,
                           const int max_iter,
                           const real_type tol)
{
    auto timer = CustTimer();
    const Eigen::Index nb_total_bus = _grid_model.total_bus();
    if(Vinit.size() != nb_total_bus){
        std::ostringstream exc_;
        exc_ << "TimeSeries::compute_Sbuses: Size of the Vinit should be the same as the total number of buses. Currently:  ";
        exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_total_bus << " buses.";
        exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
        throw std::runtime_error(exc_.str());
    }

    // init everything
    _status = 0;
    _nb_solved = 0;
    _timer_pre_proc = 0.;
    _timer_total = 0.;
    _timer_solver = 0.;

    auto timer_preproc = CustTimer();

    // init the computations
    const auto & sn_mva = _grid_model.get_sn_mva();
    const bool ac_solver_used = _solver.ac_solver_used();
    Eigen::SparseMatrix<cplx_type> Ybus = ac_solver_used ? _grid_model.get_Ybus() : _grid_model.get_dcYbus();
    const Eigen::Index nb_buses_solver = Ybus.cols();

    const auto & id_me_to_solver = ac_solver_used ? _grid_model.id_me_to_ac_solver() :  _grid_model.id_me_to_dc_solver();
    const auto & id_solver_to_me = ac_solver_used ? _grid_model.id_ac_solver_to_me() : _grid_model.id_dc_solver_to_me();
    const auto & generators = _grid_model.get_generators_as_data();
    const auto & s_generators = _grid_model.get_static_generators_as_data();
    const auto & loads = _grid_model.get_loads_as_data();

    const Eigen::Index nb_steps = gen_p.rows();

    const Eigen::VectorXi & bus_pv = _grid_model.get_pv();
    const Eigen::VectorXi & bus_pq = _grid_model.get_pq();
    const Eigen::VectorXi & slack_ids  = ac_solver_used ? _grid_model.get_slack_ids(): _grid_model.get_slack_ids_dc();
    const RealVect & slack_weights = _grid_model.get_slack_weights();
    _solver.reset();
    _solver_control.tell_none_changed();
    _solver_control.tell_recompute_sbus();

    // now build the Sbus
    _Sbuses = CplxMat::Zero(nb_steps, nb_buses_solver);

    bool add_ = true;
    fill_SBus_real(_Sbuses, generators, gen_p, id_me_to_solver, add_);
    fill_SBus_real(_Sbuses, s_generators, sgen_p, id_me_to_solver, add_);
    add_ = false;
    fill_SBus_real(_Sbuses, loads, load_p, id_me_to_solver, add_);
    fill_SBus_imag(_Sbuses, loads, load_q, id_me_to_solver, add_);
    if(sn_mva != 1.0) _Sbuses.array() /= static_cast<cplx_type>(sn_mva);
    // TODO trafo hack for Sbus !

    // init the results matrices
    _voltages = BaseBatchSolverSynch::CplxMat::Zero(nb_steps, nb_total_bus); 
    _amps_flows = RealMat::Zero(0, n_total_);

    // extract V solver from the given V
    CplxVect Vinit_solver = extract_Vsolver_from_Vinit(Vinit, nb_buses_solver, nb_total_bus, id_me_to_solver);
    _timer_pre_proc = timer_preproc.duration();

    // compute the powerflows
    // set the "right" init vector
    CplxVect V = Vinit_solver;
    _grid_model.get_generators().set_vm(V, id_me_to_solver);

    bool conv;
    Eigen::Index step_diverge = -1;
    const real_type tol_ = tol / sn_mva; 
    // do the computation for each step
    for(Eigen::Index i = 0; i < nb_steps; ++i){
        conv = false;
        conv = compute_one_powerflow(Ybus,
                                     V, 
                                     _Sbuses.row(i),
                                     slack_ids,
                                     slack_weights,
                                     bus_pv,
                                     bus_pq,
                                     max_iter,
                                     tol_);
        if(!conv){
            _timer_total = timer.duration();
            return _status;
        }
        if(conv && step_diverge < 0) _voltages.row(i)(id_solver_to_me) = V.array();
        else step_diverge = i;
    }
    if(step_diverge > 0){
        _status = 0;
    }else{
        // If i reached there, it means it is succesfull
        _status = 1;
    }
    
    _timer_total = timer.duration();
    return _status;
}
