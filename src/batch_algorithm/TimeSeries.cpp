// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "TimeSeries.hpp"

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
    CplxVect Vinit_solver = prepare_solver_input_base(Vinit, ac_solver_used);   

    // retrieve necessary data
    const auto & generators = _grid_model.get_generators_as_data();
    const auto & s_generators = _grid_model.get_static_generators_as_data();
    const auto & loads = _grid_model.get_loads_as_data();
    const Eigen::Index nb_steps = gen_p.rows();

    // now build the Sbus
    _Sbuses = CplxMat::Zero(nb_steps, nb_buses_solver_);

    bool add_ = true;
    fill_SBus_real(_Sbuses, generators, gen_p, id_me_to_solver_, add_);
    fill_SBus_real(_Sbuses, s_generators, sgen_p, id_me_to_solver_, add_);
    add_ = false;
    fill_SBus_real(_Sbuses, loads, load_p, id_me_to_solver_, add_);
    fill_SBus_imag(_Sbuses, loads, load_q, id_me_to_solver_, add_);
    if(abs(sn_mva - 1.0) > _tol_equal_float) _Sbuses.array() /= static_cast<cplx_type>(sn_mva);
    // TODO trafo hack for Sbus !

    // init the results matrices
    _voltages = BaseBatchSolverSynch::CplxMat::Zero(nb_steps, nb_total_bus); 
    _amps_flows = RealMat::Zero(0, n_total_);

    // end of the pre processing steps
    _timer_pre_proc = timer_preproc.duration();

    // compute the powerflows
    // set the "right" init vector
    CplxVect V = Vinit_solver;
    _grid_model.get_generators().set_vm(V, id_me_to_solver_);

    Eigen::Index step_diverge = -1;
    const real_type tol_ = tol / sn_mva; 
    bool conv;
    // do the computation for each step
    _solver_control.tell_all_changed();  // recompute everything at the first iteration
    for(Eigen::Index i = 0; i < nb_steps; ++i){
        conv = false;
        conv = compute_one_powerflow(Ybus_,
                                     V, 
                                     _Sbuses.row(i),
                                     slack_ids_,
                                     slack_weights_,
                                     bus_pv_,
                                     bus_pq_,
                                     max_iter,
                                     tol_);
        // nothing changes
        _solver_control.tell_none_changed(); 
        if(!ac_solver_used) _solver_control.tell_recompute_sbus(); // we need to recompute Sbus (DC case)
        if(!conv){
            _timer_total = timer.duration();
            step_diverge = i;
            _status = 0;
            return _status;
        }
        if(conv && step_diverge < 0) _voltages.row(i)(reinterpret_cast<const std::vector<int> & >(id_solver_to_me_)) = V.array();
    }

    // If i reached there, it means it is succesfull
    _status = 1;
    _timer_total = timer.duration();
    return _status;
}
