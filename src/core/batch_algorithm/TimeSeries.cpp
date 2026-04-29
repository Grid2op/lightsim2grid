// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "TimeSeries.hpp"

#include <iostream>
#include <sstream>

namespace ls2g {

int TimeSeries::compute_Vs(Eigen::Ref<const RealMat> gen_p,
                           Eigen::Ref<const RealMat> sgen_p,
                           Eigen::Ref<const RealMat> load_p,
                           Eigen::Ref<const RealMat> load_q,
                           const CplxVect & Vinit,
                           const int max_iter,
                           const real_type tol)
{
    auto timer = CustTimer();
    auto timer_preproc = CustTimer();

    // perform some initial checks and reset timers
    size_t nb_total_bus = _reset_data_and_check_vinit(Vinit);
    _status = 0;

    // read from the grid the usefull information
    const auto & sn_mva = _grid_model.get_sn_mva();
    const bool ac_solver_used = _algo.ac_solver_used();
    size_t nb_steps = gen_p.rows();

    // prepare the gridmodel (compute Ybus, Sbus etc.)
    CplxVect Vinit_solver = prepare_solver_input_base(Vinit, ac_solver_used);   

    ///////////////////////////////////////////
    // initialize what will change (here Sbus)
    // retrieve necessary data
    const auto & generators = _grid_model.get_generators_as_data();
    const auto & s_generators = _grid_model.get_static_generators_as_data();
    const auto & loads = _grid_model.get_loads_as_data();

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
    //////////////////////////////////////////

    bool init_powerflow_has_conv = _finish_preprocessing(
        nb_steps,
        nb_total_bus,
        Vinit_solver,  // is modified if _init_from_n_powerflow is true
        max_iter,
        tol,
        timer_preproc
    );

    if(!init_powerflow_has_conv){
        _status = 0;
        return -1;
    }

    // compute the powerflows
    // set the "right" init vector (either the one provided by the user or the one
    // after the initial powerflow)
    CplxVect V = Vinit_solver;

    int step_diverge = -1;
    const real_type tol_ = tol / sn_mva; 
    bool conv;
    if(!ac_solver_used) _algo_controler.tell_recompute_sbus(); // we need to recompute Sbus (DC case)
    for(size_t i = 0; i < nb_steps; ++i){
        conv = false;
        conv = compute_one_powerflow(Ybus_,
                                     V, 
                                     _Sbuses.row(i),
                                     slack_ids_me_.as_eigen(),
                                     slack_weights_,
                                     bus_pv_.as_eigen(),
                                     bus_pq_.as_eigen(),
                                     max_iter,
                                     tol_);
        // nothing changes
        _algo_controler.tell_none_changed(); 
        if(!ac_solver_used) _algo_controler.tell_recompute_sbus(); // we need to recompute Sbus (DC case)
        if(!conv){
            _timer_total = timer.duration();
            step_diverge = i;
            _status = 0;
            return _status;
        }
        if(conv && step_diverge < 0) {
            _voltages.row(i)(id_solver_to_me_.as_eigen()) = V.array();
        }
    }

    // If i reached there, it means it is succesfull
    _status = 1;
    _timer_total = timer.duration();
    return _status;
}

} // namespace ls2g
