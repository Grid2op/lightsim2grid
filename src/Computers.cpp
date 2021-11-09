// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "Computers.h"
#include <iostream>
#include <sstream>

int Computers::compute_Vs(Eigen::Ref<const RealMat> gen_p,
                          Eigen::Ref<const RealMat> sgen_p,
                          Eigen::Ref<const RealMat> load_p,
                          Eigen::Ref<const RealMat> load_q,
                          const CplxVect & Vinit,
                          int max_iter,
                          real_type tol)
{
    auto timer = CustTimer();
    const Eigen::Index nb_total_bus = _grid_model.total_bus();
    if(Vinit.size() != nb_total_bus){
        std::ostringstream exc_;
        exc_ << "Computers::compute_Sbuses: Size of the Vinit should be the same as the total number of buses. Currently:  ";
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
    const auto & sn_mva = _grid_model.get_sn_mva();
    const auto & Ybus = _grid_model.get_Ybus(); 
    const auto & id_me_to_ac_solver = _grid_model.id_me_to_ac_solver();
    const auto & id_ac_solver_to_me = _grid_model.id_ac_solver_to_me();
    const auto & generators = _grid_model.get_generators_as_data();
    const auto & s_generators = _grid_model.get_static_generators_as_data();
    const auto & loads = _grid_model.get_loads_as_data();

    const Eigen::Index nb_steps = gen_p.rows();
    const Eigen::Index nb_buses_solver = Ybus.cols();  // which is equal to Ybus.rows();

    const Eigen::VectorXi & bus_pv = _grid_model.get_pv();
    const Eigen::VectorXi & bus_pq = _grid_model.get_pq();
    _solver.reset();

    // init the computations
    // now build the Sbus
    _Sbuses = CplxMat::Zero(nb_steps, nb_buses_solver);

    bool add_ = true;
    fill_SBus_real(_Sbuses, generators, gen_p, id_me_to_ac_solver, add_);
    fill_SBus_real(_Sbuses, s_generators, sgen_p, id_me_to_ac_solver, add_);
    add_ = false;
    fill_SBus_real(_Sbuses, loads, load_p, id_me_to_ac_solver, add_);
    fill_SBus_imag(_Sbuses, loads, load_q, id_me_to_ac_solver, add_);
    if(sn_mva != 1.0) _Sbuses.array() /= static_cast<cplx_type>(sn_mva);

    // init the results matrices
    _voltages = BaseMultiplePowerflow::CplxMat::Zero(nb_steps, nb_total_bus); 
    _amps_flows = RealMat::Zero(0, n_total_);

    // extract V solver from the given V
    CplxVect Vinit_solver = extract_Vsolver_from_Vinit(Vinit, nb_buses_solver, nb_total_bus, id_me_to_ac_solver);
    _timer_pre_proc = timer_preproc.duration();

    // compute the powerflows
    // do the computa   tion for each step
    CplxVect V = Vinit_solver;
    bool conv;
    CplxVect Sbus;
    for(Eigen::Index i = 0; i < nb_steps; ++i){
        Sbus = _Sbuses.row(i).array(); 
        conv = compute_one_powerflow(Ybus, V, Sbus,
                                     bus_pv, bus_pq,
                                     max_iter,
                                     tol / sn_mva);
        if(!conv){
            _timer_total = timer.duration();
            return _status;
        }
        _voltages.row(i)(id_ac_solver_to_me) = V.array();
    }

    // If i reached there, it means it is succesfull
    _status = 1;
    _timer_total = timer.duration();
    return _status;
}
