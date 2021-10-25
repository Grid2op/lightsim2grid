// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "Computers.h"

int Computers::compute_Sbuses(Eigen::Ref<CplxMat> Sbuses,
                              const CplxVect & Vinit,
                              int max_iter,
                              real_type tol)
{
    auto timer = CustTimer();
    const Eigen::Index nb_bus = _grid_model.total_bus();
    if(Vinit.size() != nb_bus){
        std::ostringstream exc_;
        exc_ << "Computers::compute_Sbuses: Size of the Vinit should be the same as the total number of buses. Currently:  ";
        exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_bus << " buses.";
        exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
        throw std::runtime_error(exc_.str());
    }

    std::cout << "Computers::compute_Sbuses -2" << std::endl;
    const Eigen::Index nb_steps = Sbuses.rows();
    const Eigen::Index nb_buses = Sbuses.cols();

    // init everything
    _status = 0;
    _nb_solved = 0;
    const auto & sn_mva = _grid_model.get_sn_mva();
    const auto & Ybus = _grid_model.get_Ybus(); 
    const auto & id_me_to_ac_solver = _grid_model.id_me_to_ac_solver();

    const Eigen::VectorXi bus_pv = _grid_model.get_pv();
    const Eigen::VectorXi bus_pq = _grid_model.get_pq();
    _solver.reset();
    std::cout << "Computers::compute_Sbuses -1" << std::endl;

    // init the computations
    _voltages = Computers::CplxMat(nb_steps, nb_buses);

    // TODO extract V solver from the given V
    // TODO extract S solver from the given Sbus
    std::cout << "Computers::compute_Sbuses 0" << std::endl;
    CplxVect Sbus = Sbuses.row(0); 
    CplxVect V = Vinit;
    CplxVect Vtmp = Vinit;
    // TODO check the size of Sbuses (nb cols) and size of Ybus (they should match !)
    bool conv = _solver.compute_pf(Ybus, Vtmp, Sbus, bus_pv, bus_pq, max_iter, tol / sn_mva);
    if(!conv){
        // divergence !
        _timer_total = timer.duration();
        return _status;
    }
    V.array() = 1.0 * _solver.get_V().array();
    _voltages.row(0) = 1.0 * V.array();
    ++_nb_solved;
    std::cout << "Computers::compute_Sbuses 2" << std::endl;
    // do the computa   tion for each step
    for(Eigen::Index i = 1; i < nb_steps; ++i){
        Sbus = Sbuses.row(i).array(); 
        CplxVect Vtmp2 = Vinit;
        conv = _solver.compute_pf(Ybus, Vtmp2, Sbus, bus_pv, bus_pq, max_iter, tol / sn_mva);
        if(!conv){
            // divergence !
            _timer_total = timer.duration();
            return _status;
        }
        V.array() = 1.0 * _solver.get_V().array();
        _voltages.row(i) = 1.0 * V.array();
        _timer_solver += _solver.get_computation_time(); 
        ++_nb_solved;
    }

    // If i reached there, it means it is succesfull
    _status = 1;
    _timer_total = timer.duration();
    return _status;
}
