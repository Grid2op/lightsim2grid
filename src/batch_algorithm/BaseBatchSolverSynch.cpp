// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BaseBatchSolverSynch.h"

/**
 V is modified at each call !
**/
bool BaseBatchSolverSynch::compute_one_powerflow(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                                  CplxVect & V,
                                                  const CplxVect & Sbus,
                                                  const Eigen::VectorXi & slack_ids,
                                                  const RealVect & slack_weights,
                                                  const Eigen::VectorXi & bus_pv,
                                                  const Eigen::VectorXi & bus_pq,
                                                  int max_iter,
                                                  double tol
                                                  )
{
    _solver.tell_solver_control(_solver_control);
    bool conv = _solver.compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, bus_pv, bus_pq, max_iter, tol);
    if(conv){
        V = _solver.get_V().array();
    }
    ++_nb_solved;
    _timer_solver += _solver.get_computation_time(); 
    return conv;
}

void BaseBatchSolverSynch::compute_flows_from_Vs(bool amps)
{
    // TODO find a way to factorize that with TrafoContainer::compute_results
    // TODO and LineContainer::compute_results

    if (_voltages.size() == 0)
    {
        std::ostringstream exc_;
        exc_ << "BaseMultiplePowerflow::compute_flows_from_Vs: cannot compute the flows as the voltages are not set. Have you called compute(...) ? ";
        throw std::runtime_error(exc_.str());
    }
    if (amps) _timer_compute_A = 0.;
    else _timer_compute_P = 0.;

    auto timer_compute = CustTimer();
    const auto & sn_mva = _grid_model.get_sn_mva();
    const auto & nb_steps = _voltages.rows();

    // reset the results
    if (amps) _amps_flows = RealMat::Zero(nb_steps, n_total_);
    else _active_power_flows = RealMat::Zero(nb_steps, n_total_);
    
    // compute the flows for the powerlines
    Eigen::Index lag_id = 0;
    if (amps) compute_amps_flows(_grid_model.get_powerlines_as_data(), sn_mva, lag_id, false);
    else compute_active_power_flows(_grid_model.get_powerlines_as_data(), sn_mva, lag_id, false);

    // compute the flows for the trafos
    lag_id = n_line_;
    if (amps) compute_amps_flows(_grid_model.get_trafos_as_data(), sn_mva, lag_id, true);
    else compute_active_power_flows(_grid_model.get_trafos_as_data(), sn_mva, lag_id, true);

    if (amps) _timer_compute_A = timer_compute.duration();
    else _timer_compute_P = timer_compute.duration();
}
