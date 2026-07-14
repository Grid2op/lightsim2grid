// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BaseBatchSolverSynch.hpp"

namespace ls2g {

/**
 V is modified at each call !
 Member version: single-threaded path, uses the member solver / control / accumulators.
**/
bool BaseBatchSolverSynch::compute_one_powerflow(
    const Eigen::SparseMatrix<cplx_type> & Ybus,
    CplxVect & V,
    const CplxVect & Sbus,
    Eigen::Ref<const IntVect> slack_ids,
    const RealVect & slack_weights,
    Eigen::Ref<const IntVect> bus_pv,
    Eigen::Ref<const IntVect> bus_pq,
    int max_iter,
    double tol
)
{
    return compute_one_powerflow(
        _algo, _algo_controler, _nb_solved, _timer_solver,
        Ybus, V, Sbus, slack_ids, slack_weights, bus_pv, bus_pq, max_iter, tol);
}

/**
 V is modified at each call !
 Explicit version: operates on the passed solver / control / accumulators so it
 can run concurrently (one solver per thread). The member Bbus_ is read-only here.
**/
bool BaseBatchSolverSynch::compute_one_powerflow(
    AlgorithmSelector & algo,
    AlgoControl & control,
    int & nb_solved,
    double & timer_solver,
    const Eigen::SparseMatrix<cplx_type> & Ybus,
    CplxVect & V,
    const CplxVect & Sbus,
    Eigen::Ref<const IntVect> slack_ids,
    const RealVect & slack_weights,
    Eigen::Ref<const IntVect> bus_pv,
    Eigen::Ref<const IntVect> bus_pq,
    int max_iter,
    double tol
)
{
    algo.tell_solver_control(control);
    bool conv;
    if(algo.ac_solver_used()){
        conv = algo.compute_pf(
            Ybus,
            V,
            Sbus,
            slack_ids,
            slack_weights,
            bus_pv,
            bus_pq,
            max_iter,
            tol);
    } else {
        // native real DC entry point: Bbus_ is the (constant) real admittance built in
        // prepare_solver_input_base; Pbus is the real part of the (possibly per-step)
        // injection. ContingencyAnalysis leaves Sbus_ empty in DC (it relies on the
        // member Pbus_ built once), so fall back to Pbus_ when no complex Sbus is given.
        // Split into two branches (rather than a ternary) so the common/no-Sbus case
        // binds Pbus_ directly instead of unifying both branches into a fresh copy.
        if(Sbus.size() == 0){
            conv = algo.compute_pf_dc(Bbus_, V, Pbus_, slack_ids, slack_weights, bus_pv, bus_pq);
        } else {
            const RealVect Pbus = Sbus.real();
            conv = algo.compute_pf_dc(Bbus_, V, Pbus, slack_ids, slack_weights, bus_pv, bus_pq);
        }
    }
    if(conv){
        V = algo.get_V().array();
    }
    ++nb_solved;
    timer_solver += algo.get_computation_time();
    return conv;
}

bool BaseBatchSolverSynch::warmup_solver(
    AlgorithmSelector & algo,
    AlgoControl & control,
    CplxVect Vinit_solver,
    int max_iter,
    real_type tol)
{
    algo.reset();
    control.tell_all_changed();
    algo.tell_solver_control(control);
    bool conv;
    if(algo.ac_solver_used()){
        conv = algo.compute_pf(
            Ybus_,
            Vinit_solver,
            Sbus_,
            slack_ids_solver_.as_eigen(),
            slack_weights_,
            bus_pv_.as_eigen(),
            bus_pq_.as_eigen(),
            max_iter,
            tol);
    } else {
        conv = algo.compute_pf_dc(
            Bbus_,
            Vinit_solver,
            Pbus_,
            slack_ids_solver_.as_eigen(),
            slack_weights_,
            bus_pv_.as_eigen(),
            bus_pq_.as_eigen());
    }
    // subsequent per-contingency solves reuse the factorization
    control.tell_none_changed();
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
    size_t lag_id = 0;
    if (amps) compute_amps_flows(_grid_model.get_powerlines_as_data(), sn_mva, lag_id, false);
    else compute_active_power_flows(_grid_model.get_powerlines_as_data(), sn_mva, lag_id, false);

    // compute the flows for the trafos
    lag_id = n_line_;
    if (amps) compute_amps_flows(_grid_model.get_trafos_as_data(), sn_mva, lag_id, true);
    else compute_active_power_flows(_grid_model.get_trafos_as_data(), sn_mva, lag_id, true);

    if (amps) _timer_compute_A = timer_compute.duration();
    else _timer_compute_P = timer_compute.duration();
}

} // namespace ls2g
