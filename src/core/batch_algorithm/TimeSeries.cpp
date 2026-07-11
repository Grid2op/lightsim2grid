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

    // Validate the shapes of the caller-supplied matrices before fill_SBus_* uses them.
    // fill_SBus_* iterates el_id up to the grid's element count and reads
    // temporal_data.col(el_id) (unchecked Eigen .col()), so a matrix with too few
    // columns over-reads; and _Sbuses is sized with gen_p.rows(), so mismatched row
    // counts turn the `Sbuses.col(...) += tmp` into an out-of-bounds read/write.
    const auto check_mat = [nb_steps](Eigen::Ref<const RealMat> mat, Eigen::Index nb_expected_cols,
                                      const std::string & name){
        if(static_cast<size_t>(mat.rows()) != nb_steps){
            std::ostringstream exc_;
            exc_ << "TimeSeries::compute_Vs: '" << name << "' has " << mat.rows()
                 << " rows (time steps) while 'gen_p' has " << nb_steps
                 << ". All injection matrices must share the same number of rows.";
            throw std::runtime_error(exc_.str());
        }
        if(mat.cols() != nb_expected_cols){
            std::ostringstream exc_;
            exc_ << "TimeSeries::compute_Vs: '" << name << "' has " << mat.cols()
                 << " columns while the grid counts " << nb_expected_cols
                 << " such elements. The number of columns must match the number of elements.";
            throw std::runtime_error(exc_.str());
        }
    };
    check_mat(gen_p, generators.nb(), "gen_p");
    check_mat(sgen_p, s_generators.nb(), "sgen_p");
    check_mat(load_p, loads.nb(), "load_p");
    check_mat(load_q, loads.nb(), "load_q");

    // now build the Sbus
    _Sbuses = CplxMat::Zero(nb_steps, nb_buses_solver_);

    bool add_ = true;
    fill_SBus_real(_Sbuses, generators, gen_p, id_me_to_solver_, add_);
    fill_SBus_real(_Sbuses, s_generators, sgen_p, id_me_to_solver_, add_);
    add_ = false;
    fill_SBus_real(_Sbuses, loads, load_p, id_me_to_solver_, add_);
    fill_SBus_imag(_Sbuses, loads, load_q, id_me_to_solver_, add_);

    // add the (constant accross the steps) hvdc injections: fixed-setpoint
    // lines, station reactive setpoints / LCC consumptions and, in dc, the
    // saturated droop lines. The angle-droop flows of the linear-mode lines
    // are theta dependent: they are handled by the solver itself (Hvdc
    // extension of the NR system in ac, dc matrix term in dc).
    CplxVect hvdc_sbus = CplxVect::Zero(nb_buses_solver_);
    _grid_model.get_dclines_as_data().fillSbus(hvdc_sbus, id_me_to_solver_, ac_solver_used);
    _Sbuses.rowwise() += hvdc_sbus.transpose();

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
                                     slack_ids_solver_.as_eigen(),
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
