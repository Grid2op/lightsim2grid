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

int TimeSeries::compute_Vs(const Eigen::Ref<const RealMat> & gen_p,
                           const Eigen::Ref<const RealMat> & sgen_p,
                           const Eigen::Ref<const RealMat> & load_p,
                           const Eigen::Ref<const RealMat> & load_q,
                           const Eigen::Ref<const CplxVect> & Vinit,
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
    const auto check_mat = [nb_steps](const Eigen::Ref<const RealMat> & mat, Eigen::Index nb_expected_cols,
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

    if(abs(sn_mva - 1.0) > _tol_equal_float) _Sbuses.array() /= static_cast<cplx_type>(sn_mva);

    // ... and then everything else the gridmodel stamps into Sbus, which the four
    // matrices above do NOT cover. They carry generator / static-generator ACTIVE
    // power and load active+reactive power, and nothing more -- while
    // LSGrid::fillSbus_me also stamps storage units, REACTIVE_POWER-mode SVCs, the
    // reactive setpoint of a generator whose voltage regulation is off, a static
    // generator's reactive setpoint, the hvdc injections, and (dc only) the phase
    // shifter hack. None of those is a per-step input, so all of them are constant
    // across the steps and must come from the grid.
    //
    // This is derived rather than enumerated on purpose: an enumeration falling
    // behind is exactly the bug being fixed here (only the hvdc term was listed, so
    // storages, SVCs and every reactive setpoint were silently dropped -- worth up
    // to ~0.08 pu of voltage error on a 4-bus test). Instead take the COMPLETE
    // injection the gridmodel just built and subtract the share the per-step
    // matrices are responsible for, evaluated at the gridmodel's OWN injection
    // values. What remains is by construction everything else, so an element type
    // added to fillSbus_me later is picked up here for free.
    _Sbuses.rowwise() += _constant_sbus_pu(ac_solver_used).transpose();
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

CplxVect TimeSeries::_constant_sbus_pu(bool ac_solver_used) const
{
    // the complete per-unit injection, straight from the gridmodel
    // (`prepare_solver_input_base` filled exactly one of these two)
    CplxVect res;
    if(ac_solver_used){
        res = Sbus_;
    } else {
        // dc keeps everything real; the imaginary part is unused downstream
        // (compute_one_powerflow takes Sbus.real()) but must be zero so that
        // subtracting the reconstruction below cannot invent one
        res = Pbus_.cast<cplx_type>();
    }
    if(static_cast<int>(res.size()) != nb_buses_solver_){
        std::ostringstream exc_;
        exc_ << "TimeSeries::_constant_sbus_pu: the gridmodel injection has "
             << res.size() << " entries while the solver has " << nb_buses_solver_
             << " buses. prepare_solver_input_base() must run first.";
        throw std::runtime_error(exc_.str());
    }

    // ... minus the share compute_Vs rebuilds per step, evaluated at the gridmodel's
    // own target values, so that the two cancel and only the rest is left
    const auto & generators = _grid_model.get_generators_as_data();
    const auto & s_generators = _grid_model.get_static_generators_as_data();
    const auto & loads = _grid_model.get_loads_as_data();

    // The gridmodel's target vectors, presented to fill_SBus_* as the single-step
    // matrices it expects, WITHOUT copying them: a 1 x n row-major matrix has exactly
    // the memory layout of a contiguous n-vector, so an Eigen::Map over the vector's
    // own storage is a view. Handing `.transpose()` over instead would not be -- the
    // transpose of a column vector is not a RealMat, so binding it to the
    // Eigen::Ref<const RealMat> parameter has to materialise it. The Refs are held in
    // named locals rather than being called inline: that way the Map never depends on
    // a temporary Ref's lifetime, whatever the accessors return.
    const Eigen::Ref<const RealVect> gen_p_v = _grid_model.get_gen_target_p();
    const Eigen::Ref<const RealVect> sgen_p_v = _grid_model.get_sgen_target_p();
    const Eigen::Ref<const RealVect> load_p_v = _grid_model.get_load_target_p();
    const Eigen::Ref<const RealVect> load_q_v = loads.get_target_q();
    const Eigen::Map<const RealMat> gen_p(gen_p_v.data(), 1, gen_p_v.size());
    const Eigen::Map<const RealMat> sgen_p(sgen_p_v.data(), 1, sgen_p_v.size());
    const Eigen::Map<const RealMat> load_p(load_p_v.data(), 1, load_p_v.size());
    const Eigen::Map<const RealMat> load_q(load_q_v.data(), 1, load_q_v.size());

    // `accounted` does have to be its own buffer: it is written by fill_SBus_* and
    // then scaled by sn_mva, which must NOT touch `res` (already per-unit). One
    // vector-sized allocation per compute_Vs call, not per step.
    CplxMat accounted = CplxMat::Zero(1, nb_buses_solver_);
    bool add_ = true;
    fill_SBus_real(accounted, generators, gen_p, id_me_to_solver_, add_);
    fill_SBus_real(accounted, s_generators, sgen_p, id_me_to_solver_, add_);
    add_ = false;
    fill_SBus_real(accounted, loads, load_p, id_me_to_solver_, add_);
    fill_SBus_imag(accounted, loads, load_q, id_me_to_solver_, add_);

    const real_type sn_mva = _grid_model.get_sn_mva();
    if(abs(sn_mva - 1.0) > _tol_equal_float) accounted.array() /= static_cast<cplx_type>(sn_mva);
    res -= accounted.row(0).transpose();
    return res;
}

} // namespace ls2g
