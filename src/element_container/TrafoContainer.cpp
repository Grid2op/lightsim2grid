// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "TrafoContainer.hpp"

#include <iostream>
#include <sstream>

void TrafoContainer::init(
    const RealVect & trafo_r,
    const RealVect & trafo_x,
    const CplxVect & trafo_b,
    const RealVect & trafo_tap_step_pct,
    const RealVect & trafo_tap_pos,
    const RealVect & trafo_shift_degree,
    const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
    const Eigen::VectorXi & trafo_hv_id,
    const Eigen::VectorXi & trafo_lv_id,
    bool ignore_tap_side_for_shift
) {
    /**
    INPUT DATA ARE ALREADY PAIR UNIT !!
    DOES NOT WORK WITH POWERLINES
    **/
    const int size = static_cast<int>(trafo_r.size());
    GenericContainer::check_size(trafo_tap_step_pct, size, "trafo_tap_step_pct");
    GenericContainer::check_size(trafo_tap_pos, size, "trafo_tap_pos");
    //TODO "parrallel" in the pandapower dataframe, like for lines, are not handled. Handle it python side!

    RealVect ratio = my_one_ + 0.01 * trafo_tap_step_pct.array() * trafo_tap_pos.array();

    init(trafo_r, trafo_x, trafo_b, ratio, trafo_shift_degree, 
        trafo_tap_hv, trafo_hv_id, trafo_lv_id, ignore_tap_side_for_shift);
}

void TrafoContainer::init(
    const RealVect & trafo_r,
    const RealVect & trafo_x,
    const CplxVect & trafo_b,
    const RealVect & trafo_ratio,
    const RealVect & trafo_shift_degree,
    const std::vector<bool> & trafo_tap_side1,  // is tap on high voltage (true) or low voltate
    const Eigen::VectorXi & trafo_hv_id,
    const Eigen::VectorXi & trafo_lv_id,
    bool ignore_tap_side_for_shift
) {
    const int size = static_cast<int>(trafo_r.size());
    GenericContainer::check_size(trafo_r, size, "trafo_r");
    GenericContainer::check_size(trafo_x, size, "trafo_x");
    GenericContainer::check_size(trafo_b, size, "trafo_b");
    GenericContainer::check_size(trafo_ratio, size, "trafo_ratio");
    GenericContainer::check_size(trafo_shift_degree, size, "trafo_shift_degree");
    GenericContainer::check_size(trafo_tap_side1, static_cast<std::vector<bool>::size_type>(size), "trafo_tap_hv");
    GenericContainer::check_size(trafo_hv_id, size, "trafo_hv_id");
    GenericContainer::check_size(trafo_lv_id, size, "trafo_lv_id");

    r_ = trafo_r;
    x_ = trafo_x;
    h_side_1_ = 0.5 * trafo_b;
    h_side_2_ = 0.5 * trafo_b;
    ratio_ = trafo_ratio;
    shift_ = trafo_shift_degree / my_180_pi_;  // do not forget conversion degree / rad here !
    is_tap_side1_ = trafo_tap_side1;
    ignore_tap_side_for_shift_ = ignore_tap_side_for_shift;
    // bus_hv_id_ = trafo_hv_id;
    // bus_lv_id_ = trafo_lv_id;
    // status_ = std::vector<bool>(trafo_r.size(), true);
    init_tsc(trafo_hv_id, trafo_lv_id, "trafo");
    _update_model_coeffs();
    reset_results();
}

TrafoContainer::StateRes TrafoContainer::get_state() const
{
     std::vector<real_type> ratio(ratio_.begin(), ratio_.end());
     std::vector<real_type> shift(shift_.begin(), shift_.end());
     std::vector<bool> is_tap_hv_side = is_tap_side1_;
     TrafoContainer::StateRes res(
        get_tsc_rxha_state(),
        ratio,
        is_tap_hv_side,
        shift,
        ignore_tap_side_for_shift_);
     return res;
}

void TrafoContainer::set_state(TrafoContainer::StateRes & my_state)
{
    set_tsc_rxha_state(std::get<0>(my_state));

    std::vector<real_type> & ratio = std::get<1>(my_state);
    std::vector<bool> & is_tap_side1 = std::get<2>(my_state);
    std::vector<real_type> & shift = std::get<3>(my_state);

    auto size = nb();
    GenericContainer::check_size(ratio, size, "ratio");
    GenericContainer::check_size(is_tap_side1, size, "is_tap_side1");
    GenericContainer::check_size(shift, size, "shift");

    ratio_  = RealVect::Map(&ratio[0], size);
    shift_  = RealVect::Map(&shift[0], size);
    is_tap_side1_ = is_tap_side1;
    ignore_tap_side_for_shift_ = std::get<4>(my_state);
    _update_model_coeffs();
    reset_results();
}

void TrafoContainer::_update_model_coeffs()
{
    const Eigen::Index my_size = r_.size();

    yac_11_ = CplxVect::Zero(my_size);
    yac_12_ = CplxVect::Zero(my_size);
    yac_21_ = CplxVect::Zero(my_size);
    yac_22_ = CplxVect::Zero(my_size);

    ydc_11_ = CplxVect::Zero(my_size);
    ydc_12_ = CplxVect::Zero(my_size);
    ydc_21_ = CplxVect::Zero(my_size);
    ydc_22_ = CplxVect::Zero(my_size);
    dc_x_tau_shift_ = RealVect::Zero(my_size);
    for(Eigen::Index i = 0; i < my_size; ++i)
    {
        _update_model_coeffs_one_el(i);
    }
}

void TrafoContainer::_update_model_coeffs_one_el(int el_id)
{
    // for AC
    // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.2
    const cplx_type ys = 1. / cplx_type(r_(el_id), x_(el_id));
    real_type tau = ratio_(el_id);
    real_type theta_shift = shift_(el_id);
    if(!is_tap_side1_[el_id]){
        tau = my_one_ / tau;

        // pnadapower uses tap_side only for ratio, not for
        // phase shift apparently
        if (!ignore_tap_side_for_shift_) theta_shift = -theta_shift;
    }
    cplx_type eitheta_shift  = {my_one_, my_zero_};  // exp(j  * alpha)
    cplx_type emitheta_shift = {my_one_, my_zero_};  // exp(-j * alpha)
    if(std::abs(theta_shift) > _tol_equal_float)
    {
        real_type cos_theta = std::cos(theta_shift);
        real_type sin_theta = std::sin(theta_shift);
        eitheta_shift = {cos_theta, sin_theta};
        emitheta_shift = {cos_theta, -sin_theta};
    }
    real_type _1_tau = my_one_ / tau; // 1 / tau
    yac_11_(el_id) = (ys + h_side_1_(el_id)) * _1_tau * _1_tau;  // (ys + h1) / tau**2
    yac_12_(el_id) = -ys * _1_tau * eitheta_shift;  // -ys / (tau * exp(-j.theta_shift))
    
    yac_21_(el_id) = -ys * _1_tau * emitheta_shift;  // -ys / (tau * exp(j.theta_shift))
    yac_22_(el_id) = (ys + h_side_2_(el_id));  // ys + h2

    // for DC
    // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.21
    // except here I only care about the real part, so I remove the "1/j"
    cplx_type tmp = 1. / x_(el_id) * _1_tau;
    ydc_11_(el_id) = tmp;
    ydc_22_(el_id) = tmp;
    ydc_21_(el_id) = -tmp;
    ydc_12_(el_id) = -tmp;

    dc_x_tau_shift_(el_id) = -std::real(tmp) * theta_shift;
    // if(ignore_tap_side_for_shift_)
    // {
    //     dc_x_tau_shift_(el_id) = -std::real(tmp) * theta_shift;
    // } else {
    //     if(is_tap_side1_[el_id]) dc_x_tau_shift_(el_id) = -std::real(tmp) * theta_shift;
    //     else dc_x_tau_shift_(el_id) = std::real(tmp) * theta_shift;
    // }
}

void TrafoContainer::hack_Sbus_for_dc_phase_shifter(
    CplxVect & Sbus,
    bool ac,
    const std::vector<SolverBusId> & id_grid_to_solver)
{
    if(ac) return;

    // return;
    const int nb_trafo = nb();
    const std::vector<bool> & status1 = side_1_.get_status();
    const std::vector<bool> & status2 = side_2_.get_status();
    GlobalBusId bus_id_me;
    SolverBusId bus_id_solver_hv, bus_id_solver_lv;
    // cplx_type tmp;
    for(int trafo_id = 0; trafo_id < nb_trafo; ++trafo_id){
        //  i don't do anything if the load is disconnected
        if(!status_global_[trafo_id]) continue;
        if(!status1[trafo_id]) continue;
        if(!status2[trafo_id]) continue;

        if(abs(dc_x_tau_shift_[trafo_id]) < _tol_equal_float) continue; // nothing to do if the trafo is not concerned (no phase shifter)
        
        bus_id_me = get_bus_side_2(trafo_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: (GridModelId) the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (side 2) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver_lv = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver_lv.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: (SolverId) the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (side 2) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }

        bus_id_me = get_bus_side_1(trafo_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: (GridModelId) the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (side 1) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver_hv = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver_hv.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: (SolverId) the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (side 1) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        Sbus.coeffRef(bus_id_solver_hv.cast_int()) -= dc_x_tau_shift_[trafo_id];
        Sbus.coeffRef(bus_id_solver_lv.cast_int()) += dc_x_tau_shift_[trafo_id];
    }
}

TrafoContainer::FDPFCoeffs TrafoContainer::get_fdpf_coeffs(int tr_id, FDPFMethod xb_or_bx) const{
    TrafoContainer::FDPFCoeffs res;
    // get the coefficients
    // tau is needed for Bpp
    double tau_bpp = ratio_(tr_id);
    // for Bp we need shift
    real_type theta_shift = shift_(tr_id);
    if(!is_tap_side1_[tr_id]){
        tau_bpp = 1. / ratio_(tr_id);
        
        // pnadapower uses tap_side only for ratio, not for
        // phase shift apparently
        if (!ignore_tap_side_for_shift_) theta_shift = -shift_(tr_id); 
    }
    cplx_type eitheta_shift_bp  = {my_one_, my_zero_};  // exp(j  * alpha)
    cplx_type emitheta_shift_bp = {my_one_, my_zero_};  // exp(-j * alpha)
    if(abs(theta_shift) > _tol_equal_float)
    {
        const real_type cos_theta = std::cos(theta_shift);
        const real_type sin_theta = std::sin(theta_shift);
        eitheta_shift_bp = {cos_theta, sin_theta};
        emitheta_shift_bp = {cos_theta, -sin_theta};
    }

    // depending on XB or BX we define the y differently
    cplx_type ys_bp, ys_bpp;
    if(xb_or_bx==FDPFMethod::XB){
        ys_bp = 1. / (0. + my_i * x_(tr_id));
        ys_bpp = 1. / (r_(tr_id) + my_i * x_(tr_id));
    }else if (xb_or_bx==FDPFMethod::BX){
        ys_bp = 1. / (r_(tr_id) + my_i * x_(tr_id));
        ys_bpp = 1. / (0. + my_i * x_(tr_id));
    }else{
        std::ostringstream exc_;
        exc_ << "TrafoContainer::fillBp_Bpp: unknown method for the FDPF powerflow for trafo id ";
        exc_ << tr_id;
        throw std::runtime_error(exc_.str());            
    }

    const real_type ys_bp_r = std::imag(ys_bp); 
    res.yff_bp = ys_bp_r;
    res.ytt_bp = ys_bp_r;
    res.ytf_bp = -std::imag(ys_bp * emitheta_shift_bp);
    res.yft_bp = -std::imag(ys_bp * eitheta_shift_bp);
    const real_type ys_bpp_r = std::imag(ys_bpp); 
    res.yff_bpp = (ys_bpp_r + std::imag(h_side_1_(tr_id))) / (tau_bpp * tau_bpp);
    res.ytt_bpp = ys_bpp_r + std::imag(h_side_2_(tr_id));
    res.ytf_bpp = -ys_bpp_r / tau_bpp;
    res.yft_bpp = -ys_bpp_r / tau_bpp;
    return res;
}
