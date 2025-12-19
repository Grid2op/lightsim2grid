// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "TrafoContainer.hpp"

#include <iostream>
#include <sstream>

void TrafoContainer::init(const RealVect & trafo_r,
                          const RealVect & trafo_x,
                          const CplxVect & trafo_b,
                          const RealVect & trafo_tap_step_pct,
                          const RealVect & trafo_tap_pos,
                          const RealVect & trafo_shift_degree,
                          const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                          const Eigen::VectorXi & trafo_hv_id,
                          const Eigen::VectorXi & trafo_lv_id
                          )
{
    /**
    INPUT DATA ARE ALREADY PAIR UNIT !!
    DOES NOT WORK WITH POWERLINES
    **/
    const int size = static_cast<int>(trafo_r.size());
    GenericContainer::check_size(trafo_tap_step_pct, size, "trafo_tap_step_pct");
    GenericContainer::check_size(trafo_tap_pos, size, "trafo_tap_pos");
    //TODO "parrallel" in the pandapower dataframe, like for lines, are not handled. Handle it python side!

    RealVect ratio = my_one_ + 0.01 * trafo_tap_step_pct.array() * trafo_tap_pos.array();


    init(trafo_r, trafo_x, trafo_b, ratio, trafo_shift_degree, trafo_tap_hv, trafo_hv_id, trafo_lv_id);
}

void TrafoContainer::init(const RealVect & trafo_r,
                          const RealVect & trafo_x,
                          const CplxVect & trafo_b,
                          const RealVect & trafo_ratio,
                          const RealVect & trafo_shift_degree,
                          const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                          const Eigen::VectorXi & trafo_hv_id,
                          const Eigen::VectorXi & trafo_lv_id
                          )
{
    const int size = static_cast<int>(trafo_r.size());
    GenericContainer::check_size(trafo_r, size, "trafo_r");
    GenericContainer::check_size(trafo_x, size, "trafo_x");
    GenericContainer::check_size(trafo_b, size, "trafo_b");
    GenericContainer::check_size(trafo_ratio, size, "trafo_ratio");
    GenericContainer::check_size(trafo_shift_degree, size, "trafo_shift_degree");
    GenericContainer::check_size(trafo_tap_hv, static_cast<std::vector<bool>::size_type>(size), "trafo_tap_hv");
    GenericContainer::check_size(trafo_hv_id, size, "trafo_hv_id");
    GenericContainer::check_size(trafo_lv_id, size, "trafo_lv_id");

    r_ = trafo_r;
    x_ = trafo_x;
    h_side_1_ = 0.5 * trafo_b;
    h_side_2_ = 0.5 * trafo_b;
    ratio_ = trafo_ratio;
    shift_ = trafo_shift_degree / my_180_pi_;  // do not forget conversion degree / rad here !
    is_tap_hv_side_ = trafo_tap_hv;
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
     std::vector<bool> is_tap_hv_side = is_tap_hv_side_;
     TrafoContainer::StateRes res(
        get_tsc_rxha_state(),
        ratio,
        is_tap_hv_side,
        shift);
     return res;
}

void TrafoContainer::set_state(TrafoContainer::StateRes & my_state)
{
    set_tsc_rxha_state(std::get<0>(my_state));

    std::vector<real_type> & ratio = std::get<1>(my_state);
    std::vector<bool> & is_tap_hv_side = std::get<2>(my_state);
    std::vector<real_type> & shift = std::get<3>(my_state);

    auto size = nb();
    GenericContainer::check_size(ratio, size, "ratio");
    GenericContainer::check_size(is_tap_hv_side, size, "is_tap_hv_side");
    GenericContainer::check_size(shift, size, "shift");

    // input data
    ratio_  = RealVect::Map(&ratio[0], size);
    shift_  = RealVect::Map(&shift[0], size);
    is_tap_hv_side_ = is_tap_hv_side;
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
    // const cplx_type h = h_side_1_(i);
    real_type tau = ratio_(el_id);
    real_type theta_shift = shift_(el_id);
    if(!is_tap_hv_side_[el_id]){
        tau = my_one_ / tau;
        theta_shift = -theta_shift;
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

    yac_11_(el_id) = (ys + h_side_1_(el_id)) / (tau * tau);
    yac_22_(el_id) = (ys + h_side_2_(el_id));
    yac_21_(el_id) = -ys / tau * emitheta_shift ;
    yac_12_(el_id) = -ys / tau * eitheta_shift;

    // for DC
    // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.21
    // except here I only care about the real part, so I remove the "1/j"
    cplx_type tmp = 1. / (tau * x_(el_id));
    ydc_11_(el_id) = tmp;
    ydc_22_(el_id) = tmp;
    ydc_21_(el_id) = -tmp;
    ydc_12_(el_id) = -tmp;
    if(!is_tap_hv_side_[el_id]) dc_x_tau_shift_(el_id) = -std::real(tmp) * theta_shift;
    else dc_x_tau_shift_(el_id) = std::real(tmp) * theta_shift;
}

void TrafoContainer::hack_Sbus_for_dc_phase_shifter(CplxVect & Sbus,
                                                    bool ac,
                                                    const std::vector<SolverBusId> & id_grid_to_solver)
{
    if(ac) return;
    // return;
    const int nb_trafo = nb();
    int bus_id_me, bus_id_solver_hv, bus_id_solver_lv;
    // cplx_type tmp;
    for(int trafo_id = 0; trafo_id < nb_trafo; ++trafo_id){
        //  i don't do anything if the load is disconnected
        if(!status_global_[trafo_id]) continue;
        if(abs(dc_x_tau_shift_[trafo_id]) < _tol_equal_float) continue; // nothing to do if the trafo is not concerned (no phase shifter)
        bus_id_me = get_bus_side_2(trafo_id);
        bus_id_solver_lv = id_grid_to_solver[bus_id_me];
        if(bus_id_solver_lv == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (lv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        bus_id_me = get_bus_side_1(trafo_id);
        bus_id_solver_hv = id_grid_to_solver[bus_id_me];
        if(bus_id_solver_hv == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (hv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        Sbus.coeffRef(bus_id_solver_hv) += dc_x_tau_shift_[trafo_id];
        Sbus.coeffRef(bus_id_solver_lv) -= dc_x_tau_shift_[trafo_id];
    }
}

void TrafoContainer::fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                                std::vector<Eigen::Triplet<real_type> > & Bpp,
                                const std::vector<SolverBusId> & id_grid_to_solver,
                                real_type sn_mva,
                                FDPFMethod xb_or_bx) const
{

    // For Bp
    // temp_branch[:, BR_B] = zeros(nl)           ## zero out line charging shunts
    // temp_branch[:, TAP] = ones(nl)             ## cancel out taps
    // if alg == 2:                               ## if XB method
    //    temp_branch[:, BR_R] = zeros(nl)       ## zero out line resistance

    // For Bpp
    // temp_branch[:, SHIFT] = zeros(nl)          ## zero out phase shifters
    // if alg == 3:                               ## if BX method
    //     temp_branch[:, BR_R] = zeros(nl)    ## zero out line resistance
    const Eigen::Index nb_trafo = nb();
    real_type yft_bp, ytf_bp, yff_bp, ytt_bp;
    real_type yft_bpp, ytf_bpp, yff_bpp, ytt_bpp;

    for(Eigen::Index tr_id=0; tr_id < nb_trafo; ++tr_id){
        // i only add this if the powerline is connected
        if(!status_global_[tr_id]) continue;

        // get the from / to bus id
        int bus_or_id_me = get_bus_side_1(tr_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::fillBp_Bpp: the trafo with id ";
            exc_ << tr_id;
            exc_ << " is connected (hv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_ex_id_me = get_bus_side_2(tr_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::fillBp_Bpp: the trafo with id ";
            exc_ << tr_id;
            exc_ << " is connected (lv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }

        // get the coefficients
        // tau is needed for Bpp
        double tau_bpp = ratio_(tr_id);
        // for Bp we need shift
        real_type theta_shift = shift_(tr_id);
        if(!is_tap_hv_side_[tr_id]){
            tau_bpp = 1. / ratio_(tr_id);
            theta_shift = -shift_(tr_id); 
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
        yff_bp = ys_bp_r;
        ytt_bp = ys_bp_r;
        ytf_bp = -std::imag(ys_bp * emitheta_shift_bp);
        yft_bp = -std::imag(ys_bp * eitheta_shift_bp);
        const real_type ys_bpp_r = std::imag(ys_bpp); 
        yff_bpp = (ys_bpp_r + std::imag(h_side_1_(tr_id))) / (tau_bpp * tau_bpp);
        ytt_bpp = ys_bpp_r + std::imag(h_side_2_(tr_id));
        ytf_bpp = -ys_bpp_r / tau_bpp;
        yft_bpp = -ys_bpp_r / tau_bpp;

        // and now add them
        Bp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id, bus_ex_solver_id, -yft_bp));
        Bp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id, bus_or_solver_id, -ytf_bp));
        Bp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id, bus_or_solver_id, -yff_bp));
        Bp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id, bus_ex_solver_id, -ytt_bp));

        Bpp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id, bus_ex_solver_id, -yft_bpp));
        Bpp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id, bus_or_solver_id, -ytf_bpp));
        Bpp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id, bus_or_solver_id, -yff_bpp));
        Bpp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id, bus_ex_solver_id, -ytt_bpp));

    }
}