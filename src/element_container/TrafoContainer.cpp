// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "TrafoContainer.h"

#include <iostream>
#include <sstream>

void TrafoContainer::init(const RealVect & trafo_r,
                          const RealVect & trafo_x,
                          const CplxVect & trafo_b,
                          const RealVect & trafo_tap_step_pct,
            //                       const RealVect & trafo_tap_step_degree,
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
    GenericContainer::check_size(trafo_r, size, "trafo_r");
    GenericContainer::check_size(trafo_x, size, "trafo_x");
    GenericContainer::check_size(trafo_b, size, "trafo_b");
    GenericContainer::check_size(trafo_tap_step_pct, size, "trafo_tap_step_pct");
    GenericContainer::check_size(trafo_tap_pos, size, "trafo_tap_pos");
    GenericContainer::check_size(trafo_shift_degree, size, "trafo_shift_degree");
    GenericContainer::check_size(trafo_tap_hv, static_cast<std::vector<bool>::size_type>(size), "trafo_tap_hv");
    GenericContainer::check_size(trafo_hv_id, size, "trafo_hv_id");
    GenericContainer::check_size(trafo_lv_id, size, "trafo_lv_id");

    //TODO "parrallel" in the pandapower dataframe, like for lines, are not handled. Handle it python side!

    RealVect ratio = my_one_ + 0.01 * trafo_tap_step_pct.array() * trafo_tap_pos.array();

    r_ = trafo_r;
    x_ = trafo_x;
    h_ = trafo_b;
    ratio_ = ratio;
    shift_ = trafo_shift_degree / my_180_pi_;  // do not forget conversion degree / rad here !
    bus_hv_id_ = trafo_hv_id;
    bus_lv_id_ = trafo_lv_id;
    is_tap_hv_side_ = trafo_tap_hv;
    status_ = std::vector<bool>(trafo_r.size(), true);
    _update_model_coeffs();
    reset_results();

}


TrafoContainer::StateRes TrafoContainer::get_state() const
{
     std::vector<real_type> branch_r(r_.begin(), r_.end());
     std::vector<real_type> branch_x(x_.begin(), x_.end());
     std::vector<cplx_type> branch_h(h_.begin(), h_.end());
     std::vector<int > bus_hv_id(bus_hv_id_.begin(), bus_hv_id_.end());
     std::vector<int > bus_lv_id(bus_lv_id_.begin(), bus_lv_id_.end());
     std::vector<bool> status = status_;
     std::vector<real_type> ratio(ratio_.begin(), ratio_.end());
     std::vector<real_type> shift(shift_.begin(), shift_.end());
     std::vector<bool> is_tap_hv_side = is_tap_hv_side_;
     TrafoContainer::StateRes res(names_, branch_r, branch_x, branch_h, bus_hv_id, bus_lv_id, status, ratio, is_tap_hv_side, shift);
     return res;
}


void TrafoContainer::set_state(TrafoContainer::StateRes & my_state)
{
    names_ = std::get<0>(my_state);
    std::vector<real_type> & branch_r = std::get<1>(my_state);
    std::vector<real_type> & branch_x = std::get<2>(my_state);
    std::vector<cplx_type> & branch_h = std::get<3>(my_state);
    std::vector<int> & bus_hv_id = std::get<4>(my_state);
    std::vector<int> & bus_lv_id = std::get<5>(my_state);
    std::vector<bool> & status = std::get<6>(my_state);
    std::vector<real_type> & ratio = std::get<7>(my_state);
    std::vector<bool> & is_tap_hv_side = std::get<8>(my_state);
    std::vector<real_type> & shift = std::get<9>(my_state);

    auto size = branch_r.size();
    GenericContainer::check_size(branch_r, size, "branch_r");
    GenericContainer::check_size(branch_x, size, "branch_x");
    GenericContainer::check_size(branch_h, size, "branch_h");
    GenericContainer::check_size(bus_hv_id, size, "bus_hv_id");
    GenericContainer::check_size(bus_lv_id, size, "bus_lv_id");
    GenericContainer::check_size(status, size, "status");
    GenericContainer::check_size(ratio, size, "ratio");
    GenericContainer::check_size(is_tap_hv_side, size, "is_tap_hv_side");
    GenericContainer::check_size(shift, size, "shift");

    // now assign the values
    r_ = RealVect::Map(&branch_r[0], size);
    x_ = RealVect::Map(&branch_x[0], size);
    h_ = CplxVect::Map(&branch_h[0], size);

    // input data
    bus_hv_id_ = Eigen::VectorXi::Map(&bus_hv_id[0], size);
    bus_lv_id_ = Eigen::VectorXi::Map(&bus_lv_id[0], size);
    status_ = status;
    ratio_  = RealVect::Map(&ratio[0], size);
    shift_  = RealVect::Map(&shift[0], size);
    is_tap_hv_side_ = is_tap_hv_side;
    _update_model_coeffs();
    reset_results();
}


void TrafoContainer::_update_model_coeffs()
{
    const Eigen::Index my_size = r_.size();

    yac_ff_ = CplxVect::Zero(my_size);
    yac_ft_ = CplxVect::Zero(my_size);
    yac_tf_ = CplxVect::Zero(my_size);
    yac_tt_ = CplxVect::Zero(my_size);

    ydc_ff_ = CplxVect::Zero(my_size);
    ydc_ft_ = CplxVect::Zero(my_size);
    ydc_tf_ = CplxVect::Zero(my_size);
    ydc_tt_ = CplxVect::Zero(my_size);
    dc_x_tau_shift_ = RealVect::Zero(my_size);
    for(Eigen::Index i = 0; i < my_size; ++i)
    {
        // for AC
        // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.2
        const cplx_type ys = 1. / (r_(i) + my_i * x_(i));
        const cplx_type h = my_i * h_(i) * 0.5;
        double tau = ratio_(i);
        if(!is_tap_hv_side_[i]) tau = my_one_ / tau;
        real_type theta_shift = shift_(i);
        cplx_type eitheta_shift  = {my_one_, my_zero_};  // exp(j  * alpha)
        cplx_type emitheta_shift = {my_one_, my_zero_};  // exp(-j * alpha)
        if(theta_shift != 0.)
        {
            real_type cos_theta = std::cos(theta_shift);
            real_type sin_theta = std::sin(theta_shift);
            eitheta_shift = {cos_theta, sin_theta};
            emitheta_shift = {cos_theta, -sin_theta};
        }

        yac_ff_(i) = (ys + h) / (tau * tau);
        yac_tt_(i) = (ys + h);
        yac_tf_(i) = -ys / tau * emitheta_shift ;
        yac_ft_(i) = -ys / tau * eitheta_shift;

        // for DC
        // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.21
        // except here I only care about the real part, so I remove the "1/j"
        cplx_type tmp = 1. / (tau * x_(i));
        ydc_ff_(i) = tmp;
        ydc_tt_(i) = tmp;
        ydc_tf_(i) = -tmp;
        ydc_ft_(i) = -tmp;
        dc_x_tau_shift_(i) = std::real(tmp) * theta_shift;
    }
}

void TrafoContainer::fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res,
                                    bool ac,
                                    const std::vector<int> & id_grid_to_solver)
{
    throw std::runtime_error("You should not use that!");
}

void TrafoContainer::fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                              bool ac,
                              const std::vector<int> & id_grid_to_solver,
                              real_type sn_mva) const
{
    //TODO merge that with fillYbusBranch!
    //TODO template here instead of "if" for ac / dc
    const Eigen::Index nb_trafo = nb();
    cplx_type yft, ytf, yff, ytt;
    for(Eigen::Index trafo_id =0; trafo_id < nb_trafo; ++trafo_id){
        // i don't do anything if the trafo is disconnected
        if(!status_[trafo_id]) continue;

        // compute from / to
        int bus_hv_id_me = bus_hv_id_(trafo_id);
        int bus_hv_solver_id = id_grid_to_solver[bus_hv_id_me];
        if(bus_hv_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::fillYbus: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (hv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_lv_id_me = bus_lv_id_(trafo_id);
        int bus_lv_solver_id = id_grid_to_solver[bus_lv_id_me];
        if(bus_lv_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::fillYbus: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (lv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        
        if(ac){
            // ac mode
            yft = yac_ft_(trafo_id);
            ytf = yac_tf_(trafo_id);
            yff = yac_ff_(trafo_id);
            ytt = yac_tt_(trafo_id);
        }else{
            // dc mode
            yft = ydc_ft_(trafo_id);
            ytf = ydc_tf_(trafo_id);
            yff = ydc_ff_(trafo_id);
            ytt = ydc_tt_(trafo_id);
        }
        res.push_back(Eigen::Triplet<cplx_type> (bus_hv_solver_id, bus_lv_solver_id, yft));
        res.push_back(Eigen::Triplet<cplx_type> (bus_lv_solver_id, bus_hv_solver_id, ytf));
        res.push_back(Eigen::Triplet<cplx_type> (bus_hv_solver_id, bus_hv_solver_id, yff));
        res.push_back(Eigen::Triplet<cplx_type> (bus_lv_solver_id, bus_lv_solver_id, ytt));
    }
}

void TrafoContainer::hack_Sbus_for_dc_phase_shifter(CplxVect & Sbus,
                                                    bool ac,
                                                    const std::vector<int> & id_grid_to_solver)
{
    if(ac) return;
    // return;
    const int nb_trafo = nb();
    int bus_id_me, bus_id_solver_hv, bus_id_solver_lv;
    // cplx_type tmp;
    for(int trafo_id = 0; trafo_id < nb_trafo; ++trafo_id){
        //  i don't do anything if the load is disconnected
        if(!status_[trafo_id]) continue;
        if(dc_x_tau_shift_[trafo_id] == 0.) continue; // nothing to do if the trafo is not concerned (no phase shifter)
        bus_id_me = bus_lv_id_(trafo_id);
        bus_id_solver_lv = id_grid_to_solver[bus_id_me];
        if(bus_id_solver_lv == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (lv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        bus_id_me = bus_hv_id_(trafo_id);
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

void TrafoContainer::compute_results(const Eigen::Ref<const RealVect> & Va,
                                     const Eigen::Ref<const RealVect> & Vm,
                                     const Eigen::Ref<const CplxVect> & V,
                                     const std::vector<int> & id_grid_to_solver,
                                     const RealVect & bus_vn_kv,
                                     real_type sn_mva,
                                     bool ac
                                     )
{
    // it needs to be initialized at 0.
    const int nb_element = nb();
    for(int trafo_id = 0; trafo_id < nb_element; ++trafo_id){
        // don't do anything if the element is disconnected
        if(!status_[trafo_id]) {
            res_p_hv_(trafo_id) = 0.0;  // in MW
            res_q_hv_(trafo_id) = 0.0;  // in MVar
            res_v_hv_(trafo_id) = 0.0;  // in kV
            res_a_hv_(trafo_id) = 0.0;  // in kA
            res_p_lv_(trafo_id) = 0.0;  // in MW
            res_q_lv_(trafo_id) = 0.0;  // in MVar
            res_v_lv_(trafo_id) = 0.0;  // in kV
            res_a_lv_(trafo_id) = 0.0;  // in kA
            res_theta_hv_(trafo_id) = 0.0;  // in degree
            res_theta_lv_(trafo_id) = 0.0;  // in degree            
            continue;
        }

        // connectivity
        int bus_hv_id_me = bus_hv_id_(trafo_id);
        int bus_hv_solver_id = id_grid_to_solver[bus_hv_id_me];
        if(bus_hv_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::compute_results: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (hv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_lv_id_me = bus_lv_id_(trafo_id);
        int bus_lv_solver_id = id_grid_to_solver[bus_lv_id_me];
        if(bus_lv_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::compute_results: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (lv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }

        // retrieve voltages magnitude in kv instead of pu
        real_type v_hv = Vm(bus_hv_solver_id);
        real_type v_lv = Vm(bus_lv_solver_id);
        real_type bus_vn_kv_hv = bus_vn_kv(bus_hv_id_me);
        real_type bus_vn_kv_lv = bus_vn_kv(bus_lv_id_me);           

        // for voltages
        res_v_hv_(trafo_id) = v_hv * bus_vn_kv_hv;
        res_v_lv_(trafo_id) = v_lv * bus_vn_kv_lv;

        res_theta_hv_(trafo_id) = Va(bus_hv_solver_id) * my_180_pi_;
        res_theta_lv_(trafo_id) = Va(bus_lv_solver_id) * my_180_pi_;

        if(ac){
            // results of the ac powerflow
            cplx_type Ehv = V(bus_hv_solver_id);
            cplx_type Elv = V(bus_lv_solver_id);

            // TODO for DC with yff, ...
            // trafo equations
            cplx_type I_hvlv =  yac_ff_(trafo_id) * Ehv + yac_ft_(trafo_id) * Elv;
            cplx_type I_lvhv =  yac_tt_(trafo_id) * Elv + yac_tf_(trafo_id) * Ehv;

            I_hvlv = std::conj(I_hvlv);
            I_lvhv = std::conj(I_lvhv);
            cplx_type s_hvlv = Ehv * I_hvlv;
            cplx_type s_lvhv = Elv * I_lvhv;

            res_p_hv_(trafo_id) = std::real(s_hvlv) * sn_mva;
            res_q_hv_(trafo_id) = std::imag(s_hvlv) * sn_mva;
            res_p_lv_(trafo_id) = std::real(s_lvhv) * sn_mva;
            res_q_lv_(trafo_id) = std::imag(s_lvhv) * sn_mva;
        }else{
            // result of the dc powerflow
            res_p_hv_(trafo_id) = (std::real(ydc_ff_(trafo_id)) * Va(bus_hv_solver_id) + std::real(ydc_ft_(trafo_id)) * Va(bus_lv_solver_id) - dc_x_tau_shift_(trafo_id) ) * sn_mva;
            res_p_lv_(trafo_id) = (std::real(ydc_tt_(trafo_id)) * Va(bus_lv_solver_id) + std::real(ydc_tf_(trafo_id)) * Va(bus_hv_solver_id) + dc_x_tau_shift_(trafo_id) ) * sn_mva; 
           
            res_q_hv_(trafo_id) = 0.;
            res_q_lv_(trafo_id) = 0.;
            
            // for voltages, because vm = 1. pu by hypothesis
            // res_v_hv_(trafo_id) = bus_vn_kv_hv;
            // res_v_lv_(trafo_id) = bus_vn_kv_lv;
        }

    }
    _get_amps(res_a_hv_, res_p_hv_, res_q_hv_, res_v_hv_);
    _get_amps(res_a_lv_, res_p_lv_, res_q_lv_, res_v_lv_);
}

void TrafoContainer::reset_results(){
    res_p_hv_ = RealVect(nb());  // in MW
    res_q_hv_ = RealVect(nb());  // in MVar
    res_v_hv_ = RealVect(nb());  // in kV
    res_a_hv_ = RealVect(nb());  // in kA
    res_p_lv_ = RealVect(nb());  // in MW
    res_q_lv_ = RealVect(nb());  // in MVar
    res_v_lv_ = RealVect(nb());  // in kV
    res_a_lv_ = RealVect(nb());  // in kA
    res_theta_hv_ = RealVect(nb());
    res_theta_lv_ = RealVect(nb());
}


void TrafoContainer::fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                                std::vector<Eigen::Triplet<real_type> > & Bpp,
                                const std::vector<int> & id_grid_to_solver,
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

    //diagonal coefficients
    for(Eigen::Index tr_id=0; tr_id < nb_trafo; ++tr_id){
        // i only add this if the powerline is connected
        if(!status_[tr_id]) continue;

        // get the from / to bus id
        int bus_or_id_me = bus_hv_id_(tr_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::fillBp_Bpp: the trafo with id ";
            exc_ << tr_id;
            exc_ << " is connected (hv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_ex_id_me = bus_lv_id_(tr_id);
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
        double tau_bpp = is_tap_hv_side_[tr_id] ? ratio_(tr_id) : 1. / ratio_(tr_id);
        // for Bp we need shift
        const real_type theta_shift = shift_(tr_id);
        cplx_type eitheta_shift_bp  = {my_one_, my_zero_};  // exp(j  * alpha)
        cplx_type emitheta_shift_bp = {my_one_, my_zero_};  // exp(-j * alpha)
        if(theta_shift != 0.)
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
            exc_ << "TrafoContainer::fillBp_Bpp: unknown method for the FDPF powerflow for line id ";
            exc_ << tr_id;
            throw std::runtime_error(exc_.str());            
        }

        const real_type ys_bp_r = std::imag(ys_bp); 
        yff_bp = ys_bp_r;
        ytt_bp = ys_bp_r;
        ytf_bp = -std::imag(ys_bp * emitheta_shift_bp);
        yft_bp = -std::imag(ys_bp * eitheta_shift_bp);
        const real_type ys_bpp_r = std::imag(ys_bpp); 
        yff_bpp = (ys_bpp_r + std::imag(0.5 * my_i * h_(tr_id))) / (tau_bpp * tau_bpp);
        ytt_bpp = ys_bpp_r + std::imag(0.5 * my_i * h_(tr_id));
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


void TrafoContainer::fillBf_for_PTDF(std::vector<Eigen::Triplet<real_type> > & Bf,
                                     const std::vector<int> & id_grid_to_solver,
                                     real_type sn_mva,
                                     int nb_powerline,
                                     bool transpose) const
{
    const Eigen::Index nb_trafo = r_.size();

    for(Eigen::Index tr_id=0; tr_id < nb_trafo; ++tr_id){
        // i only add this if the powerline is connected
        if(!status_[tr_id]) continue;

        // get the from / to bus id
        int bus_or_id_me = bus_hv_id_(tr_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::fillBf_for_PTDF: the line with id ";
            exc_ << tr_id;
            exc_ << " is connected (hv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_ex_id_me = bus_lv_id_(tr_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "tr_id::fillBf_for_PTDF: the line with id ";
            exc_ << tr_id;
            exc_ << " is connected (lv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        real_type x = x_(tr_id);
        real_type _1_tau = is_tap_hv_side_[tr_id] ? 1. / ratio_(tr_id) : ratio_(tr_id); // 1. / tau

        // TODO
        // Bf (nb_branch, nb_bus) : en dc un truc du genre 1 / x / tap for (1..nb_branch, from_bus)
        // and -1. / x / tap for (1..nb_branch, to_bus) 
        if(transpose){
            Bf.push_back(Eigen::Triplet<real_type> (bus_or_solver_id, tr_id + nb_powerline, 1. / x * _1_tau));
            Bf.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id, tr_id + nb_powerline, -1. / x * _1_tau));
        }else{
            Bf.push_back(Eigen::Triplet<real_type> (tr_id + nb_powerline, bus_or_solver_id, 1. / x * _1_tau));
            Bf.push_back(Eigen::Triplet<real_type> (tr_id + nb_powerline, bus_ex_solver_id, -1. / x * _1_tau));
        }
    }

}


void TrafoContainer::reconnect_connected_buses(std::vector<bool> & bus_status) const{

    const Eigen::Index nb_trafo = nb();
    for(Eigen::Index trafo_id = 0; trafo_id < nb_trafo; ++trafo_id){
        // don't do anything if the element is disconnected
        if(!status_[trafo_id]) continue;
        
        const auto bus_or_id_me = bus_hv_id_(trafo_id);        
        if(bus_or_id_me == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "TrafoContainer::reconnect_connected_buses: Trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (hv) to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_trafo(...)` ?.";
            throw std::runtime_error(exc_.str());
        }
        bus_status[bus_or_id_me] = true;

        const auto bus_ex_id_me = bus_lv_id_(trafo_id);        
        if(bus_ex_id_me == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "TrafoContainer::reconnect_connected_buses: Trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (lv) to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_trafo(...)` ?.";
            throw std::runtime_error(exc_.str());
        }
        bus_status[bus_ex_id_me] = true;
    }
}

void TrafoContainer::nb_line_end(std::vector<int> & res) const{

    const Eigen::Index nb_trafo = nb();
    for(Eigen::Index trafo_id = 0; trafo_id < nb_trafo; ++trafo_id){
        // don't do anything if the element is disconnected
        if(!status_[trafo_id]) continue;
        const auto bus_or = bus_hv_id_(trafo_id);
        const auto bus_ex = bus_lv_id_(trafo_id);
        res[bus_or] += 1;
        res[bus_ex] += 1;
    }
}

void TrafoContainer::get_graph(std::vector<Eigen::Triplet<real_type> > & res) const
{
    const auto my_size = nb();
    for(Eigen::Index line_id = 0; line_id < my_size; ++line_id){
        // don't do anything if the element is disconnected
        if(!status_[line_id]) continue;
        const auto bus_or = bus_hv_id_(line_id);
        const auto bus_ex = bus_lv_id_(line_id);
        res.push_back(Eigen::Triplet<real_type>(bus_or, bus_ex, 1.));
        res.push_back(Eigen::Triplet<real_type>(bus_ex, bus_or, 1.));
    }
}

void TrafoContainer::disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component)
{
    const Eigen::Index nb_line = nb();
    SolverControl unused_solver_control;
    for(Eigen::Index i = 0; i < nb_line; ++i){
        if(!status_[i]) continue;
        auto bus_or = bus_hv_id_(i);
        auto bus_ex = bus_lv_id_(i);
        if(!busbar_in_main_component[bus_or] || !busbar_in_main_component[bus_ex]){
            deactivate(i, unused_solver_control);
        }
    }
}
