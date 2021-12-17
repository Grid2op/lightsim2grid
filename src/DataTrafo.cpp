// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DataTrafo.h"
#include <iostream>
#include <sstream>

void DataTrafo::init(const RealVect & trafo_r,
                           const RealVect & trafo_x,
                           const CplxVect & trafo_b,
                           const RealVect & trafo_tap_step_pct,
            //                        const RealVect & trafo_tap_step_degree,
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
    DataGeneric::check_size(trafo_r, size, "trafo_r");
    DataGeneric::check_size(trafo_x, size, "trafo_x");
    DataGeneric::check_size(trafo_b, size, "trafo_b");
    DataGeneric::check_size(trafo_tap_step_pct, size, "trafo_tap_step_pct");
    DataGeneric::check_size(trafo_tap_pos, size, "trafo_tap_pos");
    DataGeneric::check_size(trafo_shift_degree, size, "trafo_shift_degree");
    DataGeneric::check_size(trafo_tap_hv, static_cast<std::vector<bool>::size_type>(size), "trafo_tap_hv");
    DataGeneric::check_size(trafo_hv_id, size, "trafo_hv_id");
    DataGeneric::check_size(trafo_lv_id, size, "trafo_lv_id");

    //TODO "parrallel" in the pandapower dataframe, like for lines, are not handled. Handle it python side!

    RealVect ratio = my_one_ + 0.01 * trafo_tap_step_pct.array() * trafo_tap_pos.array();

    r_ = trafo_r;
    x_ = trafo_x;
    h_ = trafo_b;
    ratio_ = ratio;
    shift_ = trafo_shift_degree / 180. * my_pi;  // do not forget conversion degree / rad here !
    bus_hv_id_ = trafo_hv_id;
    bus_lv_id_ = trafo_lv_id;
    is_tap_hv_side_ = trafo_tap_hv;
    status_ = std::vector<bool>(trafo_r.size(), true);
    _update_model_coeffs();

}

DataTrafo::StateRes DataTrafo::get_state() const
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
     DataTrafo::StateRes res(branch_r, branch_x, branch_h, bus_hv_id, bus_lv_id, status, ratio, is_tap_hv_side, shift);
     return res;
}
void DataTrafo::set_state(DataTrafo::StateRes & my_state)
{
    reset_results();

    std::vector<real_type> & branch_r = std::get<0>(my_state);
    std::vector<real_type> & branch_x = std::get<1>(my_state);
    std::vector<cplx_type> & branch_h = std::get<2>(my_state);
    std::vector<int> & bus_hv_id = std::get<3>(my_state);
    std::vector<int> & bus_lv_id = std::get<4>(my_state);
    std::vector<bool> & status = std::get<5>(my_state);
    std::vector<real_type> & ratio = std::get<6>(my_state);
    std::vector<bool> & is_tap_hv_side = std::get<7>(my_state);
    std::vector<real_type> & shift = std::get<8>(my_state);

    auto size = branch_r.size();
    DataGeneric::check_size(branch_r, size, "branch_r");
    DataGeneric::check_size(branch_x, size, "branch_x");
    DataGeneric::check_size(branch_h, size, "branch_h");
    DataGeneric::check_size(bus_hv_id, size, "bus_hv_id");
    DataGeneric::check_size(bus_lv_id, size, "bus_lv_id");
    DataGeneric::check_size(status, size, "status");
    DataGeneric::check_size(ratio, size, "ratio");
    DataGeneric::check_size(is_tap_hv_side, size, "is_tap_hv_side");
    DataGeneric::check_size(shift, size, "shift");

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
}

void DataTrafo::_update_model_coeffs()
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

void DataTrafo::fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver)
{
    throw std::runtime_error("You should not use that!");
}

void DataTrafo::fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                         bool ac,
                         const std::vector<int> & id_grid_to_solver,
                         real_type sn_mva)
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
            exc_ << "DataTrafo::fillYbus: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (hv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_lv_id_me = bus_lv_id_(trafo_id);
        int bus_lv_solver_id = id_grid_to_solver[bus_lv_id_me];
        if(bus_lv_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "DataTrafo::fillYbus: the trafo with id ";
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

void DataTrafo::hack_Sbus_for_dc_phase_shifter(CplxVect & Sbus, bool ac, const std::vector<int> & id_grid_to_solver){
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
            exc_ << "DataTrafo::fillSbus: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (lv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        bus_id_me = bus_hv_id_(trafo_id);
        bus_id_solver_hv = id_grid_to_solver[bus_id_me];
        if(bus_id_solver_hv == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "DataTrafo::fillSbus: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (hv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        Sbus.coeffRef(bus_id_solver_hv) += dc_x_tau_shift_[trafo_id];
        Sbus.coeffRef(bus_id_solver_lv) -= dc_x_tau_shift_[trafo_id];
    }
}

void DataTrafo::compute_results(const Eigen::Ref<const RealVect> & Va,
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
    res_p_hv_ = RealVect::Constant(nb_element, 0.0);  // in MW
    res_q_hv_ = RealVect::Constant(nb_element, 0.0);  // in MVar
    res_v_hv_ = RealVect::Constant(nb_element, 0.0);  // in kV
    res_a_hv_ = RealVect::Constant(nb_element, 0.0);  // in kA
    res_p_lv_ = RealVect::Constant(nb_element, 0.0);  // in MW
    res_q_lv_ = RealVect::Constant(nb_element, 0.0);  // in MVar
    res_v_lv_ = RealVect::Constant(nb_element, 0.0);  // in kV
    res_a_lv_ = RealVect::Constant(nb_element, 0.0);  // in kA
    res_theta_hv_ = RealVect::Constant(nb_element, 0.0);  // in degree
    res_theta_lv_ = RealVect::Constant(nb_element, 0.0);  // in degree
    for(int trafo_id = 0; trafo_id < nb_element; ++trafo_id){
        // don't do anything if the element is disconnected
        if(!status_[trafo_id]) continue;

        // connectivity
        int bus_hv_id_me = bus_hv_id_(trafo_id);
        int bus_hv_solver_id = id_grid_to_solver[bus_hv_id_me];
        if(bus_hv_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "DataTrafo::compute_results: the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (hv side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_lv_id_me = bus_lv_id_(trafo_id);
        int bus_lv_solver_id = id_grid_to_solver[bus_lv_id_me];
        if(bus_lv_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "DataTrafo::compute_results: the trafo with id ";
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

        res_theta_hv_(trafo_id) = Va(bus_hv_solver_id) * 180. / my_pi;
        res_theta_lv_(trafo_id) = Va(bus_lv_solver_id) * 180. / my_pi;

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

            // for voltages, because vm = 1. pu by hypothesis
            // res_v_hv_(trafo_id) = bus_vn_kv_hv;
            // res_v_lv_(trafo_id) = bus_vn_kv_lv;
        }

    }
    _get_amps(res_a_hv_, res_p_hv_, res_q_hv_, res_v_hv_);
    _get_amps(res_a_lv_, res_p_lv_, res_q_lv_, res_v_lv_);
}

void DataTrafo::reset_results(){
    res_p_hv_ = RealVect();  // in MW
    res_q_hv_ = RealVect();  // in MVar
    res_v_hv_ = RealVect();  // in kV
    res_a_hv_ = RealVect();  // in kA
    res_p_lv_ = RealVect();  // in MW
    res_q_lv_ = RealVect();  // in MVar
    res_v_lv_ = RealVect();  // in kV
    res_a_lv_ = RealVect();  // in kA
}
