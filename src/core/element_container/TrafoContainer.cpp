// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "TrafoContainer.hpp"
#include "BinaryArchive.hpp"

#include <iostream>
#include <sstream>

namespace ls2g {

void TrafoContainer::init(
    const Eigen::Ref<const RealVect> & trafo_r,
    const Eigen::Ref<const RealVect> & trafo_x,
    const Eigen::Ref<const CplxVect> & trafo_b,
    const Eigen::Ref<const RealVect> & trafo_tap_step_pct,
    const Eigen::Ref<const RealVect> & trafo_tap_pos,
    const Eigen::Ref<const RealVect> & trafo_shift_degree,
    const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
    const Eigen::Ref<const Eigen::VectorXi> & trafo_hv_id,
    const Eigen::Ref<const Eigen::VectorXi> & trafo_lv_id,
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
    const Eigen::Ref<const RealVect> & trafo_r,
    const Eigen::Ref<const RealVect> & trafo_x,
    const Eigen::Ref<const CplxVect> & trafo_b,
    const Eigen::Ref<const RealVect> & trafo_ratio,
    const Eigen::Ref<const RealVect> & trafo_shift_degree,
    const std::vector<bool> & trafo_tap_side1,  // is tap on high voltage (true) or low voltate
    const Eigen::Ref<const Eigen::VectorXi> & trafo_hv_id,
    const Eigen::Ref<const Eigen::VectorXi> & trafo_lv_id,
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
    // no alpha-dependent r/x correction by default (neutral = the given r/x); the
    // pypowsybl converter enables it afterwards via set_shift_dependent_rx
    base_r_ = trafo_r;
    base_x_ = trafo_x;
    shift_dependent_rx_ = false;
    rx_corr_alpha_ = std::vector<std::vector<real_type> >(size, std::vector<real_type>());
    rx_corr_pct_ = std::vector<std::vector<real_type> >(size, std::vector<real_type>());
    init_tsc(trafo_hv_id, trafo_lv_id, "trafo");
    _update_model_coeffs();
    reset_results();
}

TrafoContainer::StateRes TrafoContainer::get_state() const
{
     std::vector<real_type> ratio(ratio_.begin(), ratio_.end());
     std::vector<real_type> shift(shift_.begin(), shift_.end());
     std::vector<bool> is_tap_hv_side = is_tap_side1_;
     std::vector<real_type> base_r(base_r_.begin(), base_r_.end());
     std::vector<real_type> base_x(base_x_.begin(), base_x_.end());
     TrafoContainer::StateRes res(
        get_tsc_rxha_state(),
        ratio,
        is_tap_hv_side,
        shift,
        ignore_tap_side_for_shift_,
        shift_dependent_rx_,
        base_r,
        base_x,
        rx_corr_alpha_,
        rx_corr_pct_);
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

    ratio_  = RealVect::Map(ratio.data(), size);
    shift_  = RealVect::Map(shift.data(), size);
    is_tap_side1_ = is_tap_side1;
    ignore_tap_side_for_shift_ = std::get<4>(my_state);

    shift_dependent_rx_ = std::get<5>(my_state);
    std::vector<real_type> & base_r = std::get<6>(my_state);
    std::vector<real_type> & base_x = std::get<7>(my_state);
    GenericContainer::check_size(base_r, size, "base_r");
    GenericContainer::check_size(base_x, size, "base_x");
    base_r_ = RealVect::Map(base_r.data(), size);
    base_x_ = RealVect::Map(base_x.data(), size);

    // The two alpha -> r/x correction tables come straight from a pickle or a binary
    // file. `_update_model_coeffs()` right below walks el_id over [0, nb()) and reads
    // `rx_corr_alpha_[el_id]` / `rx_corr_pct_[el_id]` with an unchecked
    // std::vector::operator[], and `_shift_rx_corr_pct` then interpolates `ys` using
    // indices derived from `xs.size()`. So a state declaring fewer tables than
    // transformers builds a std::vector out of whatever the heap holds past the end
    // (and dereferences its pointers), and a state whose alpha / correction tables
    // differ in length reads past the end of the shorter one. Neither is caught by
    // check_grid(): this runs *before* it. Validate both shapes here, exactly the
    // invariant init() and set_shift_dependent_rx() maintain.
    const std::vector<std::vector<real_type> > & rx_corr_alpha = std::get<8>(my_state);
    const std::vector<std::vector<real_type> > & rx_corr_pct = std::get<9>(my_state);
    if(rx_corr_alpha.empty() && rx_corr_pct.empty()){
        // no table at all: only legal when nothing would index them
        if(shift_dependent_rx_ && (size > 0)){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::set_state: shift-dependent r/x is enabled for " << size
                 << " transformer(s) but the state carries no (alpha -> correction) table at all. "
                 << "The tables are indexed by transformer id, so this state is inconsistent.";
            throw std::runtime_error(exc_.str());
        }
    } else {
        GenericContainer::check_size(rx_corr_alpha, size, "rx_corr_alpha");
        GenericContainer::check_size(rx_corr_pct, size, "rx_corr_pct");
        for(std::size_t el_id = 0; el_id < rx_corr_alpha.size(); ++el_id){
            if(rx_corr_alpha[el_id].size() != rx_corr_pct[el_id].size()){
                std::ostringstream exc_;
                exc_ << "TrafoContainer::set_state: transformer " << el_id << " has "
                     << rx_corr_alpha[el_id].size() << " alpha sample(s) but "
                     << rx_corr_pct[el_id].size() << " r/x correction value(s). Both tables are "
                     << "read together (one correction per alpha) and must have the same length.";
                throw std::runtime_error(exc_.str());
            }
        }
    }
    rx_corr_alpha_ = rx_corr_alpha;
    rx_corr_pct_ = rx_corr_pct;

    _update_model_coeffs();
    reset_results();
}

void TrafoContainer::set_shift_dependent_rx(
    bool enable,
    const std::vector<std::vector<real_type> > & alpha_rad,
    const std::vector<std::vector<real_type> > & rx_corr_pct,
    DualAlgoControl & solver_control)
{
    const auto size = nb();
    if(alpha_rad.size() != static_cast<std::size_t>(size))
        throw std::runtime_error("TrafoContainer::set_shift_dependent_rx: alpha_rad has a wrong size");
    if(rx_corr_pct.size() != static_cast<std::size_t>(size))
        throw std::runtime_error("TrafoContainer::set_shift_dependent_rx: rx_corr_pct has a wrong size");
    shift_dependent_rx_ = enable;
    rx_corr_alpha_ = alpha_rad;
    rx_corr_pct_ = rx_corr_pct;
    // the stored r_ / x_ at this point are the neutral impedance: keep them as the base
    base_r_ = r_;
    base_x_ = x_;
    // sort each (alpha -> correction) table by ascending alpha so the interpolation
    // in _shift_rx_corr_pct is well defined
    for(size_t el_id = 0; el_id < size; ++el_id){
        auto & xs = rx_corr_alpha_[el_id];
        auto & ys = rx_corr_pct_[el_id];
        if(xs.size() != ys.size())
            throw std::runtime_error("TrafoContainer::set_shift_dependent_rx: alpha and correction tables differ in size");
        std::vector<std::size_t> order(xs.size());
        for(std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&xs](std::size_t a, std::size_t b){return xs[a] < xs[b];});
        std::vector<real_type> sx(xs.size()), sy(ys.size());
        for(std::size_t i = 0; i < order.size(); ++i){ sx[i] = xs[order[i]]; sy[i] = ys[order[i]]; }
        xs.swap(sx); ys.swap(sy);
    }
    // re-apply the (possibly corrected) impedance at the current shift
    _update_model_coeffs();
    solver_control.ac_algo_controler().tell_recompute_ybus();
    solver_control.dc_algo_controler().tell_recompute_ybus();
}

void TrafoContainer::_update_model_coeffs_one_el(int el_id)
{
    // phase-shifting transformers whose series impedance depends on the phase-shift
    // angle: refresh the effective r / x from the neutral value and the correction
    // interpolated at the current `shift_`. This makes `change_shift` / `change_ratio`
    // (which call `_update_internal_coeffs`) keep r / x in sync with no "tap" concept.
    if(shift_dependent_rx_ && !rx_corr_alpha_[el_id].empty()){
        const real_type corr = my_one_ + _shift_rx_corr_pct(el_id) / 100.;
        r_(el_id) = base_r_(el_id) * corr;
        x_(el_id) = base_x_(el_id) * corr;
    }

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
    // except here I only care about the real part (1/x), so I remove the "1/j"
    const real_type tmp = 1. / x_(el_id) * _1_tau;
    ydc_11_(el_id) = tmp;
    ydc_22_(el_id) = tmp;
    ydc_21_(el_id) = -tmp;
    ydc_12_(el_id) = -tmp;

    dc_x_tau_shift_(el_id) = -tmp * theta_shift;
}

void TrafoContainer::hack_Sbus_for_dc_phase_shifter(
    Eigen::Ref<CplxVect> Sbus,
    bool ac,
    const SolverBusIdVect & id_grid_to_solver)
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
#ifndef NDEBUG
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: (GridModelId) the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (side 2) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
#endif
        bus_id_solver_lv = id_grid_to_solver[bus_id_me.cast_int()];
#ifndef NDEBUG
        if(bus_id_solver_lv.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: (SolverId) the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (side 2) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
#endif

        bus_id_me = get_bus_side_1(trafo_id);
#ifndef NDEBUG
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: (GridModelId) the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (side 1) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
#endif
        bus_id_solver_hv = id_grid_to_solver[bus_id_me.cast_int()];
#ifndef NDEBUG
        if(bus_id_solver_hv.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "TrafoContainer::hack_Sbus_for_dc_phase_shifter: (SolverId) the trafo with id ";
            exc_ << trafo_id;
            exc_ << " is connected (side 1) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
#endif
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

void TrafoContainer::save_binary(const std::string & path, bool atomic) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
}

TrafoContainer TrafoContainer::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<TrafoContainer>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

} // namespace ls2g
