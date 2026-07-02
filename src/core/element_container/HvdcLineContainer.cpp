// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "HvdcLineContainer.hpp"
#include "BinaryArchive.hpp"

#include <cmath>
#include <sstream>

namespace ls2g {

void HvdcLineContainer::init(const Eigen::VectorXi & bus1_id,
                             const Eigen::VectorXi & bus2_id,
                             const std::vector<int> & type1,
                             const std::vector<int> & type2,
                             const RealVect & loss_factor1,
                             const RealVect & loss_factor2,
                             const std::vector<bool> & vreg1_on,
                             const std::vector<bool> & vreg2_on,
                             const RealVect & vm1_pu,
                             const RealVect & vm2_pu,
                             const RealVect & q1_setpoint_mvar,
                             const RealVect & q2_setpoint_mvar,
                             const RealVect & min_q1,
                             const RealVect & max_q1,
                             const RealVect & min_q2,
                             const RealVect & max_q2,
                             const RealVect & power_factor1,
                             const RealVect & power_factor2,
                             const std::vector<int> & converters_mode,
                             const RealVect & p_setpoint_mw,
                             const RealVect & r_ohm,
                             const RealVect & nominal_v_kv,
                             const RealVect & loss_percent,
                             const RealVect & loss_mw,
                             const std::vector<bool> & droop_enabled,
                             const RealVect & droop_p0_mw,
                             const RealVect & droop_mw_per_deg,
                             const RealVect & pmax_1to2_mw,
                             const RealVect & pmax_2to1_mw
                             )
{
    init_tsc(bus1_id, bus2_id, "hvdc lines");
    const int size = static_cast<int>(nb());
    check_size(converters_mode, size, "converters_mode");
    check_size(p_setpoint_mw, size, "p_setpoint_mw");
    check_size(r_ohm, size, "r_ohm");
    check_size(nominal_v_kv, size, "nominal_v_kv");
    check_size(loss_percent, size, "loss_percent");
    check_size(loss_mw, size, "loss_mw");
    check_size(droop_enabled, size, "droop_enabled");
    check_size(droop_p0_mw, size, "droop_p0_mw");
    check_size(droop_mw_per_deg, size, "droop_mw_per_deg");
    check_size(pmax_1to2_mw, size, "pmax_1to2_mw");
    check_size(pmax_2to1_mw, size, "pmax_2to1_mw");

    side_1_.init(type1, loss_factor1, vreg1_on, vm1_pu, q1_setpoint_mvar,
                 min_q1, max_q1, power_factor1, bus1_id);
    side_2_.init(type2, loss_factor2, vreg2_on, vm2_pu, q2_setpoint_mvar,
                 min_q2, max_q2, power_factor2, bus2_id);

    loss_percent_ = loss_percent;
    loss_mw_ = loss_mw;
    r_ohm_ = r_ohm;
    nominal_v_kv_ = nominal_v_kv;
    p_setpoint_mw_ = p_setpoint_mw;
    converters_mode_ = IntVect(size);
    status_droop_ = IntVect::Zero(size);
    droop_enabled_ = droop_enabled;
    p0_mw_ = droop_p0_mw;
    // MW / degree -> MW / radian (my_180_pi_ = 180 / pi)
    k_mw_per_rad_ = droop_mw_per_deg * BaseConstants::my_180_pi_;
    pmax_1to2_mw_ = pmax_1to2_mw;
    pmax_2to1_mw_ = pmax_2to1_mw;

    DualAlgoControl unused_solver_control;
    for(int hvdc_id = 0; hvdc_id < size; ++hvdc_id){
        if((converters_mode[hvdc_id] != ConvertersMode::SIDE_1_RECTIFIER) &&
           (converters_mode[hvdc_id] != ConvertersMode::SIDE_2_RECTIFIER)){
            std::ostringstream exc_;
            exc_ << "HvdcLineContainer::init: unknown converters_mode ";
            exc_ << converters_mode[hvdc_id];
            exc_ << " for hvdc line ";
            exc_ << hvdc_id;
            exc_ << " (0 = side 1 rectifier, 1 = side 2 rectifier).";
            throw std::runtime_error(exc_.str());
        }
        converters_mode_(hvdc_id) = converters_mode[hvdc_id];
        if(p_setpoint_mw(hvdc_id) < 0.){
            std::ostringstream exc_;
            exc_ << "HvdcLineContainer::init: p_setpoint_mw should be >= 0 (the direction is given by converters_mode) for hvdc line ";
            exc_ << hvdc_id;
            throw std::runtime_error(exc_.str());
        }
        if(droop_enabled[hvdc_id]){
            if(side_1_.is_lcc(hvdc_id) || side_2_.is_lcc(hvdc_id)){
                std::ostringstream exc_;
                exc_ << "HvdcLineContainer::init: the angle-droop (AC emulation) is only available for VSC - VSC hvdc lines (check hvdc line ";
                exc_ << hvdc_id;
                exc_ << ").";
                throw std::runtime_error(exc_.str());
            }
            if((abs(loss_percent(hvdc_id)) > _tol_equal_float) || (abs(loss_mw(hvdc_id)) > _tol_equal_float)){
                std::ostringstream exc_;
                exc_ << "HvdcLineContainer::init: the legacy pandapower losses (loss_percent / loss_mw) cannot be combined with the angle-droop ";
                exc_ << "(use the station loss factors and the line resistance instead, check hvdc line ";
                exc_ << hvdc_id;
                exc_ << ").";
                throw std::runtime_error(exc_.str());
            }
        }
        update_station_targets(hvdc_id, unused_solver_control);
    }
}

void HvdcLineContainer::init_legacy(const Eigen::VectorXi & branch_from_id,
                                    const Eigen::VectorXi & branch_to_id,
                                    const RealVect & p_mw,
                                    const RealVect & loss_percent,
                                    const RealVect & loss_mw,
                                    const RealVect & vm_or_pu,
                                    const RealVect & vm_ex_pu,
                                    const RealVect & min_q_or,
                                    const RealVect & max_q_or,
                                    const RealVect & min_q_ex,
                                    const RealVect & max_q_ex)
{
    const int size = static_cast<int>(p_mw.size());
    const std::vector<int> type_vsc(size, ConverterStationContainer::ConverterType::VSC);
    const RealVect zeros = RealVect::Zero(size);
    const RealVect ones = RealVect::Constant(size, 1.);
    const std::vector<bool> vreg_on(size, true);
    const std::vector<bool> no_droop(size, false);
    const std::vector<int> mode_side_1(size, ConvertersMode::SIDE_1_RECTIFIER);

    init(branch_from_id, branch_to_id,
         type_vsc, type_vsc,
         zeros, zeros,          // station loss factors (the legacy loss is line-level)
         vreg_on, vreg_on,
         vm_or_pu, vm_ex_pu,
         zeros, zeros,          // q setpoints (unused: regulating)
         min_q_or, max_q_or, min_q_ex, max_q_ex,
         ones, ones,            // power factors (unused: VSC)
         mode_side_1, zeros,    // overwritten right below by update_targets_from_p1
         zeros, zeros,          // r_ohm, nominal_v_kv
         loss_percent, loss_mw,
         no_droop, zeros, zeros, zeros, zeros);

    DualAlgoControl unused_solver_control;
    for(int hvdc_id = 0; hvdc_id < size; ++hvdc_id){
        update_targets_from_p1(hvdc_id, p_mw(hvdc_id), unused_solver_control);
    }
}

HvdcLineContainer::StateRes HvdcLineContainer::get_state() const
{
    std::vector<real_type> loss_percent(loss_percent_.begin(), loss_percent_.end());
    std::vector<real_type> loss_mw(loss_mw_.begin(), loss_mw_.end());
    std::vector<int> converters_mode(converters_mode_.begin(), converters_mode_.end());
    std::vector<real_type> p_setpoint_mw(p_setpoint_mw_.begin(), p_setpoint_mw_.end());
    std::vector<real_type> r_ohm(r_ohm_.begin(), r_ohm_.end());
    std::vector<real_type> nominal_v_kv(nominal_v_kv_.begin(), nominal_v_kv_.end());
    std::vector<real_type> p0_mw(p0_mw_.begin(), p0_mw_.end());
    std::vector<real_type> k_mw_per_rad(k_mw_per_rad_.begin(), k_mw_per_rad_.end());
    std::vector<real_type> pmax_1to2_mw(pmax_1to2_mw_.begin(), pmax_1to2_mw_.end());
    std::vector<real_type> pmax_2to1_mw(pmax_2to1_mw_.begin(), pmax_2to1_mw_.end());
    std::vector<int> status_droop(status_droop_.begin(), status_droop_.end());
    HvdcLineContainer::StateRes res(get_tsc_state(),
                                    loss_percent,
                                    loss_mw,
                                    converters_mode,
                                    p_setpoint_mw,
                                    r_ohm,
                                    nominal_v_kv,
                                    droop_enabled_,
                                    p0_mw,
                                    k_mw_per_rad,
                                    pmax_1to2_mw,
                                    pmax_2to1_mw,
                                    status_droop);
    return res;
}

void HvdcLineContainer::set_state(HvdcLineContainer::StateRes & my_state)
{
    set_tsc_state(std::get<StateResIdx::TSC_STATE>(my_state));
    std::vector<real_type> & loss_percent = std::get<StateResIdx::LOSS_PERCENT>(my_state);
    std::vector<real_type> & loss_mw = std::get<StateResIdx::LOSS_MW>(my_state);
    std::vector<int> & converters_mode = std::get<StateResIdx::CONVERTERS_MODE>(my_state);
    std::vector<real_type> & p_setpoint_mw = std::get<StateResIdx::P_SETPOINT_MW>(my_state);
    std::vector<real_type> & r_ohm = std::get<StateResIdx::R_OHM>(my_state);
    std::vector<real_type> & nominal_v_kv = std::get<StateResIdx::NOMINAL_V_KV>(my_state);
    std::vector<bool> & droop_enabled = std::get<StateResIdx::DROOP_ENABLED>(my_state);
    std::vector<real_type> & p0_mw = std::get<StateResIdx::P0_MW>(my_state);
    std::vector<real_type> & k_mw_per_rad = std::get<StateResIdx::K_MW_PER_RAD>(my_state);
    std::vector<real_type> & pmax_1to2_mw = std::get<StateResIdx::PMAX_1TO2_MW>(my_state);
    std::vector<real_type> & pmax_2to1_mw = std::get<StateResIdx::PMAX_2TO1_MW>(my_state);
    std::vector<int> & status_droop = std::get<StateResIdx::STATUS_DROOP>(my_state);

    const auto size = nb();
    check_size(loss_percent, size, "loss_percent");
    check_size(loss_mw, size, "loss_mw");
    check_size(converters_mode, size, "converters_mode");
    check_size(p_setpoint_mw, size, "p_setpoint_mw");
    check_size(r_ohm, size, "r_ohm");
    check_size(nominal_v_kv, size, "nominal_v_kv");
    check_size(droop_enabled, size, "droop_enabled");
    check_size(p0_mw, size, "p0_mw");
    check_size(k_mw_per_rad, size, "k_mw_per_rad");
    check_size(pmax_1to2_mw, size, "pmax_1to2_mw");
    check_size(pmax_2to1_mw, size, "pmax_2to1_mw");
    check_size(status_droop, size, "status_droop");

    loss_percent_ = RealVect::Map(&loss_percent[0], loss_percent.size());
    loss_mw_ = RealVect::Map(&loss_mw[0], loss_mw.size());
    converters_mode_ = IntVect::Map(&converters_mode[0], converters_mode.size());
    p_setpoint_mw_ = RealVect::Map(&p_setpoint_mw[0], p_setpoint_mw.size());
    r_ohm_ = RealVect::Map(&r_ohm[0], r_ohm.size());
    nominal_v_kv_ = RealVect::Map(&nominal_v_kv[0], nominal_v_kv.size());
    droop_enabled_ = droop_enabled;
    p0_mw_ = RealVect::Map(&p0_mw[0], p0_mw.size());
    k_mw_per_rad_ = RealVect::Map(&k_mw_per_rad[0], k_mw_per_rad.size());
    pmax_1to2_mw_ = RealVect::Map(&pmax_1to2_mw[0], pmax_1to2_mw.size());
    pmax_2to1_mw_ = RealVect::Map(&pmax_2to1_mw[0], pmax_2to1_mw.size());
    status_droop_ = IntVect::Map(&status_droop[0], status_droop.size());
    reset_results();
}

real_type HvdcLineContainer::recv_mw(int hvdc_id, real_type p_ctrl_abs_mw, bool side_1_is_controller) const
{
    const real_type lf_ctrl = side_1_is_controller ? side_1_.get_loss_factor(hvdc_id) : side_2_.get_loss_factor(hvdc_id);
    const real_type lf_recv = side_1_is_controller ? side_2_.get_loss_factor(hvdc_id) : side_1_.get_loss_factor(hvdc_id);
    const real_type line_in = (1. - lf_ctrl) * (1. - 0.01 * loss_percent_(hvdc_id)) * p_ctrl_abs_mw;
    real_type line_loss = 0.;
    const real_type v_nom = nominal_v_kv_(hvdc_id);
    if((v_nom > 0.) && (r_ohm_(hvdc_id) != 0.)){
        line_loss = r_ohm_(hvdc_id) * line_in * line_in / (v_nom * v_nom);
    }
    return (1. - lf_recv) * (line_in - line_loss) - loss_mw_(hvdc_id);
}

void HvdcLineContainer::droop_flows_mw(int hvdc_id, real_type raw_mw, real_type & p1_flow_mw, real_type & p2_flow_mw) const
{
    const int status = status_droop_(hvdc_id);
    if(status == 0){
        // linear regime: the flow follows the angle-droop equation
        if(raw_mw >= 0.){
            p1_flow_mw = raw_mw;
            p2_flow_mw = -recv_mw(hvdc_id, raw_mw, true);
        }else{
            p1_flow_mw = -recv_mw(hvdc_id, -raw_mw, false);
            p2_flow_mw = -raw_mw;
        }
    }else if(status > 0){
        // saturated 1 -> 2: pinned at the limit, no theta dependence
        const real_type pmax = pmax_1to2_mw_(hvdc_id);
        p1_flow_mw = pmax;
        p2_flow_mw = -recv_mw(hvdc_id, pmax, true);
    }else{
        // saturated 2 -> 1
        const real_type pmax = pmax_2to1_mw_(hvdc_id);
        p1_flow_mw = -recv_mw(hvdc_id, pmax, false);
        p2_flow_mw = pmax;
    }
}

void HvdcLineContainer::update_station_targets(int hvdc_id, DualAlgoControl & solver_control)
{
    const real_type setpoint = p_setpoint_mw_(hvdc_id);
    real_type t1, t2;
    if(converters_mode_(hvdc_id) == ConvertersMode::SIDE_1_RECTIFIER){
        t1 = -setpoint;
        t2 = recv_mw(hvdc_id, setpoint, true);
    }else{
        t2 = -setpoint;
        t1 = recv_mw(hvdc_id, setpoint, false);
    }
    side_1_.set_station_p(hvdc_id, t1, solver_control);
    side_2_.set_station_p(hvdc_id, t2, solver_control);
}

void HvdcLineContainer::update_targets_from_p1(int hvdc_id, real_type p1_mw, DualAlgoControl & solver_control)
{
    if(p1_mw < 0.){
        // side 1 consumes |p1|: side 1 is the rectifier
        converters_mode_(hvdc_id) = ConvertersMode::SIDE_1_RECTIFIER;
        p_setpoint_mw_(hvdc_id) = -p1_mw;
    }else{
        // side 1 receives p1 (this includes p1 == 0, to match the legacy
        // DCLineContainer::get_to_mw branch `from_mw >= 0` exactly): side 2
        // is the rectifier and its draw is the inverse of the loss chain
        converters_mode_(hvdc_id) = ConvertersMode::SIDE_2_RECTIFIER;
        const real_type lf_inv = side_1_.get_loss_factor(hvdc_id);
        const real_type lf_rect = side_2_.get_loss_factor(hvdc_id);
        // line_in - r' * line_in^2 = (p1 + loss_mw) / (1 - lf_inv)
        const real_type rhs = (p1_mw + loss_mw_(hvdc_id)) / (1. - lf_inv);
        const real_type v_nom = nominal_v_kv_(hvdc_id);
        const real_type r_prime = ((v_nom > 0.) && (r_ohm_(hvdc_id) != 0.)) ?
                                  r_ohm_(hvdc_id) / (v_nom * v_nom) : 0.;
        real_type line_in;
        if(r_prime == 0.){
            line_in = rhs;
        }else{
            const real_type disc = 1. - 4. * r_prime * rhs;
            if(disc < 0.){
                std::ostringstream exc_;
                exc_ << "HvdcLineContainer::change_p: the requested receive power ";
                exc_ << p1_mw;
                exc_ << " MW cannot be transmitted through the dc line resistance of hvdc line ";
                exc_ << hvdc_id;
                throw std::runtime_error(exc_.str());
            }
            line_in = (1. - std::sqrt(disc)) / (2. * r_prime);
        }
        p_setpoint_mw_(hvdc_id) = line_in / ((1. - lf_rect) * (1. - 0.01 * loss_percent_(hvdc_id)));
    }
    update_station_targets(hvdc_id, solver_control);
}

void HvdcLineContainer::change_p(int hvdc_id, real_type new_p, DualAlgoControl & solver_control)
{
    if(!status_global_.at(hvdc_id))
    {
        std::ostringstream exc_;
        exc_ << "HvdcLineContainer::change_p: Impossible to change the active value of a disconnected hvdc line (check hvdc line id ";
        exc_ << hvdc_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    update_targets_from_p1(hvdc_id, new_p, solver_control);
}

void HvdcLineContainer::set_status_droop(int hvdc_id, int status, DualAlgoControl & solver_control)
{
    _check_in_range(hvdc_id, status_droop_, "set_status_droop");
    if((status < -1) || (status > 1)){
        std::ostringstream exc_;
        exc_ << "HvdcLineContainer::set_status_droop: status should be -1 (saturated 2->1), 0 (linear) or 1 (saturated 1->2), found ";
        exc_ << status;
        throw std::runtime_error(exc_.str());
    }
    if(!droop_enabled_[hvdc_id]){
        std::ostringstream exc_;
        exc_ << "HvdcLineContainer::set_status_droop: hvdc line ";
        exc_ << hvdc_id;
        exc_ << " has no angle-droop (AC emulation) enabled.";
        throw std::runtime_error(exc_.str());
    }
    if(status_droop_(hvdc_id) != status){
        status_droop_(hvdc_id) = status;
        // the sparsity pattern of the jacobian does not change (the droop
        // entries are declared whatever the regime), only the values /
        // injections do: the symbolic factorization is fully reused.
        solver_control.ac_algo_controler().tell_recompute_sbus();
        solver_control.dc_algo_controler().tell_recompute_sbus();
    }
}

void HvdcLineContainer::fillSbus(CplxVect & Sbus, const SolverBusIdVect & id_grid_to_solver, bool ac) const
{
    const int nb_hvdc = static_cast<int>(nb());
    // for droop lines, the active power is NOT a fixed injection:
    //  - in AC it is handled by the Hvdc extension of the NR system (all regimes);
    //  - in DC the linear regime is handled by the DC algorithm and the
    //    saturated regimes are stamped below as fixed injections.
    std::vector<bool> skip_p(nb_hvdc, false);
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        if(droop_enabled_[hvdc_id] && status_global_[hvdc_id]) skip_p[hvdc_id] = true;
    }
    side_1_.fillSbus_station(Sbus, id_grid_to_solver, ac, skip_p);
    side_2_.fillSbus_station(Sbus, id_grid_to_solver, ac, skip_p);

    if(!ac){
        for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
            if(!droop_enabled_[hvdc_id] || !status_global_[hvdc_id]) continue;
            if(status_droop_(hvdc_id) == 0) continue;  // linear: theta-dependent, handled by the DC algorithm
            real_type p1_flow, p2_flow;
            droop_flows_mw(hvdc_id, 0., p1_flow, p2_flow);  // raw unused when saturated
            const GlobalBusId bus_1 = get_bus_side_1(hvdc_id);
            const GlobalBusId bus_2 = get_bus_side_2(hvdc_id);
            const SolverBusId bus_1_solver = id_grid_to_solver[bus_1.cast_int()];
            const SolverBusId bus_2_solver = id_grid_to_solver[bus_2.cast_int()];
            if((bus_1_solver.cast_int() == _deactivated_bus_id) ||
               (bus_2_solver.cast_int() == _deactivated_bus_id)){
                std::ostringstream exc_;
                exc_ << "HvdcLineContainer::fillSbus: hvdc line with id ";
                exc_ << hvdc_id;
                exc_ << " is connected to a disconnected bus while being connected to the grid.";
                throw std::runtime_error(exc_.str());
            }
            // the flows LEAVE the buses: injection = -flow
            Sbus.coeffRef(bus_1_solver.cast_int()) += cplx_type(-p1_flow, 0.);
            Sbus.coeffRef(bus_2_solver.cast_int()) += cplx_type(-p2_flow, 0.);
        }
    }
}

void HvdcLineContainer::compute_results(const Eigen::Ref<const RealVect> & Va,
                                        const Eigen::Ref<const RealVect> & Vm,
                                        const Eigen::Ref<const CplxVect> & V,
                                        const SolverBusIdVect & id_grid_to_solver,
                                        const RealVect & bus_vn_kv,
                                        real_type sn_mva,
                                        bool ac)
{
    side_1_.compute_results(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
    side_2_.compute_results(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);

    // droop lines: the active power is theta-dependent, recompute it from the
    // solved angles (the stations' "target" p only stores the fixed-setpoint
    // fallback values)
    const int nb_hvdc = static_cast<int>(nb());
    Eigen::Ref<RealVect> res_p_1 = get_res_p_side_1();
    Eigen::Ref<RealVect> res_p_2 = get_res_p_side_2();
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        if(!droop_enabled_[hvdc_id] || !status_global_[hvdc_id]) continue;
        const GlobalBusId bus_1 = get_bus_side_1(hvdc_id);
        const GlobalBusId bus_2 = get_bus_side_2(hvdc_id);
        const SolverBusId bus_1_solver = id_grid_to_solver[bus_1.cast_int()];
        const SolverBusId bus_2_solver = id_grid_to_solver[bus_2.cast_int()];
        if((bus_1_solver.cast_int() == _deactivated_bus_id) ||
           (bus_2_solver.cast_int() == _deactivated_bus_id)) continue;
        const real_type theta_1 = Va(bus_1_solver.cast_int());
        const real_type theta_2 = Va(bus_2_solver.cast_int());
        const real_type raw = p0_mw_(hvdc_id) + k_mw_per_rad_(hvdc_id) * (theta_1 - theta_2);
        real_type p1_flow, p2_flow;
        if(ac){
            droop_flows_mw(hvdc_id, raw, p1_flow, p2_flow);
        }else{
            if(status_droop_(hvdc_id) == 0){
                // DC linear regime: lossless, mirror of pyloadflow dc.py
                p1_flow = raw;
                p2_flow = -raw;
            }else{
                droop_flows_mw(hvdc_id, raw, p1_flow, p2_flow);
            }
        }
        // res_p is the station injection (generator convention) = -flow
        res_p_1(hvdc_id) = -p1_flow;
        res_p_2(hvdc_id) = -p2_flow;
    }
}

void HvdcLineContainer::save_binary(const std::string & path) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

HvdcLineContainer HvdcLineContainer::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<HvdcLineContainer>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

} // namespace ls2g
