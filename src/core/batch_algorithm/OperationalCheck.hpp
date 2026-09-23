// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LS2G_OPERATIONAL_CHECK_H
#define LS2G_OPERATIONAL_CHECK_H

// The OPERATIONAL limit checks of a solved voltage vector: bus voltages against their
// [vmin, vmax], branch currents (both sides) against their thermal limits. Shared by
// the batch algorithms (`compute_limit_violations`, one call per row and one for the
// base case) and by LSGrid::get_violations (the grid's own last solve). Moved here
// verbatim from BaseBatchSweep.hpp so that LSGrid.cpp can use them without the batch
// template; the namespace is kept for the existing callers.

// isnan / isfinite / sqrt: the C header, not <cmath>. These are used unqualified below
// (isnan(...), not std::isnan(...)), and only math.h guarantees them in the global
// namespace across compilers -- cmath only guarantees std::isnan / std::isfinite.
#include <math.h>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "BaseConstants.hpp"
#include "SubstationContainer.hpp"
#include "TaggedIdVec.hpp"
#include "Utils.hpp"
#include "LimitViolation.hpp"

namespace ls2g {

namespace batch_sweep_detail {

// appends BUS limit violations for a single, already-solved voltage vector `V` (solver
// numbering). `masked_solver_ids` is nullptr for the base ("n") case; in "handle
// disconnected grid" mode it is this contingency's masked solver-bus-id list -- those
// buses are forced to V=0 post-solve and must never be checked. `subs` resolves the
// violating bus's substation name. `threshold` (in ]0., 1.]) tightens both checks; 1.
// is the "report exactly at the configured limit" default. Moved verbatim (as a free
// function, `inline` so it can live in this header without an ODR violation) from the
// pre-refactor ContingencyAnalysis.cpp's anonymous namespace.
inline void check_bus_voltage_violations(
    const Eigen::Ref<const CplxVect> & V,
    const SolverBusIdVect & id_me_to_solver,
    const Eigen::Ref<const RealVect> & bus_vmin_kv,
    const Eigen::Ref<const RealVect> & bus_vmax_kv,
    const Eigen::Ref<const RealVect> & bus_vn_kv,
    const SubstationContainer & subs,
    real_type threshold,
    const std::vector<int> * masked_solver_ids,
    std::vector<LimitViolation> & out)
{
    const Eigen::Index nb_bus = bus_vmin_kv.size();
    if(nb_bus == 0) return;  // voltage limits never configured on this grid

    const std::vector<std::string> & sub_names = subs.get_sub_names();  // empty if never set

    std::vector<bool> is_masked;
    if(masked_solver_ids != nullptr && !masked_solver_ids->empty()){
        is_masked.assign(static_cast<size_t>(V.size()), false);
        for(int b : *masked_solver_ids){
            if(b >= 0 && static_cast<size_t>(b) < is_masked.size()) is_masked[static_cast<size_t>(b)] = true;
        }
    }

    for(Eigen::Index grid_id = 0; grid_id < nb_bus; ++grid_id){
        const real_type vmin = bus_vmin_kv(grid_id);
        const real_type vmax = bus_vmax_kv(grid_id);
        if(isnan(vmin) && isnan(vmax)) continue;  // no limit configured for this specific bus

        const int solver_id = id_me_to_solver(static_cast<int>(grid_id)).cast_int();
        if(solver_id == BaseConstants::_deactivated_bus_id) continue;
        if(!is_masked.empty() && is_masked[static_cast<size_t>(solver_id)]) continue;

        const real_type vm_kv = std::abs(V(solver_id)) * bus_vn_kv(grid_id);
        const bool has_min = !isnan(vmin);
        const bool has_max = !isnan(vmax);

        const real_type vnom = bus_vn_kv(grid_id);
        real_type anchor = vnom;
        if(has_min && (anchor < vmin)) anchor = vmin;
        if(has_max && (anchor > vmax)) anchor = vmax;
        const real_type shrink = 1. - threshold;
        const real_type low_eff = vmin + shrink * (anchor - vmin);
        const real_type high_eff = vmax - shrink * (vmax - anchor);

        if(has_min && has_max && low_eff > high_eff){
            std::ostringstream exc_;
            exc_ << "ContingencyAnalysis: bus " << grid_id << " has inconsistent voltage "
                    "limits: its effective minimum (" << low_eff << " kV, from vmin = "
                 << vmin << " kV) is above its effective maximum (" << high_eff
                 << " kV, from vmax = " << vmax << " kV) at violation_threshold = "
                 << threshold << ". Check the vmin_kv / vmax_kv passed to "
                    "LSGrid::set_bus_voltage_limits (they must satisfy vmin <= vmax, and "
                    "should bracket the bus nominal voltage of " << vnom << " kV).";
            throw std::runtime_error(exc_.str());
        }

        if(has_min && vm_kv <= low_eff){
            const int sub_id = subs.sub_id_of_bus(static_cast<int>(grid_id));
            const std::string sub_name = static_cast<size_t>(sub_id) < sub_names.size() ? sub_names[sub_id] : std::string();
            out.push_back(LimitViolation{ViolationElementType::BUS, static_cast<int>(grid_id), 0,
                                          LimitViolationType::LOW_VOLTAGE, vm_kv, vmin, sub_name});
        } else if(has_max && vm_kv >= high_eff){
            const int sub_id = subs.sub_id_of_bus(static_cast<int>(grid_id));
            const std::string sub_name = static_cast<size_t>(sub_id) < sub_names.size() ? sub_names[sub_id] : std::string();
            out.push_back(LimitViolation{ViolationElementType::BUS, static_cast<int>(grid_id), 0,
                                          LimitViolationType::HIGH_VOLTAGE, vm_kv, vmax, sub_name});
        }
    }
}

// appends CURRENT limit violations (both sides) for a single, already-solved voltage
// vector `V` (solver numbering), for every element of `structure_data` (LineContainer
// or TrafoContainer). `skip_ids` are the LOCAL (own-type) element ids disconnected BY
// THIS CONTINGENCY. Moved verbatim from the pre-refactor
// ContingencyAnalysis.cpp's anonymous namespace.
template<class T>
inline void check_current_violations(
    const T & structure_data,
    ViolationElementType el_type,
    const Eigen::Ref<const CplxVect> & V,
    const SolverBusIdVect & id_me_to_solver,
    const Eigen::Ref<const RealVect> & bus_vn_kv,
    bool ac_solver_used,
    real_type sn_mva,
    const Eigen::Ref<const RealVect> & limit1,
    const Eigen::Ref<const RealVect> & limit2,
    real_type threshold,
    const std::vector<int> & skip_ids,
    std::vector<LimitViolation> & out)
{
    if(limit1.size() == 0 && limit2.size() == 0) return;  // thermal limits never configured

    const auto & el_status = structure_data.get_status_global();
    const auto & status1 = structure_data.get_status_side_1();
    const auto & status2 = structure_data.get_status_side_2();
    const GlobalBusIdVect & bus_from = structure_data.get_bus_id_side_1();
    const GlobalBusIdVect & bus_to = structure_data.get_bus_id_side_2();
    const std::vector<std::string> & el_names = structure_data.get_names();  // empty if never set

    Eigen::Ref<const CplxVect> yac_eff_11 = structure_data.yac_eff_11();
    Eigen::Ref<const CplxVect> yac_eff_12 = structure_data.yac_eff_12();
    Eigen::Ref<const CplxVect> yac_eff_21 = structure_data.yac_eff_21();
    Eigen::Ref<const CplxVect> yac_eff_22 = structure_data.yac_eff_22();
    Eigen::Ref<const RealVect> ydc_11 = structure_data.ydc_11();
    Eigen::Ref<const RealVect> ydc_12 = structure_data.ydc_12();
    Eigen::Ref<const RealVect> ydc_21 = structure_data.ydc_21();
    Eigen::Ref<const RealVect> ydc_22 = structure_data.ydc_22();
    Eigen::Ref<const RealVect> dc_x_tau_shift = structure_data.dc_x_tau_shift();
    const bool has_tau_shift = dc_x_tau_shift.size() > 0;  // trafo only, empty for lines

    const size_t nb_el = structure_data.nb();
    const real_type sqrt_3 = sqrt(3.);

    std::vector<bool> skip;
    if(!skip_ids.empty()){
        skip.assign(nb_el, false);
        for(int id : skip_ids){
            if(id >= 0 && static_cast<size_t>(id) < nb_el) skip[static_cast<size_t>(id)] = true;
        }
    }

    for(size_t el_id = 0; el_id < nb_el; ++el_id){
        if(!el_status[el_id]) continue;
        if(!skip.empty() && skip[el_id]) continue;

        const Eigen::Index el_idx = static_cast<Eigen::Index>(el_id);
        const bool has_lim1 = limit1.size() > 0 && !isnan(limit1(el_idx));
        const bool has_lim2 = limit2.size() > 0 && !isnan(limit2(el_idx));
        if(!has_lim1 && !has_lim2) continue;

        const bool s1 = status1[el_id];
        const bool s2 = status2[el_id];
        int bus_from_me = BaseConstants::_deactivated_bus_id;
        int bus_to_me = BaseConstants::_deactivated_bus_id;
        int solver_from = BaseConstants::_deactivated_bus_id;
        int solver_to = BaseConstants::_deactivated_bus_id;
        if(s1){
            bus_from_me = bus_from(static_cast<int>(el_id)).cast_int();
            solver_from = id_me_to_solver(bus_from_me).cast_int();
            if(solver_from == BaseConstants::_deactivated_bus_id) continue;
        }
        if(s2){
            bus_to_me = bus_to(static_cast<int>(el_id)).cast_int();
            solver_to = id_me_to_solver(bus_to_me).cast_int();
            if(solver_to == BaseConstants::_deactivated_bus_id) continue;
        }

        const cplx_type Efrom = s1 ? V(solver_from) : cplx_type(0., 0.);
        const cplx_type Eto = s2 ? V(solver_to) : cplx_type(0., 0.);
        const real_type v_from_kv = s1 ? std::abs(Efrom) * bus_vn_kv(bus_from_me) : real_type(1.);
        const real_type v_to_kv = s2 ? std::abs(Eto) * bus_vn_kv(bus_to_me) : real_type(1.);

        real_type amps1 = 0.;
        real_type amps2 = 0.;
        if(ac_solver_used){
            const cplx_type I_from = std::conj(yac_eff_11(el_idx) * Efrom + yac_eff_12(el_idx) * Eto);
            const cplx_type S_from = Efrom * I_from;
            amps1 = std::abs(S_from) * sn_mva / (sqrt_3 * v_from_kv);

            const cplx_type I_to = std::conj(yac_eff_22(el_idx) * Eto + yac_eff_21(el_idx) * Efrom);
            const cplx_type S_to = Eto * I_to;
            amps2 = std::abs(S_to) * sn_mva / (sqrt_3 * v_to_kv);
        } else if(s1 && s2){
            const real_type theta_from = std::arg(Efrom);
            const real_type theta_to = std::arg(Eto);
            real_type p_from = (ydc_11(el_idx) * theta_from + ydc_12(el_idx) * theta_to) * sn_mva;
            if(has_tau_shift) p_from -= dc_x_tau_shift(el_idx);
            amps1 = std::abs(p_from) / (sqrt_3 * v_from_kv);

            real_type p_to = (ydc_22(el_idx) * theta_to + ydc_21(el_idx) * theta_from) * sn_mva;
            if(has_tau_shift) p_to += dc_x_tau_shift(el_idx);
            amps2 = std::abs(p_to) / (sqrt_3 * v_to_kv);
        }

        // the name is copied only into a violation: one std::string per element
        // examined was a malloc per branch per row on a grid with long names
        if(has_lim1 && amps1 >= threshold * limit1(el_idx)){
            out.push_back(LimitViolation{el_type, static_cast<int>(el_id), 1,
                                          LimitViolationType::CURRENT, amps1, limit1(el_idx),
                                          el_id < el_names.size() ? el_names[el_id] : std::string()});
        }
        if(has_lim2 && amps2 >= threshold * limit2(el_idx)){
            out.push_back(LimitViolation{el_type, static_cast<int>(el_id), 2,
                                          LimitViolationType::CURRENT, amps2, limit2(el_idx),
                                          el_id < el_names.size() ? el_names[el_id] : std::string()});
        }
    }
}

}  // namespace batch_sweep_detail

}  // namespace ls2g

#endif  // LS2G_OPERATIONAL_CHECK_H
