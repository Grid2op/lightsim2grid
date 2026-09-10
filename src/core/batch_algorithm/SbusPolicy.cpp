// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SbusPolicy.hpp"

namespace ls2g {

namespace {
// the solver bus of every active element of `structure_data`, -1 for an inactive
// one; an active element on a disconnected bus is the error fill_SBus_* raised
template<class T>
std::vector<int> route(const T & structure_data,
                       const SolverBusIdVect & id_me_to_solver,
                       const char * algo_name,
                       const char * fun_name)
{
    const size_t nb_el = structure_data.nb();
    const auto & el_status = structure_data.get_status();
    const auto & el_bus_id = structure_data.get_bus_id();
    std::vector<int> res(nb_el, -1);
    for(size_t el_id = 0; el_id < nb_el; ++el_id){
        if(!el_status[el_id]) continue;
        const GlobalBusId bus_id_me = el_bus_id(el_id);
        int bus_id_solver = BaseConstants::_deactivated_bus_id;
        if(bus_id_me.cast_int() != BaseConstants::_deactivated_bus_id){
            bus_id_solver = id_me_to_solver[bus_id_me.cast_int()].cast_int();
        }
        if(bus_id_solver == BaseConstants::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << algo_name << "::" << fun_name << ": the element with id ";
            exc_ << el_id;
            exc_ << " is connected to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        res[el_id] = bus_id_solver;
    }
    return res;
}
}  // anonymous namespace

void SbusPolicy::Vary::prepare(const LSGrid & grid_model,
                               bool ac_solver_used,
                               int nb_buses_solver,
                               const SolverBusIdVect & id_me_to_solver,
                               const Eigen::Ref<const CplxVect> & complete_sbus_pu,
                               real_type sn_mva,
                               Eigen::Index nb_steps,
                               const char * algo_name)
{
    const auto & generators = grid_model.get_generators_as_data();
    const auto & s_generators = grid_model.get_static_generators_as_data();
    const auto & loads = grid_model.get_loads_as_data();

    gen_bus_ = route(generators, id_me_to_solver, algo_name, "fill_SBus_real");
    sgen_bus_ = route(s_generators, id_me_to_solver, algo_name, "fill_SBus_real");
    load_bus_ = route(loads, id_me_to_solver, algo_name, "fill_SBus_real");

    gen_target_p_ = grid_model.get_gen_target_p();
    sgen_target_p_ = grid_model.get_sgen_target_p();
    load_target_p_ = grid_model.get_load_target_p();
    load_target_q_ = loads.get_target_q();
    const int nb_gen = static_cast<int>(generators.nb());
    gen_target_q_ = RealVect::Zero(nb_gen);
    gen_vreg_.assign(static_cast<size_t>(nb_gen), 0);
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        gen_target_q_(gen_id) = generators.get_target_q_mvar(gen_id);
        gen_vreg_[static_cast<size_t>(gen_id)] = generators.get_voltage_regulator_on(gen_id) ? 1 : 0;
    }

    constant_pu_ = constant_sbus_pu(grid_model, complete_sbus_pu, nb_buses_solver, id_me_to_solver, algo_name);
    sn_mva_ = sn_mva;
    nb_buses_solver_ = nb_buses_solver;
    nb_steps_ = nb_steps;
    sbuses = CplxMat();   // a previous compute()'s matrix, if anyone asked for one
    (void)ac_solver_used;  // kept for signature symmetry with the caller's other prep calls
}

void SbusPolicy::Vary::fill_row(Eigen::Index i, CplxVect & row) const
{
    // The whole-matrix build this replaces added each element's column into the
    // bus's column: for one row, that is one complex addition per element, in
    // element order, container after container. Same operands, same order, same
    // operations (a real value enters as the complex (p, 0), a reactive one as
    // my_i * (q, 0), exactly as fill_SBus_real / fill_SBus_imag formed them).
    row.resize(nb_buses_solver_);
    row.setZero();

    const bool own_gen_p = gen_p.rows() > 0;
    for(size_t g = 0; g < gen_bus_.size(); ++g){
        if(gen_bus_[g] < 0) continue;
        const real_type p = own_gen_p ? gen_p(i, static_cast<Eigen::Index>(g)) : gen_target_p_(static_cast<Eigen::Index>(g));
        row(gen_bus_[g]) += cplx_type(p, 0.);
    }
    const bool own_sgen_p = sgen_p.rows() > 0;
    for(size_t g = 0; g < sgen_bus_.size(); ++g){
        if(sgen_bus_[g] < 0) continue;
        const real_type p = own_sgen_p ? sgen_p(i, static_cast<Eigen::Index>(g)) : sgen_target_p_(static_cast<Eigen::Index>(g));
        row(sgen_bus_[g]) += cplx_type(p, 0.);
    }
    const bool own_load_p = load_p.rows() > 0;
    for(size_t l = 0; l < load_bus_.size(); ++l){
        if(load_bus_[l] < 0) continue;
        const real_type p = own_load_p ? load_p(i, static_cast<Eigen::Index>(l)) : load_target_p_(static_cast<Eigen::Index>(l));
        row(load_bus_[l]) -= cplx_type(p, 0.);
    }
    const bool own_load_q = load_q.rows() > 0;
    for(size_t l = 0; l < load_bus_.size(); ++l){
        if(load_bus_[l] < 0) continue;
        const real_type q = own_load_q ? load_q(i, static_cast<Eigen::Index>(l)) : load_target_q_(static_cast<Eigen::Index>(l));
        row(load_bus_[l]) -= BaseConstants::my_i * cplx_type(q, 0.);
    }

    // generator contingencies: take back, for this row, exactly what the gen pass
    // above put in for a generator this row disconnects -- and, for one that does
    // not regulate voltage, the reactive setpoint that reached the row through
    // constant_pu_. A generator this row keeps subtracted a zero in the whole-matrix
    // build, which leaves every value as it is: skipped here.
    if(gen_off.rows() > 0 && i < gen_off.rows()){
        const Eigen::Index nb_cols = gen_off.cols();
        for(size_t g = 0; g < gen_bus_.size(); ++g){
            const Eigen::Index gen_id = static_cast<Eigen::Index>(g);
            if(gen_bus_[g] < 0 || gen_id >= nb_cols || !gen_off(i, gen_id)) continue;
            const real_type p = own_gen_p ? gen_p(i, gen_id) : gen_target_p_(gen_id);
            row(gen_bus_[g]) -= cplx_type(p, 0.);
        }
        for(size_t g = 0; g < gen_bus_.size(); ++g){
            const Eigen::Index gen_id = static_cast<Eigen::Index>(g);
            if(gen_bus_[g] < 0 || gen_id >= nb_cols || !gen_off(i, gen_id)) continue;
            if(gen_vreg_[g]) continue;
            row(gen_bus_[g]) -= BaseConstants::my_i * cplx_type(gen_target_q_(gen_id), 0.);
        }
    }

    // the same Eigen expression the whole-matrix build applied (a complex division,
    // vectorized): written this way rather than as a scalar division so the
    // rounding is the one it always was
    if(abs(sn_mva_ - 1.0) > BaseConstants::_tol_equal_float) row.array() /= static_cast<cplx_type>(sn_mva_);

    // ... and then everything else the gridmodel stamps into Sbus that the four
    // matrices above do not cover (storage/SVC/HVDC/..., see constant_sbus_pu), same
    // for every step
    row += constant_pu_;
}

const SbusPolicy::Vary::CplxMat & SbusPolicy::Vary::materialize() const
{
    if(sbuses.rows() == nb_steps_ && sbuses.cols() == nb_buses_solver_ && nb_steps_ > 0) return sbuses;
    sbuses = CplxMat(nb_steps_, nb_buses_solver_);
    CplxVect row;
    for(Eigen::Index i = 0; i < nb_steps_; ++i){
        fill_row(i, row);
        sbuses.row(i) = row.transpose();
    }
    return sbuses;
}

CplxVect SbusPolicy::Vary::constant_sbus_pu(const LSGrid & grid_model,
                                            const Eigen::Ref<const CplxVect> & complete_sbus_pu,
                                            int nb_buses_solver,
                                            const SolverBusIdVect & id_me_to_solver,
                                            const char * algo_name) const
{
    // the complete per-unit injection, straight from the gridmodel (the caller's own
    // ac_cache_.inj (ac) / dc_cache_.inj.cast<cplx_type>() (dc))
    CplxVect res = complete_sbus_pu;
    if(static_cast<int>(res.size()) != nb_buses_solver){
        std::ostringstream exc_;
        exc_ << algo_name << "::constant_sbus_pu: the gridmodel injection has "
             << res.size() << " entries while the solver has " << nb_buses_solver
             << " buses. prepare_solver_input_base() must run first.";
        throw std::runtime_error(exc_.str());
    }

    // ... minus the share the per-step matrices rebuild, evaluated at the gridmodel's
    // own target values, so that the two cancel and only the rest is left
    const auto & generators = grid_model.get_generators_as_data();
    const auto & s_generators = grid_model.get_static_generators_as_data();
    const auto & loads = grid_model.get_loads_as_data();

    // The gridmodel's target vectors, presented to fill_SBus_* as the single-step
    // matrices it expects, WITHOUT copying them: a 1 x n row-major matrix has exactly
    // the memory layout of a contiguous n-vector, so an Eigen::Map over the vector's
    // own storage is a view. The Refs are held in named locals rather than being
    // called inline: that way the Map never depends on a temporary Ref's lifetime.
    const Eigen::Ref<const RealVect> gen_p_v = grid_model.get_gen_target_p();
    const Eigen::Ref<const RealVect> sgen_p_v = grid_model.get_sgen_target_p();
    const Eigen::Ref<const RealVect> load_p_v = grid_model.get_load_target_p();
    const Eigen::Ref<const RealVect> load_q_v = loads.get_target_q();
    const Eigen::Map<const RealMat> gen_p(gen_p_v.data(), 1, gen_p_v.size());
    const Eigen::Map<const RealMat> sgen_p(sgen_p_v.data(), 1, sgen_p_v.size());
    const Eigen::Map<const RealMat> load_p(load_p_v.data(), 1, load_p_v.size());
    const Eigen::Map<const RealMat> load_q(load_q_v.data(), 1, load_q_v.size());

    // `accounted` does have to be its own buffer: it is written by fill_SBus_* and
    // then scaled by sn_mva, which must NOT touch `res` (already per-unit). One
    // vector-sized allocation per call, not per step.
    CplxMat accounted = CplxMat::Zero(1, nb_buses_solver);
    bool add_ = true;
    fill_SBus_real(accounted, generators, gen_p, id_me_to_solver, add_, algo_name);
    fill_SBus_real(accounted, s_generators, sgen_p, id_me_to_solver, add_, algo_name);
    add_ = false;
    fill_SBus_real(accounted, loads, load_p, id_me_to_solver, add_, algo_name);
    fill_SBus_imag(accounted, loads, load_q, id_me_to_solver, add_, algo_name);

    const real_type sn_mva = grid_model.get_sn_mva();
    if(abs(sn_mva - 1.0) > BaseConstants::_tol_equal_float) accounted.array() /= static_cast<cplx_type>(sn_mva);
    res -= accounted.row(0).transpose();
    return res;
}

} // namespace ls2g
