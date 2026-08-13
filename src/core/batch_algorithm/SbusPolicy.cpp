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
// a nb_steps x target.size() matrix with every row equal to `target` -- the default
// injection for an axis the caller never set via modify_*.
SbusPolicy::Vary::RealMat broadcast_default(const Eigen::Ref<const RealVect> & target, Eigen::Index nb_steps)
{
    SbusPolicy::Vary::RealMat res(nb_steps, target.size());
    res.rowwise() = target.transpose();
    return res;
}
}  // anonymous namespace

void SbusPolicy::Vary::assemble(const LSGrid & grid_model,
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

    const RealMat gen_p_use = gen_p.rows() > 0 ? gen_p : broadcast_default(grid_model.get_gen_target_p(), nb_steps);
    const RealMat sgen_p_use = sgen_p.rows() > 0 ? sgen_p : broadcast_default(grid_model.get_sgen_target_p(), nb_steps);
    const RealMat load_p_use = load_p.rows() > 0 ? load_p : broadcast_default(grid_model.get_load_target_p(), nb_steps);
    const RealMat load_q_use = load_q.rows() > 0 ? load_q : broadcast_default(loads.get_target_q(), nb_steps);

    sbuses = CplxMat::Zero(nb_steps, nb_buses_solver);
    bool add_ = true;
    fill_SBus_real(sbuses, generators, gen_p_use, id_me_to_solver, add_, algo_name);
    fill_SBus_real(sbuses, s_generators, sgen_p_use, id_me_to_solver, add_, algo_name);
    add_ = false;
    fill_SBus_real(sbuses, loads, load_p_use, id_me_to_solver, add_, algo_name);
    fill_SBus_imag(sbuses, loads, load_q_use, id_me_to_solver, add_, algo_name);

    if(abs(sn_mva - 1.0) > BaseConstants::_tol_equal_float) sbuses.array() /= static_cast<cplx_type>(sn_mva);

    // ... and then everything else the gridmodel stamps into Sbus that the four
    // matrices above do not cover (storage/SVC/HVDC/..., see constant_sbus_pu), same
    // for every step
    sbuses.rowwise() += constant_sbus_pu(grid_model, complete_sbus_pu, nb_buses_solver, id_me_to_solver, algo_name).transpose();
    (void)ac_solver_used;  // kept for signature symmetry with the caller's other prep calls
}

CplxVect SbusPolicy::Vary::constant_sbus_pu(const LSGrid & grid_model,
                                            const Eigen::Ref<const CplxVect> & complete_sbus_pu,
                                            int nb_buses_solver,
                                            const SolverBusIdVect & id_me_to_solver,
                                            const char * algo_name) const
{
    // the complete per-unit injection, straight from the gridmodel (the caller's own
    // Sbus_ (ac) / Pbus_.cast<cplx_type>() (dc))
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
