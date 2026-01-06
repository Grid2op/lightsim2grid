// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SGenContainer.hpp"

#include <iostream>

void SGenContainer::init(const RealVect & sgen_p,
                         const RealVect & sgen_q,
                         const RealVect & sgen_pmin,
                         const RealVect & sgen_pmax,
                         const RealVect & sgen_qmin,
                         const RealVect & sgen_qmax,
                         const Eigen::VectorXi & sgen_bus_id)
{
    init_osc_pq(sgen_p, sgen_q, sgen_bus_id, "static_generators");
    
    int size = nb();
    check_size(sgen_pmin, size, "sgen_pmin");
    check_size(sgen_pmax, size, "sgen_pmax");
    check_size(sgen_qmin, size, "sgen_qmin");
    check_size(sgen_qmax, size, "sgen_qmax");

    p_min_mw_ = sgen_pmin;
    p_max_mw_ = sgen_pmax;
    q_min_mvar_ = sgen_qmin;
    q_max_mvar_ = sgen_qmax;
    reset_results();
}

SGenContainer::StateRes SGenContainer::get_state() const
{
     std::vector<real_type> p_min(p_min_mw_.begin(), p_min_mw_.end());
     std::vector<real_type> p_max(p_max_mw_.begin(), p_max_mw_.end());
     std::vector<real_type> q_min(q_min_mvar_.begin(), q_min_mvar_.end());
     std::vector<real_type> q_max(q_max_mvar_.begin(), q_max_mvar_.end());
     SGenContainer::StateRes res(get_osc_pq_state(), p_min, p_max, q_min, q_max);
     return res;
}

void SGenContainer::set_state(SGenContainer::StateRes & my_state )
{    
    set_osc_pq_state(std::get<0>(my_state));

    std::vector<real_type> & p_min = std::get<1>(my_state);
    std::vector<real_type> & p_max = std::get<2>(my_state);
    std::vector<real_type> & q_min = std::get<3>(my_state);
    std::vector<real_type> & q_max = std::get<4>(my_state);
    const auto size = nb();

    GenericContainer::check_size(p_min, size, "p_min");
    GenericContainer::check_size(p_max, size, "p_max");
    GenericContainer::check_size(q_min, size, "q_min");
    GenericContainer::check_size(q_max, size, "q_max");

    p_min_mw_ = RealVect::Map(&p_min[0], size);
    p_max_mw_ = RealVect::Map(&p_max[0], size);
    q_min_mvar_ = RealVect::Map(&q_min[0], size);
    q_max_mvar_ = RealVect::Map(&q_max[0], size);
    reset_results();
}

void SGenContainer::fillSbus(CplxVect & Sbus, const std::vector<SolverBusId> & id_grid_to_solver, bool ac) const {
    const int nb_sgen = nb();
    GlobalBusId bus_id_me;
    SolverBusId bus_id_solver;
    cplx_type tmp;
    for(int sgen_id = 0; sgen_id < nb_sgen; ++sgen_id){
        //  i don't do anything if the static generator is disconnected
        if(!status_[sgen_id]) continue;

        bus_id_me = bus_id_(sgen_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "SGenContainer::fillSbus: Static Generator with id ";
            exc_ << sgen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "SGenContainer::fillSbus: Static Generator with id ";
            exc_ << sgen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        tmp = {target_p_mw_(sgen_id), target_q_mvar_(sgen_id)};
        Sbus.coeffRef(bus_id_solver.cast_int()) += tmp;
    }
}
