// Copyright (c) 2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DCLineContainer.h"

#include <iostream>
#include <sstream>

DCLineContainer::StateRes DCLineContainer::get_state() const
{
    std::vector<real_type> loss_percent(loss_percent_.begin(), loss_percent_.end());
    std::vector<real_type> loss_mw(loss_mw_.begin(), loss_mw_.end());
    std::vector<bool> status = status_;
    DCLineContainer::StateRes res(names_,
                                  from_gen_.get_state(),
                                  to_gen_.get_state(),
                                  loss_percent,
                                  loss_mw,
                                  status);
    return res;
}

void DCLineContainer::set_state(DCLineContainer::StateRes & my_state){
    reset_results();
    names_ = std::get<0>(my_state);
    from_gen_.set_state(std::get<1>(my_state));
    to_gen_.set_state(std::get<2>(my_state));
    std::vector<real_type> & loss_percent = std::get<3>(my_state);
    std::vector<real_type> & loss_mw = std::get<4>(my_state);
    std::vector<bool> & status = std::get<5>(my_state);
    status_ = status;
    loss_percent_ = RealVect::Map(&loss_percent[0], loss_percent.size());
    loss_mw_ = RealVect::Map(&loss_mw[0], loss_percent.size());
    reset_results();
}

void DCLineContainer::init(const Eigen::VectorXi & branch_from_id,
                           const Eigen::VectorXi & branch_to_id,
                           const RealVect & p_mw,
                           const RealVect & loss_percent,
                           const RealVect & loss_mw,
                           const RealVect & vm_or_pu,
                           const RealVect & vm_ex_pu,
                           const RealVect & min_q_or,
                           const RealVect & max_q_or,
                           const RealVect & min_q_ex,
                           const RealVect & max_q_ex){
    loss_percent_ = loss_percent;
    loss_mw_ = loss_mw;
    status_ = std::vector<bool>(branch_from_id.size(), true);

    from_gen_.init(p_mw, vm_or_pu, min_q_or, max_q_or, branch_from_id);
    RealVect p_ex = p_mw;
    Eigen::Index size_ = p_mw.size();
    for(Eigen::Index i = 0; i < size_; ++i){
        p_ex(i) = get_to_mw(i, p_ex(i));
    }
    to_gen_.init(p_ex, vm_ex_pu, min_q_ex, max_q_ex, branch_to_id);
}

void DCLineContainer::nb_line_end(std::vector<int> & res) const
{
    const Eigen::Index nb = from_gen_.nb();
    const auto & bus_or_id = get_bus_id_or();
    const auto & bus_ex_id = get_bus_id_ex();
    for(Eigen::Index i = 0; i < nb; ++i){
        if(!status_[i]) continue;
        auto bus_or = bus_or_id(i);
        auto bus_ex = bus_ex_id(i);
        res[bus_or] += 1;
        res[bus_ex] += 1;
    }
}

// TODO DC LINE: one side might be in the connected comp and not the other !
void DCLineContainer::disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component)
{
    const Eigen::Index nb = from_gen_.nb();
    const auto & bus_or_id = get_bus_id_or();
    const auto & bus_ex_id = get_bus_id_ex(); 
    SolverControl unused_solver_control;
    for(Eigen::Index i = 0; i < nb; ++i){
        if(!status_[i]) continue;
        auto bus_or = bus_or_id(i);
        auto bus_ex = bus_ex_id(i);
        if(!busbar_in_main_component[bus_or]) {
            from_gen_.deactivate(i, unused_solver_control);
        }
        if(!busbar_in_main_component[bus_ex]) {
            to_gen_.deactivate(i, unused_solver_control);
        }
        // if(!busbar_in_main_component[bus_or] || !busbar_in_main_component[bus_ex]){
        //     bool tmp = false;
        //     deactivate(i, tmp);
        // }
    }
}
