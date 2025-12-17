// Copyright (c) 2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DCLineContainer.hpp"

#include <iostream>
#include <sstream>

DCLineInfo::DCLineInfo(const DCLineContainer & r_data_dcline, int my_id):
    TwoSidesContainer<GeneratorContainer>::TwoSidesInfo(r_data_dcline, my_id),
    target_p_1_mw(0.),
    p_2_mw(0.),
    target_vm_1_pu(0.),
    target_vm_2_pu(0.),
    loss_pct(0.),
    loss_mw(0.),
    gen_side_1(r_data_dcline.side_1_, my_id),
    gen_side_2(r_data_dcline.side_2_, my_id)
{
    if (my_id < 0) return;
    if (my_id >= r_data_dcline.nb()) return;
    loss_pct = r_data_dcline.loss_percent_(my_id);
    loss_mw = r_data_dcline.loss_mw_(my_id);

    target_p_1_mw = gen_side_1.target_p_mw;
    p_2_mw = gen_side_2.target_p_mw;

    target_vm_1_pu = gen_side_1.target_vm_pu;
    target_vm_2_pu = gen_side_2.target_vm_pu;
}

DCLineContainer::StateRes DCLineContainer::get_state() const
{
    std::vector<real_type> loss_percent(loss_percent_.begin(), loss_percent_.end());
    std::vector<real_type> loss_mw(loss_mw_.begin(), loss_mw_.end());
    DCLineContainer::StateRes res(get_tsc_state(),
                                  loss_percent,
                                  loss_mw);
    return res;
}

void DCLineContainer::set_state(DCLineContainer::StateRes & my_state){
    set_tsc_state(std::get<0>(my_state));
    std::vector<real_type> & loss_percent = std::get<1>(my_state);
    std::vector<real_type> & loss_mw = std::get<2>(my_state);
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
    init_tsc(branch_from_id, branch_to_id, "DC lines");

    side_1_.init(p_mw, vm_or_pu, min_q_or, max_q_or, branch_from_id);
    RealVect p_ex = p_mw;
    Eigen::Index size_ = p_mw.size();
    for(Eigen::Index i = 0; i < size_; ++i){
        p_ex(i) = get_to_mw(i, p_ex(i));
    }
    side_2_.init(p_ex, vm_ex_pu, min_q_ex, max_q_ex, branch_to_id);
}
