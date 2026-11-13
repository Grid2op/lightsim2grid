// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "NRSystem.hpp"

// out-of-line on purpose: NRSystem.hpp only forward-declares LSGrid (it is
// included BY LSGrid.hpp), the full type is needed to pull the droop data
#include "LSGrid.hpp"

#include <algorithm>

namespace ls2g {

void Hvdc::update_state(
    const Base                       * /*nr_system_base_ptr*/,
    const LSGrid                     * lsgrid_ptr,
    const EigenRefConstCplxSpMat     & /*Ybus*/,
    const Eigen::Ref<const CplxVect> & /*Sbus*/,
    const Eigen::Ref<const RealVect> & /*slack_weights*/
)
{
    data_.clear();
    if(lsgrid_ptr != nullptr) lsgrid_ptr->fill_hvdc_droop_solver_data(data_, true);
    my_size_ = data_.size();
    off_.assign(static_cast<std::size_t>(my_size_), 0);
    for(int k = 0; k < my_size_; ++k){
        if(std::find(tripped_lines_.begin(), tripped_lines_.end(), data_.hvdc_id[static_cast<std::size_t>(k)]) != tripped_lines_.end())
            off_[static_cast<std::size_t>(k)] = 1;
    }
}

} // namespace ls2g
