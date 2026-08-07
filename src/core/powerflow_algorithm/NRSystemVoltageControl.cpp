// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "NRSystem.hpp"

// out-of-line on purpose: NRSystem.hpp only forward-declares LSGrid (it is
// included BY LSGrid.hpp), the full type is needed to pull the controller data
#include "LSGrid.hpp"

namespace ls2g {

void VoltageControl::update_state(
    const Base                       * /*nr_system_base_ptr*/,
    const LSGrid                     * lsgrid_ptr,
    const EigenRefConstCplxSpMat     & /*Ybus*/,
    const Eigen::Ref<const CplxVect> & /*Sbus*/,
    const Eigen::Ref<const RealVect> & /*slack_weights*/,
    const AlgoControl                & solver_control
)
{
    // The controller data is expensive to rebuild (fill_voltage_control_solver_data
    // groups the controllers in a pass quadratic in their number), and can only
    // change between solves -- the grid is const during one. Skip it when
    // AlgoControl reports nothing relevant changed.
    if(nr_extension_data_is_stale(solver_control) || !data_cached_){
        data_.clear();
        if(lsgrid_ptr != nullptr) lsgrid_ptr->fill_voltage_control_solver_data(data_, true);
        my_size_ = data_.n_controllers();
        data_cached_ = true;
    }
    // NOT cached: the reactive injection is solver STATE, not grid data, and its
    // per-solve initial guess is 0 (generator convention) whatever happened to
    // the grid. setZero() rather than a fresh vector so an unchanged controller
    // count does not reallocate.
    if(q_.size() != my_size_) q_ = RealVect::Zero(my_size_);
    else q_.setZero();
}

} // namespace ls2g
