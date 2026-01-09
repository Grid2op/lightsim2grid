// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "LineContainer.hpp"

#include <sstream>

void LineContainer::init(
    const RealVect & branch_r,
    const RealVect & branch_x,
    const CplxVect & branch_h,
    const Eigen::VectorXi & branch_from_id,
    const Eigen::VectorXi & branch_to_id
)
{
    /**
    This method initialize the Ybus matrix from the branch matrix.
    It has to be called once when the solver is initialized. Afterwards, a call to
    "updateYbus" should be made for performance optimiaztion instead. //TODO
    **/

   init(
    branch_r, 
    branch_x,
    0.5 * branch_h,
    0.5 * branch_h,
    branch_from_id,
    branch_to_id);
}

void LineContainer::init(const RealVect & branch_r,
                         const RealVect & branch_x,
                         const CplxVect & branch_h_or,
                         const CplxVect & branch_h_ex,
                         const Eigen::VectorXi & branch_from_id,
                         const Eigen::VectorXi & branch_to_id
                         )
{
    /**
    This method initialize the Ybus matrix from the branch matrix.
    It has to be called once when the solver is initialized. Afterwards, a call to
    "updateYbus" should be made for performance optimiaztion instead. //TODO
    **/

    // TODO check what can be checked: branch_* have same size, no id in branch_to_id that are
    // TODO not in [0, .., buv_vn_kv.size()] etc.

    //TODO consistency with trafo: have a converter methods to convert this value into pu, and store the pu
    // in this method

    int size = static_cast<int>(branch_r.size());
    GenericContainer::check_size(branch_r, size, "branch_r");
    GenericContainer::check_size(branch_x, size, "branch_x");
    GenericContainer::check_size(branch_h_or, size, "branch_h_or");
    GenericContainer::check_size(branch_h_ex, size, "branch_h_ex");
    GenericContainer::check_size(branch_from_id, size, "branch_from_id");
    GenericContainer::check_size(branch_to_id, size, "branch_to_id");

    h_side_1_ = branch_h_or;
    h_side_2_ = branch_h_ex;
    r_ = branch_r;
    x_ = branch_x;
    init_tsc(branch_from_id, branch_to_id, "trafo");
    _update_model_coeffs();
    reset_results();
}
