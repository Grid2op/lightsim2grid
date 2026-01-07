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
    // init_tsc_rxha();
}

// void LineContainer::_update_model_coeffs()
// {
//     const auto my_size = r_.size();

//     yac_11_ = CplxVect::Zero(my_size);
//     yac_12_ = CplxVect::Zero(my_size);
//     yac_21_ = CplxVect::Zero(my_size);
//     yac_22_ = CplxVect::Zero(my_size);

//     ydc_11_ = CplxVect::Zero(my_size);
//     ydc_12_ = CplxVect::Zero(my_size);
//     ydc_21_ = CplxVect::Zero(my_size);
//     ydc_22_ = CplxVect::Zero(my_size);
//     for(int i = 0; i < my_size; ++i)
//     {
//         // for AC
//         // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.2
//         const cplx_type ys = 1. / cplx_type(r_(i), x_(i));
//         const cplx_type h_or = h_side_1_(i);
//         const cplx_type h_ex = h_side_2_(i);
//         yac_11_(i) = (ys + h_or);
//         yac_22_(i) = (ys + h_ex);
//         yac_12_(i) = -ys;
//         yac_21_(i) = -ys;

//         // for DC
//         // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.21
//         // except here I only care about the real part, so I remove the "1/j"
//         cplx_type tmp = 1. / cplx_type(x_(i), 0.);
//         ydc_11_(i) = tmp;
//         ydc_22_(i) = tmp;
//         ydc_21_(i) = -tmp;
//         ydc_12_(i) = -tmp;
//     }
// }
