// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BaseConstants.hpp"

const cplx_type BaseConstants::my_i = {0., 1.};
const real_type BaseConstants::my_pi = M_PI;
const real_type BaseConstants::my_one_ = 1.0;
const real_type BaseConstants::my_two_ = 2.0;
const real_type BaseConstants::my_half_ = 0.5;
const real_type BaseConstants::my_zero_ = 0.;
const real_type BaseConstants::my_180_pi_ = 180. / M_PI;
const real_type BaseConstants::v_disco_el_ = 0.0;
const real_type BaseConstants::theta_disco_el_ = 0.0;
const int BaseConstants::_deactivated_bus_id = -1;
const real_type BaseConstants::_tol_equal_float = 1e-7;  // two floats with a difference less than this are equal
const real_type BaseConstants::_1_sqrt_3 = 1.0 / std::sqrt((real_type) 3.); 
