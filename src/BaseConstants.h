// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASECONSTANTS_H
#define BASECONSTANTS_H

#include "Utils.h"

/**
Definition of some basic constants, because sometimes Eigen cannot deduce types.
Eg if I type "1.0" then Eigen cast it to "double" and i cannot use it with real_type = float for example
**/
class BaseConstants
{
    public:
        static const cplx_type my_i;
        static const real_type my_pi;
        static const real_type my_half_;
        static const real_type my_one_;
        static const real_type my_two_;
        static const real_type my_zero_;
        static const real_type my_180_pi_;
};

enum class FDPFMethod {XB, BX};  // Different type of FDPF powerflow
// FDPFMethod::XB => alg = 2 in pypower / pandapower
// FDPFMethod::BX => alg = 3 in pypower / pandapower

#endif // BASECONSTANTS_H
