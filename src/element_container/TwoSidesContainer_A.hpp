// Copyright (c) 2025, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TWO_SIDES_CONTAINER_A_H
#define TWO_SIDES_CONTAINER_A_H

#include "TwoSidesContainer.hpp"
/**
 * Type of container to represent a line or a transformer.
 * 
 * It has results in amps (A), and some physical properties (r, x and h = g+j.b)
 */
template<class OneSideType>
class TwoSidesContainer_A: public TwoSidesContainer<OneSideType>
{

    protected:
        // physical properties
        RealVect r_;
        RealVect x_;
        CplxVect h_;
};
#endif  // TWO_SIDES_CONTAINER_A_H
