// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef CUSTTIMER_H
#define CUSTTIMER_H

#include <chrono>

/**

This class presents a basic timer that is used in KLUSolver to know on which part of the solver
most time were taken.

**/
class CustTimer{
    public:
        CustTimer():start_(std::chrono::steady_clock::now()), end_(start_){};

        double duration(){
            end_ = std::chrono::steady_clock::now();
            std::chrono::duration<double> res = end_ - start_;
            return res.count();
        }
    private:
        std::chrono::time_point<std::chrono::steady_clock> start_;
        std::chrono::time_point<std::chrono::steady_clock> end_;
};

#endif //CUSTTIMER_H
