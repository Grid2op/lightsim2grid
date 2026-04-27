// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ALGO_CONFIG_H
#define ALGO_CONFIG_H

#include <vector>

namespace ls2g {

/**
 * Serialisable configuration blob for any BaseAlgo subclass.
 * Used by get_config() / set_config() to support pickling of algorithm parameters.
 */
struct AlgoConfig {
    std::vector<int>    int_params;
    std::vector<double> real_params;
};

} // namespace ls2g

#endif // ALGO_CONFIG_H
