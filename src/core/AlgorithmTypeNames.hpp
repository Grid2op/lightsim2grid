// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ALGORITHM_TYPE_NAMES_H
#define ALGORITHM_TYPE_NAMES_H

#include <string>

namespace ls2g {

// Forward-declare to avoid including AlgorithmSelector.hpp here.
// AlgorithmType is defined in AlgorithmSelector.hpp; callers must include it first.
enum class AlgorithmType;

// Returns the registry name for a known AlgorithmType enum value.
// Throws std::runtime_error for AlgorithmType::Custom or unknown values.
std::string algo_type_to_name(AlgorithmType t);

// Returns the AlgorithmType for a known registry name.
// Returns AlgorithmType::Custom for names not in the known enum set.
AlgorithmType name_to_algo_type(const std::string& name);


} // namespace ls2g

#endif // ALGORITHM_TYPE_NAMES_H
