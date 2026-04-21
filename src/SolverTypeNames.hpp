// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SOLVER_TYPE_NAMES_H
#define SOLVER_TYPE_NAMES_H

#include <string>

// Forward-declare to avoid including ChooseSolver.hpp here.
// SolverType is defined in ChooseSolver.hpp; callers must include it first.
enum class SolverType;

// Returns the registry name for a known SolverType enum value.
// Throws std::runtime_error for SolverType::Custom or unknown values.
std::string solver_type_to_name(SolverType t);

// Returns the SolverType for a known registry name.
// Returns SolverType::Custom for names not in the known enum set.
SolverType name_to_solver_type(const std::string& name);

#endif // SOLVER_TYPE_NAMES_H
