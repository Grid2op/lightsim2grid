// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BaseNRSolver.h"
#include "BaseNRSolverSingleSlack.h"

#include "SparseLUSolver.h"
#include "KLUSolver.h"
#include "NICSLUSolver.h"

typedef BaseNRSolver<SparseLULinearSolver> SparseLUSolver;
typedef BaseNRSolverSingleSlack<SparseLULinearSolver> SparseLUSolverSingleSlack;

#ifdef KLU_SOLVER_AVAILABLE
typedef BaseNRSolver<KLULinearSolver> KLUSolver;
typedef BaseNRSolverSingleSlack<KLULinearSolver> KLUSolverSingleSlack;
#endif  // KLU_SOLVER_AVAILABLE

#ifdef NICSLU_SOLVER_AVAILABLE
typedef BaseNRSolver<NICSLULinearSolver> NICSLUSolver;
typedef BaseNRSolverSingleSlack<NICSLULinearSolver> NICSLUSolverSingleSlack;
#endif  // NICSLU_SOLVER_AVAILABLE

