// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BaseNRSolver.h"
#include "BaseNRSolverSingleSlack.h"
#include "DCSolver.h"

#include "SparseLUSolver.h"
#include "KLUSolver.h"
#include "NICSLUSolver.h"

/** Solver based on Newton Raphson, using the SparseLU decomposition of Eigen**/
typedef BaseNRSolver<SparseLULinearSolver> SparseLUSolver;
/** Solver based on Newton Raphson, using the SparseLU decomposition of Eigen, do not consider multiple slack bus**/
typedef BaseNRSolverSingleSlack<SparseLULinearSolver> SparseLUSolverSingleSlack;
/** Solver based on Newton Raphson, using the SparseLU decomposition of Eigen, only suitable for the DC approximation**/
typedef BaseDCSolver<SparseLULinearSolver> DCSolver;

#ifdef KLU_SOLVER_AVAILABLE
/** Solver based on Newton Raphson, using the KLU linear solver**/
typedef BaseNRSolver<KLULinearSolver> KLUSolver;
/** Solver based on Newton Raphson, using the KLU linear solver, do not consider multiple slack bus**/
typedef BaseNRSolverSingleSlack<KLULinearSolver> KLUSolverSingleSlack;
/** Solver based on Newton Raphson, using the KLU linear solver, only suitable for the DC approximation**/
typedef BaseDCSolver<KLULinearSolver> KLUDCSolver;
#endif  // KLU_SOLVER_AVAILABLE

#ifdef NICSLU_SOLVER_AVAILABLE
/** Solver based on Newton Raphson, using the NICSLU linear solver (needs a specific license)**/
typedef BaseNRSolver<NICSLULinearSolver> NICSLUSolver;
/** Solver based on Newton Raphson, using the NICSLU linear solver (needs a specific license), do not consider multiple slack bus**/
typedef BaseNRSolverSingleSlack<NICSLULinearSolver> NICSLUSolverSingleSlack;
/** Solver based on Newton Raphson, using the NICSLU linear solver (needs a specific license), only suitable for the DC approximation**/
typedef BaseDCSolver<NICSLULinearSolver> NICSLUDCSolver;
#endif  // NICSLU_SOLVER_AVAILABLE
