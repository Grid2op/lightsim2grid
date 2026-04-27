// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SOLVERS_H
#define SOLVERS_H

#include "powerflow_algorithm/BaseDCAlgo.hpp"
#include "powerflow_algorithm/NRAlgo.hpp"
#include "powerflow_algorithm/BaseFDPFAlgo.hpp"
#include "powerflow_algorithm/GaussSeidelSynchAlgo.hpp"
#include "powerflow_algorithm/GaussSeidelAlgo.hpp"

#include "linear_solvers/SparseLUSolver.hpp"
#include "linear_solvers/KLUSolver.hpp"
#include "linear_solvers/NICSLUSolver.hpp"
#include "linear_solvers/CKTSOSolver.hpp"

namespace ls2g {

/** Solver based on Newton Raphson, using the SparseLU decomposition of Eigen **/
using SparseLUSolver = NRAlgo<SparseLULinearSolver, MultiSlackNRSystem>;
/** Solver based on Newton Raphson, using the SparseLU decomposition of Eigen, do not consider multiple slack bus **/
using SparseLUSolverSingleSlack = NRAlgo<SparseLULinearSolver, SingleSlackNRSystem>;
/** Solver based on Newton Raphson, using the SparseLU decomposition of Eigen, only suitable for the DC approximation **/
using DCSolver = BaseDCAlgo<SparseLULinearSolver>;
/** Solver based on Fast Decoupled, using the SparseLU decomposition of Eigen **/
using FDPF_XB_SparseLUSolver = BaseFDPFAlgo<SparseLULinearSolver, FDPFMethod::XB>;
using FDPF_BX_SparseLUSolver = BaseFDPFAlgo<SparseLULinearSolver, FDPFMethod::BX>;

#ifdef KLU_SOLVER_AVAILABLE
    /** Solver based on Newton Raphson, using the KLU linear solver **/
    using KLUSolver = NRAlgo<KLULinearSolver, MultiSlackNRSystem>;
    /** Solver based on Newton Raphson, using the KLU linear solver, do not consider multiple slack bus **/
    using KLUSolverSingleSlack = NRAlgo<KLULinearSolver, SingleSlackNRSystem>;
    /** Solver based on Newton Raphson, using the KLU linear solver, only suitable for the DC approximation **/
    using KLUDCSolver = BaseDCAlgo<KLULinearSolver>;
    /** Solver based on Fast Decoupled, using the KLU linear solver **/
    using FDPF_XB_KLUSolver = BaseFDPFAlgo<KLULinearSolver, FDPFMethod::XB>;
    using FDPF_BX_KLUSolver = BaseFDPFAlgo<KLULinearSolver, FDPFMethod::BX>;
#elif defined(_READ_THE_DOCS)
    using KLUSolver = SparseLUSolver;
    using KLUSolverSingleSlack = SparseLUSolverSingleSlack;
    using KLUDCSolver = DCSolver;
    using FDPF_XB_KLUSolver = FDPF_XB_SparseLUSolver;
    using FDPF_BX_KLUSolver = FDPF_BX_SparseLUSolver;
#endif  // KLU_SOLVER_AVAILABLE

#ifdef NICSLU_SOLVER_AVAILABLE
    /** Solver based on Newton Raphson, using the NICSLU linear solver (needs a specific license) **/
    using NICSLUSolver = NRAlgo<NICSLULinearSolver, MultiSlackNRSystem>;
    /** Solver based on Newton Raphson, using the NICSLU linear solver (needs a specific license), do not consider multiple slack bus **/
    using NICSLUSolverSingleSlack = NRAlgo<NICSLULinearSolver, SingleSlackNRSystem>;
    /** Solver based on Newton Raphson, using the NICSLU linear solver (needs a specific license), only suitable for the DC approximation **/
    using NICSLUDCSolver = BaseDCAlgo<NICSLULinearSolver>;
    /** Solver based on Fast Decoupled, using the NICSLU linear solver (needs a specific license) **/
    using FDPF_XB_NICSLUSolver = BaseFDPFAlgo<NICSLULinearSolver, FDPFMethod::XB>;
    using FDPF_BX_NICSLUSolver = BaseFDPFAlgo<NICSLULinearSolver, FDPFMethod::BX>;
#elif defined(_READ_THE_DOCS)
    using NICSLUSolver = SparseLUSolver;
    using NICSLUSolverSingleSlack = SparseLUSolverSingleSlack;
    using NICSLUDCSolver = DCSolver;
    using FDPF_XB_NICSLUSolver = FDPF_XB_SparseLUSolver;
    using FDPF_BX_NICSLUSolver = FDPF_BX_SparseLUSolver;
#endif  // NICSLU_SOLVER_AVAILABLE

#ifdef CKTSO_SOLVER_AVAILABLE
    /** Solver based on Newton Raphson, using the CKTSO linear solver (needs a specific license) **/
    using CKTSOSolver = NRAlgo<CKTSOLinearSolver, MultiSlackNRSystem>;
    /** Solver based on Newton Raphson, using the CKTSO linear solver (needs a specific license), do not consider multiple slack bus **/
    using CKTSOSolverSingleSlack = NRAlgo<CKTSOLinearSolver, SingleSlackNRSystem>;
    /** Solver based on Newton Raphson, using the CKTSO linear solver (needs a specific license), only suitable for the DC approximation **/
    using CKTSODCSolver = BaseDCAlgo<CKTSOLinearSolver>;
    /** Solver based on Fast Decoupled, using the CKTSO linear solver (needs a specific license) **/
    using FDPF_XB_CKTSOSolver = BaseFDPFAlgo<CKTSOLinearSolver, FDPFMethod::XB>;
    using FDPF_BX_CKTSOSolver = BaseFDPFAlgo<CKTSOLinearSolver, FDPFMethod::BX>;
#elif defined(_READ_THE_DOCS)
    using CKTSOSolver = SparseLUSolver;
    using CKTSOSolverSingleSlack = SparseLUSolverSingleSlack;
    using CKTSODCSolver = DCSolver;
    using FDPF_XB_CKTSOSolver = FDPF_XB_SparseLUSolver;
    using FDPF_BX_CKTSOSolver = FDPF_BX_SparseLUSolver;
#endif  // CKTSO_SOLVER_AVAILABLE


#ifndef LS2G_BUILDING_CORE
    extern template class LS2G_API NRAlgo<SparseLULinearSolver, MultiSlackNRSystem>;
    extern template class LS2G_API NRAlgo<SparseLULinearSolver, SingleSlackNRSystem>;
    extern template class LS2G_API BaseDCAlgo<SparseLULinearSolver>;
    extern template class LS2G_API BaseFDPFAlgo<SparseLULinearSolver, FDPFMethod::XB>;
    extern template class LS2G_API BaseFDPFAlgo<SparseLULinearSolver, FDPFMethod::BX>;

#ifdef KLU_SOLVER_AVAILABLE
    extern template class LS2G_API NRAlgo<KLULinearSolver, MultiSlackNRSystem>;
    extern template class LS2G_API NRAlgo<KLULinearSolver, SingleSlackNRSystem>;
    extern template class LS2G_API BaseDCAlgo<KLULinearSolver>;
    extern template class LS2G_API BaseFDPFAlgo<KLULinearSolver, FDPFMethod::XB>;
    extern template class LS2G_API BaseFDPFAlgo<KLULinearSolver, FDPFMethod::BX>;
#endif

#ifdef NICSLU_SOLVER_AVAILABLE
    extern template class LS2G_API NRAlgo<NICSLULinearSolver, MultiSlackNRSystem>;
    extern template class LS2G_API NRAlgo<NICSLULinearSolver, SingleSlackNRSystem>;
    extern template class LS2G_API BaseDCAlgo<NICSLULinearSolver>;
    extern template class LS2G_API BaseFDPFAlgo<NICSLULinearSolver, FDPFMethod::XB>;
    extern template class LS2G_API BaseFDPFAlgo<NICSLULinearSolver, FDPFMethod::BX>;
#endif

#ifdef CKTSO_SOLVER_AVAILABLE
    extern template class LS2G_API NRAlgo<CKTSOLinearSolver, MultiSlackNRSystem>;
    extern template class LS2G_API NRAlgo<CKTSOLinearSolver, SingleSlackNRSystem>;
    extern template class LS2G_API BaseDCAlgo<CKTSOLinearSolver>;
    extern template class LS2G_API BaseFDPFAlgo<CKTSOLinearSolver, FDPFMethod::XB>;
    extern template class LS2G_API BaseFDPFAlgo<CKTSOLinearSolver, FDPFMethod::BX>;
#endif
#endif // LS2G_BUILDING_CORE


} // namespace ls2g

#endif // SOLVERS_H
