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

/** Newton-Raphson (multi-slack) with Eigen SparseLU linear solver **/
using NR_SparseLU = NRAlgo<SparseLULinearSolver, MultiSlackNRSystem>;
/** Newton-Raphson (single-slack) with Eigen SparseLU linear solver **/
using NRSing_SparseLU = NRAlgo<SparseLULinearSolver, SingleSlackNRSystem>;
/** DC approximation with Eigen SparseLU linear solver **/
using DC_SparseLU = BaseDCAlgo<SparseLULinearSolver>;
/** Fast-Decoupled Power Flow (XB variant) with Eigen SparseLU linear solver **/
using FDPF_XB_SparseLU = BaseFDPFAlgo<SparseLULinearSolver, FDPFMethod::XB>;
/** Fast-Decoupled Power Flow (BX variant) with Eigen SparseLU linear solver **/
using FDPF_BX_SparseLU = BaseFDPFAlgo<SparseLULinearSolver, FDPFMethod::BX>;

#ifdef KLU_SOLVER_AVAILABLE
    /** Newton-Raphson (multi-slack) with KLU linear solver **/
    using NR_KLU = NRAlgo<KLULinearSolver, MultiSlackNRSystem>;
    /** Newton-Raphson (single-slack) with KLU linear solver **/
    using NRSing_KLU = NRAlgo<KLULinearSolver, SingleSlackNRSystem>;
    /** DC approximation with KLU linear solver **/
    using DC_KLU = BaseDCAlgo<KLULinearSolver>;
    /** Fast-Decoupled Power Flow (XB variant) with KLU linear solver **/
    using FDPF_XB_KLU = BaseFDPFAlgo<KLULinearSolver, FDPFMethod::XB>;
    /** Fast-Decoupled Power Flow (BX variant) with KLU linear solver **/
    using FDPF_BX_KLU = BaseFDPFAlgo<KLULinearSolver, FDPFMethod::BX>;
#elif defined(_READ_THE_DOCS)
    using NR_KLU = NR_SparseLU;
    using NRSing_KLU = NRSing_SparseLU;
    using DC_KLU = DC_SparseLU;
    using FDPF_XB_KLU = FDPF_XB_SparseLU;
    using FDPF_BX_KLU = FDPF_BX_SparseLU;
#endif  // KLU_SOLVER_AVAILABLE

#ifdef NICSLU_SOLVER_AVAILABLE
    /** Newton-Raphson (multi-slack) with NICSLU linear solver (requires license) **/
    using NR_NICSLU = NRAlgo<NICSLULinearSolver, MultiSlackNRSystem>;
    /** Newton-Raphson (single-slack) with NICSLU linear solver (requires license) **/
    using NRSing_NICSLU = NRAlgo<NICSLULinearSolver, SingleSlackNRSystem>;
    /** DC approximation with NICSLU linear solver (requires license) **/
    using DC_NICSLU = BaseDCAlgo<NICSLULinearSolver>;
    /** Fast-Decoupled Power Flow (XB variant) with NICSLU linear solver (requires license) **/
    using FDPF_XB_NICSLU = BaseFDPFAlgo<NICSLULinearSolver, FDPFMethod::XB>;
    /** Fast-Decoupled Power Flow (BX variant) with NICSLU linear solver (requires license) **/
    using FDPF_BX_NICSLU = BaseFDPFAlgo<NICSLULinearSolver, FDPFMethod::BX>;
#elif defined(_READ_THE_DOCS)
    using NR_NICSLU = NR_SparseLU;
    using NRSing_NICSLU = NRSing_SparseLU;
    using DC_NICSLU = DC_SparseLU;
    using FDPF_XB_NICSLU = FDPF_XB_SparseLU;
    using FDPF_BX_NICSLU = FDPF_BX_SparseLU;
#endif  // NICSLU_SOLVER_AVAILABLE

#ifdef CKTSO_SOLVER_AVAILABLE
    /** Newton-Raphson (multi-slack) with CKTSO linear solver (requires license) **/
    using NR_CKTSO = NRAlgo<CKTSOLinearSolver, MultiSlackNRSystem>;
    /** Newton-Raphson (single-slack) with CKTSO linear solver (requires license) **/
    using NRSing_CKTSO = NRAlgo<CKTSOLinearSolver, SingleSlackNRSystem>;
    /** DC approximation with CKTSO linear solver (requires license) **/
    using DC_CKTSO = BaseDCAlgo<CKTSOLinearSolver>;
    /** Fast-Decoupled Power Flow (XB variant) with CKTSO linear solver (requires license) **/
    using FDPF_XB_CKTSO = BaseFDPFAlgo<CKTSOLinearSolver, FDPFMethod::XB>;
    /** Fast-Decoupled Power Flow (BX variant) with CKTSO linear solver (requires license) **/
    using FDPF_BX_CKTSO = BaseFDPFAlgo<CKTSOLinearSolver, FDPFMethod::BX>;
#elif defined(_READ_THE_DOCS)
    using NR_CKTSO = NR_SparseLU;
    using NRSing_CKTSO = NRSing_SparseLU;
    using DC_CKTSO = DC_SparseLU;
    using FDPF_XB_CKTSO = FDPF_XB_SparseLU;
    using FDPF_BX_CKTSO = FDPF_BX_SparseLU;
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
