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
#include "linear_solvers/LinearSolverPolicy.hpp"
#include "linear_solvers/RefactorRetryLinearSolver.hpp"

namespace ls2g {

// Every built-in solver below is wrapped in LinearSolverPolicy<...>: a transparent
// pass-through that counts and times each analyze/factorize/refactorize/solve call
// (see LinearSolverPolicy.hpp). NRAlgo/BaseDCAlgo/BaseFDPFAlgo source their own
// get_timers_jacobian()/get_linear_solver_stats() from it instead of keeping separate
// bookkeeping.

/** Newton-Raphson (multi-slack) with Eigen SparseLU linear solver **/
using NR_SparseLU = NRAlgo<LinearSolverPolicy<SparseLULinearSolver>, MultiSlackNRSystem>;
/** Newton-Raphson (single-slack) with Eigen SparseLU linear solver **/
using NRSing_SparseLU = NRAlgo<LinearSolverPolicy<SparseLULinearSolver>, SingleSlackNRSystem>;
/** DC approximation with Eigen SparseLU linear solver **/
using DC_SparseLU = BaseDCAlgo<LinearSolverPolicy<SparseLULinearSolver>>;
/** Fast-Decoupled Power Flow (XB variant) with Eigen SparseLU linear solver **/
using FDPF_XB_SparseLU = BaseFDPFAlgo<LinearSolverPolicy<SparseLULinearSolver>, FDPFMethod::XB>;
/** Fast-Decoupled Power Flow (BX variant) with Eigen SparseLU linear solver **/
using FDPF_BX_SparseLU = BaseFDPFAlgo<LinearSolverPolicy<SparseLULinearSolver>, FDPFMethod::BX>;
// NB: SparseLU's factorize()==refactorize() already (Eigen has no cheaper "reuse pivot
// order" refactor), so there is deliberately no NRRefactorRetry_SparseLU -- wrapping it
// would be a no-op retry.

#if defined(KLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    // Under _READ_THE_DOCS without the real library, KLULinearSolver (see
    // linear_solvers/KLUSolver.hpp) is a doc-only, no-op stand-in with the same name
    // and interface as the real one -- so every alias below is still a genuine,
    // distinct, fully-functional-looking C++ type, not a copy of NR_SparseLU.
    /** Newton-Raphson (multi-slack) with KLU linear solver **/
    using NR_KLU = NRAlgo<LinearSolverPolicy<KLULinearSolver>, MultiSlackNRSystem>;
    /** Newton-Raphson (single-slack) with KLU linear solver **/
    using NRSing_KLU = NRAlgo<LinearSolverPolicy<KLULinearSolver>, SingleSlackNRSystem>;
    /** DC approximation with KLU linear solver **/
    using DC_KLU = BaseDCAlgo<LinearSolverPolicy<KLULinearSolver>>;
    /** Fast-Decoupled Power Flow (XB variant) with KLU linear solver **/
    using FDPF_XB_KLU = BaseFDPFAlgo<LinearSolverPolicy<KLULinearSolver>, FDPFMethod::XB>;
    /** Fast-Decoupled Power Flow (BX variant) with KLU linear solver **/
    using FDPF_BX_KLU = BaseFDPFAlgo<LinearSolverPolicy<KLULinearSolver>, FDPFMethod::BX>;
    /** Newton-Raphson (multi-slack) with KLU linear solver, retrying a failed refactor
     *  with a full factorize() before giving up **/
    using NRRefactorRetry_KLU = NRAlgo<RefactorRetryLinearSolver<KLULinearSolver>, MultiSlackNRSystem>;
#endif  // KLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

#if defined(NICSLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    // See the KLU comment above: NICSLULinearSolver is a doc-only no-op stand-in under
    // _READ_THE_DOCS without the real (licensed) library.
    /** Newton-Raphson (multi-slack) with NICSLU linear solver (requires license) **/
    using NR_NICSLU = NRAlgo<LinearSolverPolicy<NICSLULinearSolver>, MultiSlackNRSystem>;
    /** Newton-Raphson (single-slack) with NICSLU linear solver (requires license) **/
    using NRSing_NICSLU = NRAlgo<LinearSolverPolicy<NICSLULinearSolver>, SingleSlackNRSystem>;
    /** DC approximation with NICSLU linear solver (requires license) **/
    using DC_NICSLU = BaseDCAlgo<LinearSolverPolicy<NICSLULinearSolver>>;
    /** Fast-Decoupled Power Flow (XB variant) with NICSLU linear solver (requires license) **/
    using FDPF_XB_NICSLU = BaseFDPFAlgo<LinearSolverPolicy<NICSLULinearSolver>, FDPFMethod::XB>;
    /** Fast-Decoupled Power Flow (BX variant) with NICSLU linear solver (requires license) **/
    using FDPF_BX_NICSLU = BaseFDPFAlgo<LinearSolverPolicy<NICSLULinearSolver>, FDPFMethod::BX>;
    /** Newton-Raphson (multi-slack) with NICSLU linear solver, retrying a failed refactor
     *  with a full factorize() before giving up (requires license) **/
    using NRRefactorRetry_NICSLU = NRAlgo<RefactorRetryLinearSolver<NICSLULinearSolver>, MultiSlackNRSystem>;
#endif  // NICSLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

#if defined(CKTSO_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    // See the KLU comment above: CKTSOLinearSolver is a doc-only no-op stand-in under
    // _READ_THE_DOCS without the real (licensed) library.
    /** Newton-Raphson (multi-slack) with CKTSO linear solver (requires license) **/
    using NR_CKTSO = NRAlgo<LinearSolverPolicy<CKTSOLinearSolver>, MultiSlackNRSystem>;
    /** Newton-Raphson (single-slack) with CKTSO linear solver (requires license) **/
    using NRSing_CKTSO = NRAlgo<LinearSolverPolicy<CKTSOLinearSolver>, SingleSlackNRSystem>;
    /** DC approximation with CKTSO linear solver (requires license) **/
    using DC_CKTSO = BaseDCAlgo<LinearSolverPolicy<CKTSOLinearSolver>>;
    /** Fast-Decoupled Power Flow (XB variant) with CKTSO linear solver (requires license) **/
    using FDPF_XB_CKTSO = BaseFDPFAlgo<LinearSolverPolicy<CKTSOLinearSolver>, FDPFMethod::XB>;
    /** Fast-Decoupled Power Flow (BX variant) with CKTSO linear solver (requires license) **/
    using FDPF_BX_CKTSO = BaseFDPFAlgo<LinearSolverPolicy<CKTSOLinearSolver>, FDPFMethod::BX>;
    /** Newton-Raphson (multi-slack) with CKTSO linear solver, retrying a failed refactor
     *  with a full factorize() before giving up (requires license) **/
    using NRRefactorRetry_CKTSO = NRAlgo<RefactorRetryLinearSolver<CKTSOLinearSolver>, MultiSlackNRSystem>;
#endif  // CKTSO_SOLVER_AVAILABLE (or _READ_THE_DOCS)


#ifndef LS2G_BUILDING_CORE
    extern template class LS2G_API NRAlgo<LinearSolverPolicy<SparseLULinearSolver>, MultiSlackNRSystem>;
    extern template class LS2G_API NRAlgo<LinearSolverPolicy<SparseLULinearSolver>, SingleSlackNRSystem>;
    extern template class LS2G_API BaseDCAlgo<LinearSolverPolicy<SparseLULinearSolver>>;
    extern template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<SparseLULinearSolver>, FDPFMethod::XB>;
    extern template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<SparseLULinearSolver>, FDPFMethod::BX>;

#if defined(KLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    extern template class LS2G_API NRAlgo<LinearSolverPolicy<KLULinearSolver>, MultiSlackNRSystem>;
    extern template class LS2G_API NRAlgo<LinearSolverPolicy<KLULinearSolver>, SingleSlackNRSystem>;
    extern template class LS2G_API BaseDCAlgo<LinearSolverPolicy<KLULinearSolver>>;
    extern template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<KLULinearSolver>, FDPFMethod::XB>;
    extern template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<KLULinearSolver>, FDPFMethod::BX>;
    extern template class LS2G_API NRAlgo<RefactorRetryLinearSolver<KLULinearSolver>, MultiSlackNRSystem>;
#endif

#if defined(NICSLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    extern template class LS2G_API NRAlgo<LinearSolverPolicy<NICSLULinearSolver>, MultiSlackNRSystem>;
    extern template class LS2G_API NRAlgo<LinearSolverPolicy<NICSLULinearSolver>, SingleSlackNRSystem>;
    extern template class LS2G_API BaseDCAlgo<LinearSolverPolicy<NICSLULinearSolver>>;
    extern template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<NICSLULinearSolver>, FDPFMethod::XB>;
    extern template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<NICSLULinearSolver>, FDPFMethod::BX>;
    extern template class LS2G_API NRAlgo<RefactorRetryLinearSolver<NICSLULinearSolver>, MultiSlackNRSystem>;
#endif

#if defined(CKTSO_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    extern template class LS2G_API NRAlgo<LinearSolverPolicy<CKTSOLinearSolver>, MultiSlackNRSystem>;
    extern template class LS2G_API NRAlgo<LinearSolverPolicy<CKTSOLinearSolver>, SingleSlackNRSystem>;
    extern template class LS2G_API BaseDCAlgo<LinearSolverPolicy<CKTSOLinearSolver>>;
    extern template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<CKTSOLinearSolver>, FDPFMethod::XB>;
    extern template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<CKTSOLinearSolver>, FDPFMethod::BX>;
    extern template class LS2G_API NRAlgo<RefactorRetryLinearSolver<CKTSOLinearSolver>, MultiSlackNRSystem>;
#endif
#endif // LS2G_BUILDING_CORE


} // namespace ls2g

#endif // SOLVERS_H
