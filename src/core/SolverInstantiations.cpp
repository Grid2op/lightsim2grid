// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "Solvers.hpp"
#include <stdexcept>

#include "LSGrid.hpp"   // required: fillBp_Bpp body references lsgrid_ptr_ (a GridModel*)

namespace ls2g {

// fillBp_Bpp cannot be defined in BaseFDPFAlgo.tpp because it calls GridModel methods
// and GridModel is only a forward declaration there.
template<class LinearSolver, FDPFMethod XB_BX>
void BaseFDPFAlgo<LinearSolver, XB_BX>::fillBp_Bpp(
    Eigen::SparseMatrix<real_type> & Bp,
    Eigen::SparseMatrix<real_type> & Bpp) const
{
    if (lsgrid_ptr_ == nullptr) {
        // standalone use (eg a bare FDPF_XB_SparseLU() built from python):
        // unlike NR / Gauss-Seidel, FDPF cannot work from Ybus alone -- the
        // B' / B'' matrices are built from the grid data. Without this check
        // the call below is a null-pointer dereference (segfault).
        throw std::runtime_error(
            "BaseFDPFAlgo::compute_pf: this fast-decoupled solver is not attached to an LSGrid "
            "(the B'/B'' matrices are built from the grid data, not from Ybus alone), so it "
            "cannot be used standalone: use it through LSGrid instead "
            "(change_algorithm(...) then ac_pf(...)).");
    }
    lsgrid_ptr_->fillBp_Bpp(Bp, Bpp, XB_BX);
}

// ---- SparseLU (always available) ----
template class LS2G_API NRAlgo<LinearSolverPolicy<SparseLULinearSolver>, MultiSlackNRSystem>;
template class LS2G_API NRAlgo<LinearSolverPolicy<SparseLULinearSolver>, SingleSlackNRSystem>;
template class LS2G_API BaseDCAlgo<LinearSolverPolicy<SparseLULinearSolver>>;
template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<SparseLULinearSolver>, FDPFMethod::XB>;
template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<SparseLULinearSolver>, FDPFMethod::BX>;

// ---- KLU (optional) ----
#ifdef KLU_SOLVER_AVAILABLE
template class LS2G_API NRAlgo<LinearSolverPolicy<KLULinearSolver>, MultiSlackNRSystem>;
template class LS2G_API NRAlgo<LinearSolverPolicy<KLULinearSolver>, SingleSlackNRSystem>;
template class LS2G_API BaseDCAlgo<LinearSolverPolicy<KLULinearSolver>>;
template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<KLULinearSolver>, FDPFMethod::XB>;
template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<KLULinearSolver>, FDPFMethod::BX>;
template class LS2G_API NRAlgo<RefactorRetryLinearSolver<KLULinearSolver>, MultiSlackNRSystem>;
#endif

// ---- NICSLU (optional) ----
#ifdef NICSLU_SOLVER_AVAILABLE
template class LS2G_API NRAlgo<LinearSolverPolicy<NICSLULinearSolver>, MultiSlackNRSystem>;
template class LS2G_API NRAlgo<LinearSolverPolicy<NICSLULinearSolver>, SingleSlackNRSystem>;
template class LS2G_API BaseDCAlgo<LinearSolverPolicy<NICSLULinearSolver>>;
template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<NICSLULinearSolver>, FDPFMethod::XB>;
template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<NICSLULinearSolver>, FDPFMethod::BX>;
template class LS2G_API NRAlgo<RefactorRetryLinearSolver<NICSLULinearSolver>, MultiSlackNRSystem>;
#endif

// ---- CKTSO (optional) ----
#ifdef CKTSO_SOLVER_AVAILABLE
template class LS2G_API NRAlgo<LinearSolverPolicy<CKTSOLinearSolver>, MultiSlackNRSystem>;
template class LS2G_API NRAlgo<LinearSolverPolicy<CKTSOLinearSolver>, SingleSlackNRSystem>;
template class LS2G_API BaseDCAlgo<LinearSolverPolicy<CKTSOLinearSolver>>;
template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<CKTSOLinearSolver>, FDPFMethod::XB>;
template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<CKTSOLinearSolver>, FDPFMethod::BX>;
template class LS2G_API NRAlgo<RefactorRetryLinearSolver<CKTSOLinearSolver>, MultiSlackNRSystem>;
#endif

} // namespace ls2g
