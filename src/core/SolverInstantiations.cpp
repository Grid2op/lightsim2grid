// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "Solvers.hpp"
#include "GridModel.hpp"   // required: fillBp_Bpp body references gridmodel_ptr_ (a GridModel*)

namespace ls2g {

// fillBp_Bpp cannot be defined in BaseFDPFAlgo.tpp because it calls GridModel methods
// and GridModel is only a forward declaration there.
template<class LinearSolver, FDPFMethod XB_BX>
void BaseFDPFAlgo<LinearSolver, XB_BX>::fillBp_Bpp(
    Eigen::SparseMatrix<real_type> & Bp,
    Eigen::SparseMatrix<real_type> & Bpp) const
{
    gridmodel_ptr_->fillBp_Bpp(Bp, Bpp, XB_BX);
}

// ---- SparseLU (always available) ----
template class BaseNRAlgo<SparseLULinearSolver>;
template class BaseNRSingleSlackAlgo<SparseLULinearSolver>;
template class BaseDCAlgo<SparseLULinearSolver>;
template class BaseFDPFAlgo<SparseLULinearSolver, FDPFMethod::XB>;
template class BaseFDPFAlgo<SparseLULinearSolver, FDPFMethod::BX>;

// ---- KLU (optional) ----
#ifdef KLU_SOLVER_AVAILABLE
template class BaseNRAlgo<KLULinearSolver>;
template class BaseNRSingleSlackAlgo<KLULinearSolver>;
template class BaseDCAlgo<KLULinearSolver>;
template class BaseFDPFAlgo<KLULinearSolver, FDPFMethod::XB>;
template class BaseFDPFAlgo<KLULinearSolver, FDPFMethod::BX>;
#endif

// ---- NICSLU (optional) ----
#ifdef NICSLU_SOLVER_AVAILABLE
template class BaseNRAlgo<NICSLULinearSolver>;
template class BaseNRSingleSlackAlgo<NICSLULinearSolver>;
template class BaseDCAlgo<NICSLULinearSolver>;
template class BaseFDPFAlgo<NICSLULinearSolver, FDPFMethod::XB>;
template class BaseFDPFAlgo<NICSLULinearSolver, FDPFMethod::BX>;
#endif

// ---- CKTSO (optional) ----
#ifdef CKTSO_SOLVER_AVAILABLE
template class BaseNRAlgo<CKTSOLinearSolver>;
template class BaseNRSingleSlackAlgo<CKTSOLinearSolver>;
template class BaseDCAlgo<CKTSOLinearSolver>;
template class BaseFDPFAlgo<CKTSOLinearSolver, FDPFMethod::XB>;
template class BaseFDPFAlgo<CKTSOLinearSolver, FDPFMethod::BX>;
#endif

} // namespace ls2g
