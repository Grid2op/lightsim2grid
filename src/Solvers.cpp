// Copyright (c) 2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "Solvers.h"
#include "GridModel.h"

// these functions use the _gridmodel that is a not a real type in the template class declaration.
// this is why i need to define them here for every specialization.

template<class LinearSolver, FDPFMethod XB_BX>
void BaseFDPFAlgo<LinearSolver, XB_BX>::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, Eigen::SparseMatrix<real_type> & Bpp) const
{
    _gridmodel->fillBp_Bpp(Bp, Bpp, XB_BX);
}

template void FDPF_XB_SparseLUSolver::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, Eigen::SparseMatrix<real_type> & Bpp) const;
template void FDPF_BX_SparseLUSolver::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, Eigen::SparseMatrix<real_type> & Bpp) const;

#ifdef KLU_SOLVER_AVAILABLE
template void FDPF_XB_KLUSolver::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, Eigen::SparseMatrix<real_type> & Bpp) const;
template void FDPF_BX_KLUSolver::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, Eigen::SparseMatrix<real_type> & Bpp) const;
#endif  // KLU_SOLVER_AVAILABLE

#ifdef NICSLU_SOLVER_AVAILABLE
template void FDPF_XB_NICSLUSolver::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, Eigen::SparseMatrix<real_type> & Bpp) const;
template void FDPF_BX_NICSLUSolver::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, Eigen::SparseMatrix<real_type> & Bpp) const;
#endif  // NICSLU_SOLVER_AVAILABLE

#ifdef CKTSO_SOLVER_AVAILABLE
template void FDPF_XB_CKTSOSolver::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, Eigen::SparseMatrix<real_type> & Bpp) const;
template void FDPF_BX_CKTSOSolver::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, Eigen::SparseMatrix<real_type> & Bpp) const;
#endif  // CKTSO_SOLVER_AVAILABLE
