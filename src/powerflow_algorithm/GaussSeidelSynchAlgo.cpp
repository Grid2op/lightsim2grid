// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GaussSeidelSynchAlgo.h"

void GaussSeidelSynchAlgo::one_iter(CplxVect & tmp_Sbus,
                                      const Eigen::SparseMatrix<cplx_type> & Ybus,
                                      const Eigen::VectorXi & pv,
                                      const Eigen::VectorXi & pq)
{
    // do an update with all nodes being updated at the same time (different than the original GaussSeidel)
    cplx_type tmp;

    const int n_pv = static_cast<int>(pv.size());
    const int n_pq = static_cast<int>(pq.size());

    // CplxVect tmp_YbusV;  // Ybus[k, :] * V
    // CplxVect tmp_conj_Sbus_V;  //  conj(Sbus[k] / V[k])
    CplxVect tmp_YbusV = Ybus * V_;
    CplxVect tmp_conj_Sbus_V = tmp_Sbus.array() / V_.array();
    tmp_conj_Sbus_V = tmp_conj_Sbus_V.array().conjugate();

    // update PQ buses
    for(int k_tmp=0; k_tmp<n_pq; ++k_tmp)
    {
        int k = pq.coeff(k_tmp);
        tmp = (tmp_conj_Sbus_V.coeff(k) -  tmp_YbusV.coeff(k)) / Ybus.coeff(k,k);
        V_.coeffRef(k) += tmp;
    }

    // update PV buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        // update Sbus
        tmp = tmp_YbusV.coeff(k);  // Ybus[k,:] * V
        tmp = std::conj(tmp);  // conj(Ybus[k,:] * V)
        tmp *= V_.coeff(k);  // (V[k] * conj(Ybus[k,:] * V))
        tmp = my_i * std::imag(tmp);
        tmp_Sbus.coeffRef(k) = std::real(tmp_Sbus.coeff(k)) + tmp;

        // update V
        tmp = (tmp_conj_Sbus_V(k) -  tmp_YbusV(k)) / Ybus.coeff(k,k);
        V_.coeffRef(k) += tmp;
    }

    // make sure the voltage magnitudes are not modified at pv buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        V_.coeffRef(k) *= Vm_.coeff(k) / std::abs(V_.coeff(k));
    }
}
