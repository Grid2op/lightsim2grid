// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GaussSeidelSynchAlgo.hpp"

namespace ls2g {

void GaussSeidelSynchAlgo::one_iter(
    Eigen::Ref<CplxVect>            tmp_Sbus,
    const EigenRefConstCplxSpMat    & Ybus,
    const Eigen::Ref<const IntVect> & pv,
    const Eigen::Ref<const IntVect> & pq)
{
    // do an update with all nodes being updated at the same time (different than the original GaussSeidel)
    cplx_type tmp;

    const int n_pv = static_cast<int>(pv.size());
    const int n_pq = static_cast<int>(pq.size());

    // the two working vectors are members, not locals: one_iter runs once per
    // Gauss-Seidel sweep -- hundreds to thousands of times per solve -- and
    // both are exactly nb_bus long every time. noalias(): straight into the
    // buffer, without the temporary Eigen would need for the sparse * dense
    // product. The conjugation is fused into the division rather than done as
    // a second full pass over a vector that was just written.
    tmp_YbusV_.noalias() = Ybus * V_;
    tmp_conj_Sbus_V_ = (tmp_Sbus.array() / V_.array()).conjugate();

    // update PQ buses
    for(int k_tmp=0; k_tmp<n_pq; ++k_tmp)
    {
        int k = pq.coeff(k_tmp);
        // ybus_diag_ (filled by refresh_ybus_cache) instead of Ybus.coeff(k, k),
        // which is a binary search inside column k on every single lookup
        tmp = (tmp_conj_Sbus_V_.coeff(k) -  tmp_YbusV_.coeff(k)) / ybus_diag_.coeff(k);
        V_.coeffRef(k) += tmp;
    }

    // update PV buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        // update Sbus
        tmp = tmp_YbusV_.coeff(k);  // Ybus[k,:] * V
        tmp = std::conj(tmp);  // conj(Ybus[k,:] * V)
        tmp *= V_.coeff(k);  // (V[k] * conj(Ybus[k,:] * V))
        tmp = my_i * std::imag(tmp);
        tmp_Sbus.coeffRef(k) = std::real(tmp_Sbus.coeff(k)) + tmp;

        // update V
        tmp = (tmp_conj_Sbus_V_(k) -  tmp_YbusV_(k)) / ybus_diag_.coeff(k);
        V_.coeffRef(k) += tmp;
    }

    // make sure the voltage magnitudes are not modified at pv buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        V_.coeffRef(k) *= Vm_.coeff(k) / std::abs(V_.coeff(k));
    }
}

} // namespace ls2g
