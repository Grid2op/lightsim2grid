// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GAUSSSEIDELSYNCH_ALGO_H
#define GAUSSSEIDELSYNCH_ALGO_H

#include "GaussSeidelAlgo.hpp"

namespace ls2g {

/**
The gauss seidel method, where all the updates are happening in a synchronous way, instead of
in a asynchronous way (like for standard gauss seidel)
**/
class LS2G_API GaussSeidelSynchAlgo final: public GaussSeidelAlgo
{
    public:
        GaussSeidelSynchAlgo() noexcept : GaussSeidelAlgo() {};

        ~GaussSeidelSynchAlgo() noexcept override = default;

        void reset() override {
            GaussSeidelAlgo::reset();
            tmp_YbusV_ = CplxVect();
            tmp_conj_Sbus_V_ = CplxVect();
        }

    protected:
        void one_iter(
            Eigen::Ref<CplxVect>            tmp_Sbus,
            const EigenRefConstCplxSpMat    & Ybus,
            const Eigen::Ref<const IntVect> & pv,
            const Eigen::Ref<const IntVect> & pq
        ) override;

        // per-sweep scratch, kept across sweeps (see one_iter)
        CplxVect tmp_YbusV_;        // Ybus * V
        CplxVect tmp_conj_Sbus_V_;  // conj(Sbus / V)

    private:
        // no copy allowed
        GaussSeidelSynchAlgo(const GaussSeidelSynchAlgo&) = delete;
        GaussSeidelSynchAlgo(GaussSeidelSynchAlgo&&) = delete;
        GaussSeidelSynchAlgo & operator=(GaussSeidelSynchAlgo&&) = delete;
        GaussSeidelSynchAlgo & operator=(const GaussSeidelSynchAlgo&) = delete;

};


} // namespace ls2g

#endif // GAUSSSEIDELSYNCH_ALGO_H