// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GAUSSSEIDELSYNCH_ALGO_H
#define GAUSSSEIDELSYNCH_ALGO_H

#include "GaussSeidelAlgo.h"

/**
The gauss seidel method, where all the updates are happening in a synchronous way, instead of
in a asynchronous way (like for standard gauss seidel)
**/
class GaussSeidelSynchAlgo: public GaussSeidelAlgo
{
    public:
        GaussSeidelSynchAlgo():GaussSeidelAlgo() {};

        ~GaussSeidelSynchAlgo(){}

    protected:
        void one_iter(CplxVect & tmp_Sbus,
                      const Eigen::SparseMatrix<cplx_type> & Ybus,
                      const Eigen::VectorXi & pv,
                      const Eigen::VectorXi & pq
                      );

    private:
        // no copy allowed
        GaussSeidelSynchAlgo( const GaussSeidelSynchAlgo & ) =delete;
        GaussSeidelSynchAlgo & operator=( const GaussSeidelSynchAlgo & )=delete ;

};

#endif // GAUSSSEIDELSYNCH_ALGO_H
