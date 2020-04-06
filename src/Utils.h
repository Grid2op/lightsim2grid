// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of PyKLU2Grid, PyKLU2Grid a implements a c++ backend targeting the Grid2Op platform.

#ifndef UTILS_H
#define UTILS_H

/**
Some typedef and other structures define here and used everywhere else
**/
typedef std::complex<double> cdouble;
typedef std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> tuple3d;
typedef std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> tuple4d;

#endif // UTILS_H