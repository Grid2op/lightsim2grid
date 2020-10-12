// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef UTILS_H
#define UTILS_H

/**
Some typedef and other structures define here and used everywhere else
**/
typedef std::complex<double> cdouble;
//typedef Eigen::Ref<Eigen::VectorXd> EigenPythonNumType;  // Eigen::VectorXd
typedef Eigen::VectorXd EigenPythonNumType;  // Eigen::VectorXd
typedef std::tuple<EigenPythonNumType, EigenPythonNumType, EigenPythonNumType> tuple3d;
typedef std::tuple<EigenPythonNumType, EigenPythonNumType, EigenPythonNumType, EigenPythonNumType> tuple4d;

#endif // UTILS_H
