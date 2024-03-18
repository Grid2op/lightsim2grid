// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "Utils.h"

std::ostream& operator<<(std::ostream& out, const ErrorType & error_type){
    switch (error_type)
    {
    case ErrorType::NoError:
        out << "NoError";
        break;
    case ErrorType::SingularMatrix:
        out << "SingularMatrix";
        break;
    case ErrorType::TooManyIterations:
        out << "TooManyIterations";
        break;
    case ErrorType::InifiniteValue:
        out << "InifiniteValue";
        break;
    case ErrorType::SolverAnalyze:
        out << "SolverAnalyze";
        break;
    case ErrorType::SolverFactor:
        out << "SolverFactor";
        break;
    case ErrorType::SolverReFactor:
        out << "SolverReFactor";
        break;
    case ErrorType::SolverSolve:
        out << "SolverSolve";
        break;
    case ErrorType::NotInitError:
        out << "NotInitError";
        break;
    case ErrorType::LicenseError:
        out << "LicenseError";
        break;
    default:
        out << "unknown error (check utils.cpp)";
        break;
    }
    return out;
}
