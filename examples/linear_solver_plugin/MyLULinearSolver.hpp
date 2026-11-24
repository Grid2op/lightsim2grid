// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// MyLULinearSolver: lightsim2grid's linear-solver interface on top of the mylu C API (see
// docs/add_linear_solver.rst for the contract). Header-only, so that the plugin is one
// translation unit; the in-tree version of the same class is split in .hpp / .cpp.

#ifndef MYLU_LINEAR_SOLVER_H
#define MYLU_LINEAR_SOLVER_H

#include <memory>

#include "Utils.hpp"   // ErrorType, RealVect, EigenRefConstRealSpMat
#include "Eigen/Core"
#include "Eigen/SparseCore"

#include "mylu.h"      // has its own extern "C" guard

namespace ls2g {

// No LS2G_API here: that macro is for lightsim2grid_core's own symbols (dllimport outside
// of it, on Windows). This class is compiled into the plugin and never crosses a library
// boundary.
class MyLULinearSolver final
{
    public:
        MyLULinearSolver() noexcept = default;
        ~MyLULinearSolver() noexcept = default;  // handle_ frees itself

        ErrorType reset(){
            handle_.reset();  // the next call is an analyze()
            return ErrorType::NoError;
        }

        ErrorType analyze(const EigenRefConstRealSpMat & J){
            mylu_handle * raw = nullptr;
            if(mylu_create(&raw) != 0) return ErrorType::SolverAnalyze;
            handle_.reset(raw);
            const int n = static_cast<int>(J.cols());
            // J is compressed column storage, 0-based, int indices: handed over as is
            if(mylu_analyze(handle_.get(), n, J.outerIndexPtr(), J.innerIndexPtr()) != 0){
                handle_.reset();
                return ErrorType::SolverAnalyze;
            }
            return ErrorType::NoError;
        }

        ErrorType factorize(const EigenRefConstRealSpMat & J){
            if(!handle_) return ErrorType::NotInitError;
            if(mylu_factor(handle_.get(), J.valuePtr()) != 0) return ErrorType::SolverFactor;
            return ErrorType::NoError;
        }

        ErrorType refactorize(const EigenRefConstRealSpMat & J){
            if(!handle_) return ErrorType::NotInitError;
            // a library with no refactorization of its own can simply `return factorize(J);` here.
            // SolverReFactor, not SolverFactor: what LinearSolverPolicy's fallback (a fresh
            // factorize on the same matrix) answers to, when it is on
            if(mylu_refactor(handle_.get(), J.valuePtr()) != 0) return ErrorType::SolverReFactor;
            return ErrorType::NoError;
        }

        ErrorType solve(Eigen::Ref<RealVect> b){
            if(!handle_) return ErrorType::NotInitError;
            if(mylu_solve(handle_.get(), b.data()) != 0) return ErrorType::SolverSolve;
            if(!b.allFinite()) return ErrorType::InifiniteValue;
            return ErrorType::NoError;
        }

        ErrorType solve_transpose(Eigen::Ref<RealVect> /*b*/){
            return ErrorType::NotImplemented;  // see CAN_SOLVE_TRANSPOSE
        }

        static constexpr bool CAN_SOLVE_MAT = false;
        static constexpr bool CAN_SOLVE_TRANSPOSE = false;

        MyLULinearSolver(const MyLULinearSolver &) = delete;
        MyLULinearSolver(MyLULinearSolver &&) = delete;
        MyLULinearSolver & operator=(const MyLULinearSolver &) = delete;
        MyLULinearSolver & operator=(MyLULinearSolver &&) = delete;

    private:
        // RAII (good practice): the unique_ptr frees the handle, in reset() and in the destructor,
        // whatever happens in between
        struct HandleDeleter {
            void operator()(mylu_handle * h) const noexcept { mylu_destroy(h); }
        };
        std::unique_ptr<mylu_handle, HandleDeleter> handle_;
};

}  // namespace ls2g

#endif  // MYLU_LINEAR_SOLVER_H
