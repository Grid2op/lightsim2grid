// Copyright (c) 2020-2024, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "ChooseSolver.h"
#include "DataConverter.h"
#include "GridModel.h"

#include "batch_algorithm/TimeSeries.h"
#include "batch_algorithm/ContingencyAnalysis.h"

#include "help_fun_msg.h"

#ifndef KLU_SOLVER_AVAILABLE
#define this_KLU_SOLVER_AVAILABLE 0
#else 
#define this_KLU_SOLVER_AVAILABLE 1
#endif
#ifndef NICSLU_SOLVER_AVAILABLE
#define this_NICSLU_SOLVER_AVAILABLE 0
#else 
#define this_NICSLU_SOLVER_AVAILABLE 1
#endif
#ifndef CKTSO_SOLVER_AVAILABLE
#define this_CKTSO_SOLVER_AVAILABLE 0
#else 
#define this_CKTSO_SOLVER_AVAILABLE 1
#endif
#ifndef __COMPILE_MARCHNATIVE
#define this__COMPILE_MARCHNATIVE 0
#else 
#define this__COMPILE_MARCHNATIVE 1
#endif
#ifndef __O3_OPTIM
#define this__O3_OPTIM 0
#else 
#define this__O3_OPTIM 1
#endif
#ifndef VERSION
#define this_VERSION "unknown"
#else 
#define this_VERSION VERSION
#endif
#ifdef NICSLU_PATH
#define this_NICSLU_PATH NICSLU_PATH
#endif
#ifdef CKTSO_PATH
#define this_CKTSO_PATH CKTSO_PATH
#endif

namespace py = pybind11;

PYBIND11_MODULE(lightsim2grid_cpp, m)
{

    // constant and compilation information
    m.attr("klu_solver_available") = py::bool_(this_KLU_SOLVER_AVAILABLE);
    m.attr("nicslu_solver_available") = py::bool_(this_NICSLU_SOLVER_AVAILABLE);
    m.attr("cktso_solver_available") = py::bool_(this_CKTSO_SOLVER_AVAILABLE);
    m.attr("compiled_march_native") = py::bool_(this__COMPILE_MARCHNATIVE);
    m.attr("compiled_o3_optim") = py::bool_(this__O3_OPTIM);
    m.attr("version") = py::str(this_VERSION);
    #ifdef NICSLU_PATH
    m.attr("nicslu_lib") = py::str(this_NICSLU_PATH);
    #endif
    #ifdef CKTSO_PATH
    m.attr("cktso_lib") = py::str(this_CKTSO_PATH);
    #endif
    
    // solver method for FDPF
    py::enum_<FDPFMethod>(m, "FDPFMethod", "This enum controls the type of method you can use for Fast Decoupled Powerflow (XB or BX)")
        .value("XB", FDPFMethod::XB, "denotes the XB method")
        .value("BX", FDPFMethod::BX, "denotes the BX method")
        .export_values();

    // solvers
    py::enum_<SolverType>(m, "SolverType", "This enum controls the solver you want to use.")
        .value("GaussSeidel", SolverType::GaussSeidel, "denotes the :class:`lightsim2grid.solver.GaussSeidelSolver`")
        .value("GaussSeidelSynch", SolverType::GaussSeidelSynch, "denotes the :class:`lightsim2grid.solver.GaussSeidelSynchSolver`")
        .value("SparseLU", SolverType::SparseLU, "denotes the :class:`lightsim2grid.solver.SparseLUSolver`")
        .value("SparseLUSingleSlack", SolverType::SparseLUSingleSlack, "denotes the :class:`lightsim2grid.solver.SparseLUSolverSingleSlack`")
        .value("DC", SolverType::DC, "denotes the :class:`lightsim2grid.solver.DCSolver`")
        .value("KLU", SolverType::KLU, "denotes the :class:`lightsim2grid.solver.KLUSolver`")
        .value("KLUSingleSlack", SolverType::KLUSingleSlack, "denotes the :class:`lightsim2grid.solver.KLUSolverSingleSlack`")
        .value("KLUDC", SolverType::KLUDC, "denotes the :class:`lightsim2grid.solver.KLUDCSolver`")
        .value("NICSLU", SolverType::NICSLU, "denotes the :class:`lightsim2grid.solver.NICSLUSolver`")
        .value("NICSLUSingleSlack", SolverType::NICSLUSingleSlack, "denotes the :class:`lightsim2grid.solver.NICSLUSolverSingleSlack`")
        .value("NICSLUDC", SolverType::NICSLUDC, "denotes the :class:`lightsim2grid.solver.NICSLUDCSolver`")
        .value("CKTSO", SolverType::CKTSO, "denotes the :class:`lightsim2grid.solver.CKTSOSolver`")
        .value("CKTSOSingleSlack", SolverType::CKTSOSingleSlack, "denotes the :class:`lightsim2grid.solver.CKTSOSolverSingleSlack`")
        .value("CKTSODC", SolverType::CKTSODC, "denotes the :class:`lightsim2grid.solver.CKTSODCSolver`")
        .value("FDPF_XB_SparseLU", SolverType::FDPF_XB_SparseLU, "denotes the :class:`lightsim2grid.solver.FDPF_XB_SparseLUSolver`")
        .value("FDPF_BX_SparseLU", SolverType::FDPF_BX_SparseLU, "denotes the :class:`lightsim2grid.solver.FDPF_BX_SparseLUSolver`")
        .value("FDPF_XB_KLU", SolverType::FDPF_XB_KLU, "denotes the :class:`lightsim2grid.solver.FDPF_XB_KLUSolver`")
        .value("FDPF_BX_KLU", SolverType::FDPF_BX_KLU, "denotes the :class:`lightsim2grid.solver.FDPF_BX_KLUSolver`")
        .value("FDPF_XB_NICSLU", SolverType::FDPF_XB_NICSLU, "denotes the :class:`lightsim2grid.solver.FDPF_XB_NICSLUSolver`")
        .value("FDPF_BX_NICSLU", SolverType::FDPF_BX_NICSLU, "denotes the :class:`lightsim2grid.solver.FDPF_BX_NICSLUSolver`")
        .value("FDPF_XB_CKTSO", SolverType::FDPF_XB_CKTSO, "denotes the :class:`lightsim2grid.solver.FDPF_XB_CKTSOSolver`")
        .value("FDPF_BX_CKTSO", SolverType::FDPF_BX_CKTSO, "denotes the :class:`lightsim2grid.solver.FDPF_BX_CKTSOSolver`")
        .export_values();

    py::enum_<ErrorType>(m, "ErrorType", "This enum controls the error encountered in the solver")
        .value("NoError", ErrorType::NoError, "No error were encountered")
        .value("SingularMatrix", ErrorType::SingularMatrix, "The Jacobian matrix was singular and could not be factorized (most likely, the grid is not connex)")
        .value("TooManyIterations", ErrorType::TooManyIterations, "The solver reached the maximum number of iterations allowed")
        .value("InifiniteValue", ErrorType::InifiniteValue, "Some infinite values were encountered in the update vector (to update Vm or Va)")
        .value("SolverAnalyze", ErrorType::SolverAnalyze, "The linear solver failed at the 'analyze' step (*eg* `analyzePattern` for Eigen, `klu_analyze` for KLU or `Initialize` for NICSLU")
        .value("SolverFactor", ErrorType::SolverFactor, "The linear solver failed to factor the jacobian matrix (*eg* `factorize` for Eigen (first call), `klu_factor` for KLU or `FactorizeMatrix` for NICSLU (first call)")
        .value("SolverReFactor", ErrorType::SolverReFactor, "The linear solver failed to (re)factor the jacobian matrix (*eg* `factorize` for Eigen (later calls), `klu_refactor` for KLU or `FactorizeMatrix` for NICSLU (later calls)")
        .value("SolverSolve", ErrorType::SolverSolve, "The linear solve failed to solve the linear system J.X = b (*eg* `solve` for Eigen, `klu_solve` for KLU or `Solve` for NICSLU")
        .value("NotInitError", ErrorType::NotInitError, "Attempt to perform some powerflow computation when the linear solver is not initiliazed")
        .value("LicenseError", ErrorType::LicenseError, "Impossible to use the linear solver as the license cannot be found (*eg* unable to locate the `nicslu.lic` file")
        .export_values();

    py::class_<SparseLUSolver>(m, "SparseLUSolver", DocSolver::SparseLUSolver.c_str())
        .def(py::init<>())
        .def("get_J", &SparseLUSolver::get_J_python, DocSolver::get_J_python.c_str())  // (get the jacobian matrix, sparse csc matrix)
        .def("get_Va", &SparseLUSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
        .def("get_Vm", &SparseLUSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
        .def("get_V", &SparseLUSolver::get_V, DocSolver::get_V.c_str()) 
        .def("get_error", &SparseLUSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &SparseLUSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
        .def("reset", &SparseLUSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
        .def("converged", &SparseLUSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
        .def("compute_pf", &SparseLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
        .def("get_timers", &SparseLUSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &SparseLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
    
    py::class_<SparseLUSolverSingleSlack>(m, "SparseLUSolverSingleSlack", DocSolver::SparseLUSolverSingleSlack.c_str())
        .def(py::init<>())
        .def("get_J", &SparseLUSolverSingleSlack::get_J_python, DocSolver::get_J_python.c_str())  // (get the jacobian matrix, sparse csc matrix)
        .def("get_Va", &SparseLUSolverSingleSlack::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
        .def("get_Vm", &SparseLUSolverSingleSlack::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
        .def("get_V", &SparseLUSolverSingleSlack::get_V, DocSolver::get_V.c_str()) 
        .def("get_error", &SparseLUSolverSingleSlack::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &SparseLUSolverSingleSlack::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
        .def("reset", &SparseLUSolverSingleSlack::reset, DocSolver::reset.c_str())  // reset the solver to its original state
        .def("converged", &SparseLUSolverSingleSlack::converged, DocSolver::converged.c_str())  // whether the solver has converged
        .def("compute_pf", &SparseLUSolverSingleSlack::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
        .def("get_timers", &SparseLUSolverSingleSlack::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &SparseLUSolverSingleSlack::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization

    py::class_<DCSolver>(m, "DCSolver", DocSolver::DCSolver.c_str())
        .def(py::init<>())
        .def("get_Va", &DCSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
        .def("get_Vm", &DCSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
        .def("get_V", &DCSolver::get_V, DocSolver::get_V.c_str()) 
        .def("get_error", &DCSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &DCSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
        .def("reset", &DCSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
        .def("converged", &DCSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
        .def("compute_pf", &DCSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
        .def("get_timers", &DCSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &DCSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization

    py::class_<FDPF_XB_SparseLUSolver>(m, "FDPF_XB_SparseLUSolver", DocSolver::FDPF_XB_SparseLUSolver.c_str())
        .def(py::init<>())
        .def("get_Va", &FDPF_XB_SparseLUSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
        .def("get_Vm", &FDPF_XB_SparseLUSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
        .def("get_V", &FDPF_XB_SparseLUSolver::get_V, DocSolver::get_V.c_str()) 
        .def("get_error", &FDPF_XB_SparseLUSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &FDPF_XB_SparseLUSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
        .def("reset", &FDPF_XB_SparseLUSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
        .def("converged", &FDPF_XB_SparseLUSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
        .def("compute_pf", &FDPF_XB_SparseLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
        .def("get_timers", &FDPF_XB_SparseLUSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &FDPF_XB_SparseLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
        .def("debug_get_Bp_python", &FDPF_XB_SparseLUSolver::debug_get_Bp_python, DocGridModel::_internal_do_not_use.c_str())  // perform the newton raphson optimization
        .def("debug_get_Bpp_python", &FDPF_XB_SparseLUSolver::debug_get_Bpp_python, DocGridModel::_internal_do_not_use.c_str());  // perform the newton raphson optimization

    py::class_<FDPF_BX_SparseLUSolver>(m, "FDPF_BX_SparseLUSolver", DocSolver::FDPF_BX_SparseLUSolver.c_str())
        .def(py::init<>())
        .def("get_Va", &FDPF_BX_SparseLUSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
        .def("get_Vm", &FDPF_BX_SparseLUSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
        .def("get_V", &FDPF_BX_SparseLUSolver::get_V, DocSolver::get_V.c_str()) 
        .def("get_error", &FDPF_BX_SparseLUSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &FDPF_BX_SparseLUSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
        .def("reset", &FDPF_BX_SparseLUSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
        .def("converged", &FDPF_BX_SparseLUSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
        .def("compute_pf", &FDPF_BX_SparseLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
        .def("get_timers", &FDPF_BX_SparseLUSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &FDPF_BX_SparseLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
        .def("debug_get_Bp_python", &FDPF_BX_SparseLUSolver::debug_get_Bp_python, DocGridModel::_internal_do_not_use.c_str())  // perform the newton raphson optimization
        .def("debug_get_Bpp_python", &FDPF_BX_SparseLUSolver::debug_get_Bpp_python, DocGridModel::_internal_do_not_use.c_str());  // perform the newton raphson optimization

    #if defined(KLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
        py::class_<KLUSolver>(m, "KLUSolver", DocSolver::KLUSolver.c_str())
            .def(py::init<>())
            .def("get_J", &KLUSolver::get_J_python, DocSolver::get_J_python.c_str())  // (get the jacobian matrix, sparse csc matrix)
            .def("get_Va", &KLUSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &KLUSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &KLUSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &KLUSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &KLUSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &KLUSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &KLUSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &KLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
            .def("get_timers", &KLUSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &KLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str() );  // perform the newton raphson optimization
        
        py::class_<KLUSolverSingleSlack>(m, "KLUSolverSingleSlack", DocSolver::KLUSolverSingleSlack.c_str())
            .def(py::init<>())
            .def("get_J", &KLUSolverSingleSlack::get_J_python, DocSolver::get_J_python.c_str())  // (get the jacobian matrix, sparse csc matrix)
            .def("get_Va", &KLUSolverSingleSlack::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &KLUSolverSingleSlack::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &KLUSolverSingleSlack::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &KLUSolverSingleSlack::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &KLUSolverSingleSlack::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &KLUSolverSingleSlack::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &KLUSolverSingleSlack::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &KLUSolverSingleSlack::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
            .def("get_timers", &KLUSolverSingleSlack::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &KLUSolverSingleSlack::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
        
        py::class_<KLUDCSolver>(m, "KLUDCSolver", DocSolver::KLUDCSolver.c_str())
            .def(py::init<>())
            .def("get_Va", &KLUDCSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &KLUDCSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &KLUDCSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &KLUDCSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &KLUDCSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &KLUDCSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &KLUDCSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &KLUDCSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
            .def("get_timers", &KLUDCSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &KLUDCSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
        
        py::class_<FDPF_XB_KLUSolver>(m, "FDPF_XB_KLUSolver", DocSolver::FDPF_XB_KLUSolver.c_str())
            .def(py::init<>())
            .def("get_Va", &FDPF_XB_KLUSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &FDPF_XB_KLUSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &FDPF_XB_KLUSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &FDPF_XB_KLUSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &FDPF_XB_KLUSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &FDPF_XB_KLUSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &FDPF_XB_KLUSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &FDPF_XB_KLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
            .def("get_timers", &FDPF_XB_KLUSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &FDPF_XB_KLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
                
        py::class_<FDPF_BX_KLUSolver>(m, "FDPF_BX_KLUSolver", DocSolver::FDPF_BX_KLUSolver.c_str())
            .def(py::init<>())
            .def("get_Va", &FDPF_BX_KLUSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &FDPF_BX_KLUSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &FDPF_BX_KLUSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &FDPF_BX_KLUSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &FDPF_BX_KLUSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &FDPF_BX_KLUSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &FDPF_BX_KLUSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &FDPF_BX_KLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
            .def("get_timers", &FDPF_BX_KLUSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &FDPF_BX_KLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
                
    #endif  // KLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

    #if defined(NICSLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
        py::class_<NICSLUSolver>(m, "NICSLUSolver", DocSolver::NICSLUSolver.c_str())
            .def(py::init<>())
            .def("get_J", &NICSLUSolver::get_J_python, DocSolver::get_J_python.c_str())  // (get the jacobian matrix, sparse csc matrix)
            .def("get_Va", &NICSLUSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &NICSLUSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &NICSLUSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &NICSLUSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &NICSLUSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &NICSLUSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &NICSLUSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &NICSLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
            .def("get_timers", &NICSLUSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &NICSLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
        
        py::class_<NICSLUSolverSingleSlack>(m, "NICSLUSolverSingleSlack", DocSolver::NICSLUSolverSingleSlack.c_str())
            .def(py::init<>())
            .def("get_J", &NICSLUSolverSingleSlack::get_J_python, DocSolver::get_J_python.c_str())  // (get the jacobian matrix, sparse csc matrix)
            .def("get_Va", &NICSLUSolverSingleSlack::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &NICSLUSolverSingleSlack::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &NICSLUSolverSingleSlack::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &NICSLUSolverSingleSlack::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &NICSLUSolverSingleSlack::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &NICSLUSolverSingleSlack::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &NICSLUSolverSingleSlack::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &NICSLUSolverSingleSlack::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
            .def("get_timers", &NICSLUSolverSingleSlack::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &NICSLUSolverSingleSlack::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
        
        py::class_<NICSLUDCSolver>(m, "NICSLUDCSolver", DocSolver::NICSLUDCSolver.c_str())
            .def(py::init<>())
            .def("get_Va", &NICSLUDCSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &NICSLUDCSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &NICSLUDCSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &NICSLUDCSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &NICSLUDCSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &NICSLUDCSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &NICSLUDCSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &NICSLUDCSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
            .def("get_timers", &NICSLUDCSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &NICSLUDCSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
        
        py::class_<FDPF_XB_NICSLUSolver>(m, "FDPF_XB_NICSLUSolver", DocSolver::FDPF_XB_NICSLUSolver.c_str())
            .def(py::init<>())
            .def("get_Va", &FDPF_XB_NICSLUSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &FDPF_XB_NICSLUSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &FDPF_XB_NICSLUSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &FDPF_XB_NICSLUSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &FDPF_XB_NICSLUSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &FDPF_XB_NICSLUSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &FDPF_XB_NICSLUSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &FDPF_XB_NICSLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
            .def("get_timers", &FDPF_XB_NICSLUSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &FDPF_XB_NICSLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
        
        py::class_<FDPF_BX_NICSLUSolver>(m, "FDPF_BX_NICSLUSolver", DocSolver::FDPF_BX_NICSLUSolver.c_str())
            .def(py::init<>())
            .def("get_Va", &FDPF_BX_NICSLUSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &FDPF_BX_NICSLUSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &FDPF_BX_NICSLUSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &FDPF_BX_NICSLUSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &FDPF_BX_NICSLUSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &FDPF_BX_NICSLUSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &FDPF_BX_NICSLUSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &FDPF_BX_NICSLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
            .def("get_timers", &FDPF_BX_NICSLUSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &FDPF_BX_NICSLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
    
    #endif  // NICSLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

    #if defined(CKTSO_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
        py::class_<CKTSOSolver>(m, "CKTSOSolver", DocSolver::CKTSOSolver.c_str())
            .def(py::init<>())
            .def("get_J", &CKTSOSolver::get_J_python, DocSolver::get_J_python.c_str())  // (get the jacobian matrix, sparse csc matrix)
            .def("get_Va", &CKTSOSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &CKTSOSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &CKTSOSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &CKTSOSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &CKTSOSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &CKTSOSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &CKTSOSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &CKTSOSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
            .def("get_timers", &CKTSOSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &CKTSOSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
        
        py::class_<CKTSOSolverSingleSlack>(m, "CKTSOSolverSingleSlack", DocSolver::CKTSOSolverSingleSlack.c_str())
            .def(py::init<>())
            .def("get_J", &CKTSOSolverSingleSlack::get_J_python, DocSolver::get_J_python.c_str())  // (get the jacobian matrix, sparse csc matrix)
            .def("get_Va", &CKTSOSolverSingleSlack::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &CKTSOSolverSingleSlack::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &CKTSOSolverSingleSlack::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &CKTSOSolverSingleSlack::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &CKTSOSolverSingleSlack::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &CKTSOSolverSingleSlack::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &CKTSOSolverSingleSlack::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &CKTSOSolverSingleSlack::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
            .def("get_timers", &CKTSOSolverSingleSlack::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &CKTSOSolverSingleSlack::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
        
        py::class_<CKTSODCSolver>(m, "CKTSODCSolver", DocSolver::CKTSODCSolver.c_str())
            .def(py::init<>())
            .def("get_Va", &CKTSODCSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &CKTSODCSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &CKTSODCSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &CKTSODCSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &CKTSODCSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &CKTSODCSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &CKTSODCSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &CKTSODCSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // perform the newton raphson optimization
            .def("get_timers", &CKTSODCSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &CKTSODCSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization

        py::class_<FDPF_XB_CKTSOSolver>(m, "FDPF_XB_CKTSOSolver", DocSolver::FDPF_XB_CKTSOSolver.c_str())
            .def(py::init<>())
            .def("get_Va", &FDPF_XB_CKTSOSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &FDPF_XB_CKTSOSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &FDPF_XB_CKTSOSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &FDPF_XB_CKTSOSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &FDPF_XB_CKTSOSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &FDPF_XB_CKTSOSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &FDPF_XB_CKTSOSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &FDPF_XB_CKTSOSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
            .def("get_timers", &FDPF_XB_CKTSOSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &FDPF_XB_CKTSOSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
    
        py::class_<FDPF_BX_CKTSOSolver>(m, "FDPF_BX_CKTSOSolver", DocSolver::FDPF_BX_CKTSOSolver.c_str())
            .def(py::init<>())
            .def("get_Va", &FDPF_BX_CKTSOSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
            .def("get_Vm", &FDPF_BX_CKTSOSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
            .def("get_V", &FDPF_BX_CKTSOSolver::get_V, DocSolver::get_V.c_str()) 
            .def("get_error", &FDPF_BX_CKTSOSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
            .def("get_nb_iter", &FDPF_BX_CKTSOSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
            .def("reset", &FDPF_BX_CKTSOSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
            .def("converged", &FDPF_BX_CKTSOSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
            .def("compute_pf", &FDPF_BX_CKTSOSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
            .def("get_timers", &FDPF_BX_CKTSOSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
            .def("solve", &FDPF_BX_CKTSOSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization
    
    #endif  // CKTSO_SOLVER_AVAILABLE (or _READ_THE_DOCS)

    py::class_<GaussSeidelAlgo>(m, "GaussSeidelSolver", DocSolver::GaussSeidelSolver.c_str())
        .def(py::init<>())
        .def("get_Va", &GaussSeidelAlgo::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
        .def("get_Vm", &GaussSeidelAlgo::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
        .def("get_V", &GaussSeidelAlgo::get_V, DocSolver::get_V.c_str()) 
        .def("get_error", &GaussSeidelAlgo::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &GaussSeidelAlgo::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
        .def("reset", &GaussSeidelAlgo::reset, DocSolver::reset.c_str())  // reset the solver to its original state
        .def("converged", &GaussSeidelAlgo::converged, DocSolver::converged.c_str())  // whether the solver has converged
        .def("compute_pf", &GaussSeidelAlgo::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
        .def("get_timers", &GaussSeidelAlgo::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &GaussSeidelAlgo::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization

    py::class_<GaussSeidelSynchAlgo>(m, "GaussSeidelSynchSolver", DocSolver::GaussSeidelSynchSolver.c_str())
        .def(py::init<>())
        .def("get_Va", &GaussSeidelSynchAlgo::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
        .def("get_Vm", &GaussSeidelSynchAlgo::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
        .def("get_V", &GaussSeidelSynchAlgo::get_V, DocSolver::get_V.c_str()) 
        .def("get_error", &GaussSeidelSynchAlgo::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &GaussSeidelSynchAlgo::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
        .def("reset", &GaussSeidelSynchAlgo::reset, DocSolver::reset.c_str())  // reset the solver to its original state
        .def("converged", &GaussSeidelSynchAlgo::converged, DocSolver::converged.c_str())  // whether the solver has converged
        .def("compute_pf", &GaussSeidelSynchAlgo::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
        .def("get_timers", &GaussSeidelSynchAlgo::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &GaussSeidelSynchAlgo::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization

    // Only "const" method are exported
    // it is so that i cannot modify the internal solver of a gridmodel python side
    py::class_<ChooseSolver>(m, "AnySolver", DocSolver::AnySolver.c_str())
        .def(py::init<>())
        .def("get_type", &ChooseSolver::get_type, DocSolver::get_type.c_str())
        // .def("change_solver", &ChooseSolver::change_solver)
        // .def("reset", &ChooseSolver::reset)
        // .def("compute_pf", &ChooseSolver::compute_pf, py::call_guard<py::gil_scoped_release>())  // compute the powerflow
        // .def("solve", &ChooseSolver::compute_pf, py::call_guard<py::gil_scoped_release>() )
        .def("get_Va", &ChooseSolver::get_Va, DocSolver::get_Va.c_str())  
        .def("get_Vm", &ChooseSolver::get_Vm, DocSolver::get_Vm.c_str()) 
        .def("get_V", &ChooseSolver::get_V, DocSolver::get_V.c_str()) 
        .def("get_J", &ChooseSolver::get_J_python, DocSolver::chooseSolver_get_J_python.c_str()) 
        .def("get_error", &ChooseSolver::get_error, DocSolver::get_V.c_str()) 
        .def("get_nb_iter", &ChooseSolver::get_nb_iter, DocSolver::get_nb_iter.c_str()) 
        .def("converged", &ChooseSolver::converged, DocSolver::converged.c_str()) 
        .def("get_computation_time", &ChooseSolver::get_computation_time, DocSolver::get_computation_time.c_str())
        .def("get_timers", &ChooseSolver::get_timers, "TODO")
        .def("get_timers_jacobian", &ChooseSolver::get_timers_jacobian, "TODO")
        .def("get_fdpf_xb_lu", &ChooseSolver::get_fdpf_xb_lu, py::return_value_policy::reference, DocGridModel::_internal_do_not_use.c_str())  // TODO this for all solver !
        .def("get_fdpf_bx_lu", &ChooseSolver::get_fdpf_bx_lu, py::return_value_policy::reference, DocGridModel::_internal_do_not_use.c_str());

    // iterator for generators
    py::class_<GeneratorContainer>(m, "GeneratorContainer", DocIterator::GeneratorContainer.c_str())
        .def("__len__", [](const GeneratorContainer & data) { return data.nb(); })
        .def("__getitem__", [](const GeneratorContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const GeneratorContainer & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<GeneratorContainer::GenInfo>(m, "GenInfo", DocIterator::GenInfo.c_str())
        .def_readonly("id", &GeneratorContainer::GenInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &GeneratorContainer::GenInfo::name, DocIterator::name.c_str())
        .def_readonly("connected", &GeneratorContainer::GenInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &GeneratorContainer::GenInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("is_slack", &GeneratorContainer::GenInfo::is_slack, DocIterator::is_slack.c_str())
        .def_readonly("slack_weight", &GeneratorContainer::GenInfo::slack_weight, DocIterator::slack_weight.c_str())
        .def_readonly("voltage_regulator_on", &GeneratorContainer::GenInfo::voltage_regulator_on, "TODO")
        .def_readonly("target_p_mw", &GeneratorContainer::GenInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_vm_pu", &GeneratorContainer::GenInfo::target_vm_pu, DocIterator::target_vm_pu.c_str())
        .def_readonly("target_q_mvar", &GeneratorContainer::GenInfo::target_q_mvar, "TODO")
        .def_readonly("min_q_mvar", &GeneratorContainer::GenInfo::min_q_mvar, DocIterator::min_q_mvar.c_str())
        .def_readonly("max_q_mvar", &GeneratorContainer::GenInfo::max_q_mvar, DocIterator::max_q_mvar.c_str())
        .def_readonly("has_res", &GeneratorContainer::GenInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &GeneratorContainer::GenInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &GeneratorContainer::GenInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &GeneratorContainer::GenInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &GeneratorContainer::GenInfo::res_v_kv, DocIterator::res_v_kv.c_str());

    // iterator for sgens
    py::class_<SGenContainer>(m, "SGenContainer", DocIterator::SGenContainer.c_str())
        .def("__len__", [](const SGenContainer & data) { return data.nb(); })
        .def("__getitem__", [](const SGenContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const SGenContainer & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<SGenContainer::SGenInfo>(m, "SGenInfo", DocIterator::SGenInfo.c_str())
        .def_readonly("id", &SGenContainer::SGenInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &SGenContainer::SGenInfo::name, DocIterator::name.c_str())
        .def_readonly("connected", &SGenContainer::SGenInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &SGenContainer::SGenInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("min_q_mvar", &SGenContainer::SGenInfo::min_q_mvar, DocIterator::min_q_mvar.c_str())
        .def_readonly("max_q_mvar", &SGenContainer::SGenInfo::max_q_mvar, DocIterator::max_q_mvar.c_str())
        .def_readonly("min_p_mw", &SGenContainer::SGenInfo::min_p_mw, DocIterator::min_p_mw.c_str())
        .def_readonly("max_p_mw", &SGenContainer::SGenInfo::max_p_mw, DocIterator::max_p_mw.c_str())
        .def_readonly("target_p_mw", &SGenContainer::SGenInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &SGenContainer::SGenInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &SGenContainer::SGenInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &SGenContainer::SGenInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &SGenContainer::SGenInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &SGenContainer::SGenInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &SGenContainer::SGenInfo::res_v_kv, DocIterator::res_v_kv.c_str());

    // iterator for loads (and storage units)
    py::class_<LoadContainer>(m, "LoadContainer", DocIterator::LoadContainer.c_str())
        .def("__len__", [](const LoadContainer & data) { return data.nb(); })
        .def("__getitem__", [](const LoadContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const LoadContainer & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<LoadContainer::LoadInfo>(m, "LoadInfo", DocIterator::LoadInfo.c_str())
        .def_readonly("id", &LoadContainer::LoadInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &LoadContainer::LoadInfo::name, DocIterator::name.c_str())
        .def_readonly("connected", &LoadContainer::LoadInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &LoadContainer::LoadInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("target_p_mw", &LoadContainer::LoadInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &LoadContainer::LoadInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &LoadContainer::LoadInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &LoadContainer::LoadInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &LoadContainer::LoadInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &LoadContainer::LoadInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &LoadContainer::LoadInfo::res_v_kv, DocIterator::res_v_kv.c_str());

    // iterator for shunts
    py::class_<ShuntContainer>(m, "ShuntContainer", DocIterator::ShuntContainer.c_str())
        .def("__len__", [](const ShuntContainer & data) { return data.nb(); })
        .def("__getitem__", [](const ShuntContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const ShuntContainer & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<ShuntContainer::ShuntInfo>(m, "ShuntInfo", DocIterator::ShuntInfo.c_str())
        .def_readonly("id", &ShuntContainer::ShuntInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &ShuntContainer::ShuntInfo::name, DocIterator::name.c_str())
        .def_readonly("connected", &ShuntContainer::ShuntInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &ShuntContainer::ShuntInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("target_p_mw", &ShuntContainer::ShuntInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &ShuntContainer::ShuntInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &ShuntContainer::ShuntInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &ShuntContainer::ShuntInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &ShuntContainer::ShuntInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &ShuntContainer::ShuntInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &ShuntContainer::ShuntInfo::res_v_kv, DocIterator::res_v_kv.c_str());

    // iterator for trafos
    py::class_<TrafoContainer>(m, "TrafoContainer", DocIterator::TrafoContainer.c_str())
        .def("__len__", [](const TrafoContainer & data) { return data.nb(); })
        .def("__getitem__", [](const TrafoContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const TrafoContainer & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<TrafoContainer::TrafoInfo>(m, "TrafoInfo", DocIterator::TrafoInfo.c_str())
        .def_readonly("id", &TrafoContainer::TrafoInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &TrafoContainer::TrafoInfo::name, DocIterator::name.c_str())
        .def_readonly("connected", &TrafoContainer::TrafoInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_hv_id", &TrafoContainer::TrafoInfo::bus_hv_id, DocIterator::bus_hv_id.c_str())
        .def_readonly("bus_lv_id", &TrafoContainer::TrafoInfo::bus_lv_id, DocIterator::bus_lv_id.c_str())
        .def_readonly("r_pu", &TrafoContainer::TrafoInfo::r_pu, DocIterator::r_pu.c_str())
        .def_readonly("x_pu", &TrafoContainer::TrafoInfo::x_pu, DocIterator::x_pu.c_str())
        .def_readonly("h_pu", &TrafoContainer::TrafoInfo::h_pu, DocIterator::h_pu.c_str())
        .def_readonly("is_tap_hv_side", &TrafoContainer::TrafoInfo::is_tap_hv_side, DocIterator::is_tap_hv_side.c_str())
        .def_readonly("ratio", &TrafoContainer::TrafoInfo::ratio, DocIterator::ratio.c_str())
        .def_readonly("shift_rad", &TrafoContainer::TrafoInfo::shift_rad, DocIterator::shift_rad.c_str())
        .def_readonly("has_res", &TrafoContainer::TrafoInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_hv_mw", &TrafoContainer::TrafoInfo::res_p_hv_mw, DocIterator::res_p_hv_mw.c_str())
        .def_readonly("res_q_hv_mvar", &TrafoContainer::TrafoInfo::res_q_hv_mvar, DocIterator::res_q_hv_mvar.c_str())
        .def_readonly("res_v_hv_kv", &TrafoContainer::TrafoInfo::res_v_hv_kv, DocIterator::res_v_hv_kv.c_str())
        .def_readonly("res_a_hv_ka", &TrafoContainer::TrafoInfo::res_a_hv_ka, DocIterator::res_a_hv_ka.c_str())
        .def_readonly("res_p_lv_mw", &TrafoContainer::TrafoInfo::res_p_lv_mw, DocIterator::res_p_lv_mw.c_str())
        .def_readonly("res_q_lv_mvar", &TrafoContainer::TrafoInfo::res_q_lv_mvar, DocIterator::res_q_lv_mvar.c_str())
        .def_readonly("res_v_lv_kv", &TrafoContainer::TrafoInfo::res_v_lv_kv, DocIterator::res_v_lv_kv.c_str())
        .def_readonly("res_a_lv_ka", &TrafoContainer::TrafoInfo::res_a_lv_ka, DocIterator::res_a_lv_ka.c_str())
        .def_readonly("res_theta_hv_deg", &TrafoContainer::TrafoInfo::res_theta_hv_deg, DocIterator::res_theta_hv_deg.c_str())
        .def_readonly("res_theta_lv_deg", &TrafoContainer::TrafoInfo::res_theta_lv_deg, DocIterator::res_theta_lv_deg.c_str());

    // iterator for trafos
    py::class_<LineContainer>(m, "LineContainer", DocIterator::LineContainer.c_str())
        .def("__len__", [](const LineContainer & data) { return data.nb(); })
        .def("__getitem__", [](const LineContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const LineContainer & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<LineContainer::LineInfo>(m, "LineInfo", DocIterator::LineInfo.c_str())
        .def_readonly("id", &LineContainer::LineInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &LineContainer::LineInfo::name, DocIterator::name.c_str())
        .def_readonly("connected", &LineContainer::LineInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_or_id", &LineContainer::LineInfo::bus_or_id, DocIterator::bus_or_id.c_str())
        .def_readonly("bus_ex_id", &LineContainer::LineInfo::bus_ex_id, DocIterator::bus_ex_id.c_str())
        .def_readonly("r_pu", &LineContainer::LineInfo::r_pu, DocIterator::r_pu.c_str())
        .def_readonly("x_pu", &LineContainer::LineInfo::x_pu, DocIterator::x_pu.c_str())
        .def_readonly("h_pu", &LineContainer::LineInfo::h_pu, DocIterator::x_pu.c_str())
        .def_readonly("h_or_pu", &LineContainer::LineInfo::h_or_pu, DocIterator::h_pu.c_str())
        .def_readonly("h_ex_pu", &LineContainer::LineInfo::h_ex_pu, DocIterator::h_pu.c_str())
        .def_readonly("has_res", &LineContainer::LineInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_or_mw", &LineContainer::LineInfo::res_p_or_mw, DocIterator::res_p_or_mw.c_str())
        .def_readonly("res_q_or_mvar", &LineContainer::LineInfo::res_q_or_mvar, DocIterator::res_q_or_mvar.c_str())
        .def_readonly("res_v_or_kv", &LineContainer::LineInfo::res_v_or_kv, DocIterator::res_v_or_kv.c_str())
        .def_readonly("res_a_or_ka", &LineContainer::LineInfo::res_a_or_ka, DocIterator::res_a_or_ka.c_str())
        .def_readonly("res_p_ex_mw", &LineContainer::LineInfo::res_p_ex_mw, DocIterator::res_p_ex_mw.c_str())
        .def_readonly("res_q_ex_mvar", &LineContainer::LineInfo::res_q_ex_mvar, DocIterator::res_q_ex_mvar.c_str())
        .def_readonly("res_v_ex_kv", &LineContainer::LineInfo::res_v_ex_kv, DocIterator::res_v_ex_kv.c_str())
        .def_readonly("res_a_ex_ka", &LineContainer::LineInfo::res_a_ex_ka, DocIterator::res_a_ex_ka.c_str())
        .def_readonly("res_theta_or_deg", &LineContainer::LineInfo::res_theta_or_deg, DocIterator::res_theta_or_deg.c_str())
        .def_readonly("res_theta_ex_deg", &LineContainer::LineInfo::res_theta_ex_deg, DocIterator::res_theta_ex_deg.c_str());

    // iterator for dc lines
    py::class_<DCLineContainer>(m, "DCLineContainer", DocIterator::DCLineContainer.c_str())
        .def("__len__", [](const DCLineContainer & data) { return data.nb(); })
        .def("__getitem__", [](const DCLineContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const DCLineContainer & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<DCLineContainer::DCLineInfo>(m, "DCLineInfo", DocIterator::DCLineInfo.c_str())
        .def_readonly("id", &DCLineContainer::DCLineInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &DCLineContainer::DCLineInfo::name, DocIterator::name.c_str())
        .def_readonly("connected", &DCLineContainer::DCLineInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_or_id", &DCLineContainer::DCLineInfo::bus_or_id, DocIterator::bus_or_id.c_str())
        .def_readonly("bus_ex_id", &DCLineContainer::DCLineInfo::bus_ex_id, DocIterator::bus_ex_id.c_str())
        .def_readonly("target_p_or_mw", &DCLineContainer::DCLineInfo::target_p_or_mw, DocIterator::target_p_or_mw.c_str())
        .def_readonly("target_vm_or_pu", &DCLineContainer::DCLineInfo::target_vm_or_pu, DocIterator::target_vm_or_pu.c_str())
        .def_readonly("target_vm_ex_pu", &DCLineContainer::DCLineInfo::target_vm_ex_pu, DocIterator::target_vm_ex_pu.c_str())
        .def_readonly("loss_pct", &DCLineContainer::DCLineInfo::loss_pct, DocIterator::loss_pct.c_str())
        .def_readonly("loss_mw", &DCLineContainer::DCLineInfo::loss_mw, DocIterator::loss_mw.c_str())
        .def_readonly("gen_or", &DCLineContainer::DCLineInfo::gen_or, DocIterator::gen_or.c_str())
        .def_readonly("gen_ex", &DCLineContainer::DCLineInfo::gen_ex, DocIterator::gen_ex.c_str())
        .def_readonly("has_res", &DCLineContainer::DCLineInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_or_mw", &DCLineContainer::DCLineInfo::res_p_or_mw, DocIterator::res_p_or_mw_dcline.c_str())
        .def_readonly("res_p_ex_mw", &DCLineContainer::DCLineInfo::res_p_ex_mw, DocIterator::res_p_ex_mw_dcline.c_str())
        .def_readonly("res_q_or_mvar", &DCLineContainer::DCLineInfo::res_q_or_mvar, DocIterator::res_q_or_mvar_dcline.c_str())
        .def_readonly("res_q_ex_mvar", &DCLineContainer::DCLineInfo::res_q_ex_mvar, DocIterator::res_q_ex_mvar_dcline.c_str())
        .def_readonly("res_v_or_kv", &DCLineContainer::DCLineInfo::res_v_or_kv, DocIterator::res_v_or_kv_dcline.c_str())
        .def_readonly("res_v_ex_kv", &DCLineContainer::DCLineInfo::res_v_ex_kv, DocIterator::res_v_ex_kv_dcline.c_str())
        .def_readonly("res_theta_or_deg", &DCLineContainer::DCLineInfo::res_theta_or_deg, DocIterator::res_theta_or_deg_dcline.c_str())
        .def_readonly("res_theta_ex_deg", &DCLineContainer::DCLineInfo::res_theta_ex_deg, DocIterator::res_theta_ex_deg_dcline.c_str())
        ;

    // converters
    py::class_<PandaPowerConverter>(m, "PandaPowerConverter")
        .def(py::init<>())
        .def("set_f_hz", &PandaPowerConverter::set_f_hz)
        .def("set_sn_mva", &PandaPowerConverter::set_sn_mva)
        .def("get_line_param", &PandaPowerConverter::get_line_param)
        .def("get_trafo_param", &PandaPowerConverter::get_trafo_param);

    py::class_<SolverControl>(m, "SolverControl", "TODO")
        .def(py::init<>())
        .def("has_dimension_changed", &SolverControl::has_dimension_changed, "TODO")
        .def("has_pv_changed", &SolverControl::has_pv_changed, "TODO")
        .def("has_pq_changed", &SolverControl::has_pq_changed, "TODO")
        .def("has_slack_participate_changed", &SolverControl::has_slack_participate_changed, "TODO")
        .def("need_reset_solver", &SolverControl::need_reset_solver, "TODO")
        .def("need_recompute_sbus", &SolverControl::need_recompute_sbus, "TODO")
        .def("need_recompute_ybus", &SolverControl::need_recompute_ybus, "TODO")
        .def("ybus_change_sparsity_pattern", &SolverControl::ybus_change_sparsity_pattern, "TODO")
        .def("has_slack_weight_changed", &SolverControl::has_slack_weight_changed, "TODO")
        .def("has_v_changed", &SolverControl::has_v_changed, "TODO")
        .def("has_ybus_some_coeffs_zero", &SolverControl::has_ybus_some_coeffs_zero, "TODO")
        ;

    py::class_<GridModel>(m, "GridModel", DocGridModel::GridModel.c_str())
        .def(py::init<>())
        .def("copy", &GridModel::copy)
        .def_property("_ls_to_orig", &GridModel::get_ls_to_orig, &GridModel::set_ls_to_orig, "remember the conversion from bus index in lightsim2grid to bus index in original file format (*eg* pandapower of pypowsybl).")
        .def_property("_orig_to_ls", &GridModel::get_orig_to_ls, &GridModel::set_orig_to_ls, "remember the conversion from bus index in original file format (*eg* pandapower of pypowsybl) to bus index in lightsim2grid.")
        .def_property("_max_nb_bus_per_sub",
                      &GridModel::get_max_nb_bus_per_sub,
                      &GridModel::set_max_nb_bus_per_sub,
                      "do not modify it after loading !")
        .def_property_readonly("timer_last_ac_pf", &GridModel::timer_last_ac_pf, "TODO")
        .def_property_readonly("timer_last_dc_pf", &GridModel::timer_last_dc_pf, "TODO")

        // pickle
        .def(py::pickle(
                        [](const GridModel &gm) { // __getstate__
                            // Return a tuple that fully encodes the state of the object
                            return py::make_tuple(gm.get_state());
                        },
                        [](py::tuple py_state) { // __setstate__
                            if (py_state.size() != 1){
                                std::cout << "GridModel.__setstate__ : state size " << py_state.size() << std::endl;
                                throw std::runtime_error("Invalid state size when loading GridModel.__setstate__");
                                }
                            // Create a new C++ instance
                            GridModel gm = GridModel();
                            // TODO check the size of the input tuple!

                            // now set the status
                            GridModel::StateRes state = py_state[0].cast<GridModel::StateRes>();
                            gm.set_state(state);
                            return gm;
        }))

        // general parameters
        // solver control
        .def("change_solver", &GridModel::change_solver, DocGridModel::change_solver.c_str())
        .def("available_solvers", &GridModel::available_solvers, DocGridModel::available_solvers.c_str())  // retrieve the solver available for your installation
        .def("get_computation_time", &GridModel::get_computation_time, DocGridModel::get_computation_time.c_str())  // get the computation time spent in the solver
        .def("get_dc_computation_time", &GridModel::get_dc_computation_time, DocGridModel::get_dc_computation_time.c_str())  // get the computation time spent in the solver
        .def("get_solver_type", &GridModel::get_solver_type, DocGridModel::get_solver_type.c_str())  // get the type of solver used
        .def("get_dc_solver_type", &GridModel::get_dc_solver_type, DocGridModel::get_dc_solver_type.c_str())  // get the type of solver used
        .def("get_solver", &GridModel::get_solver, py::return_value_policy::reference, DocGridModel::get_solver.c_str())  // get the solver (AnySolver type python side) used
        .def("get_dc_solver", &GridModel::get_dc_solver, py::return_value_policy::reference, DocGridModel::get_dc_solver.c_str())  // get the solver (AnySolver type python side) used

        // init the grid
        .def("init_bus", &GridModel::init_bus, DocGridModel::_internal_do_not_use.c_str())
        .def("init_bus_status", &GridModel::init_bus_status, DocGridModel::_internal_do_not_use.c_str())
        .def("set_init_vm_pu", &GridModel::set_init_vm_pu, DocGridModel::_internal_do_not_use.c_str())  // TODO use python "property" for that
        .def("get_init_vm_pu", &GridModel::get_init_vm_pu, DocGridModel::_internal_do_not_use.c_str())
        .def("set_sn_mva", &GridModel::set_sn_mva, DocGridModel::_internal_do_not_use.c_str())   // TODO use python "property" for that
        .def("get_sn_mva", &GridModel::get_sn_mva, DocGridModel::_internal_do_not_use.c_str())

        // init its elements
        .def("init_powerlines", &GridModel::init_powerlines, DocGridModel::_internal_do_not_use.c_str())  // TODO code the possibility to add / remove a powerline after creation
        .def("init_powerlines_full", &GridModel::init_powerlines_full, DocGridModel::_internal_do_not_use.c_str())  // TODO code the possibility to add / remove a powerline after creation
        .def("init_shunt", &GridModel::init_shunt, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_trafo", &GridModel::init_trafo, DocGridModel::_internal_do_not_use.c_str())  // same 
        .def("init_generators", &GridModel::init_generators, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_generators_full", &GridModel::init_generators_full, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_loads", &GridModel::init_loads, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_storages", &GridModel::init_storages, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_sgens", &GridModel::init_sgens, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_dclines", &GridModel::init_dclines, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("add_gen_slackbus", &GridModel::add_gen_slackbus, DocGridModel::_internal_do_not_use.c_str()) // same
        .def("remove_gen_slackbus", &GridModel::remove_gen_slackbus, DocGridModel::_internal_do_not_use.c_str())  // same

        // inspect the grid
        .def("get_lines", &GridModel::get_lines, DocGridModel::get_lines.c_str())
        .def("get_dclines", &GridModel::get_dclines, DocGridModel::get_dclines.c_str())
        .def("get_trafos", &GridModel::get_trafos, DocGridModel::get_trafos.c_str())
        .def("get_generators", &GridModel::get_generators, DocGridModel::get_generators.c_str())
        .def("get_static_generators", &GridModel::get_static_generators, DocGridModel::get_static_generators.c_str())
        .def("get_shunts", &GridModel::get_shunts, DocGridModel::get_shunts.c_str())
        .def("get_storages", &GridModel::get_storages, DocGridModel::get_storages.c_str())
        .def("get_loads", &GridModel::get_loads, DocGridModel::get_loads.c_str())
        .def("get_bus_vn_kv", &GridModel::get_bus_vn_kv, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_status", &GridModel::get_bus_status, DocGridModel::_internal_do_not_use.c_str())

        // modify the grid
        .def("turnedoff_no_pv", &GridModel::turnedoff_no_pv, "Turned off (or generators with p = 0) generators will not be pv buses, they will not maintain voltage")
        .def("turnedoff_pv", &GridModel::turnedoff_pv, "Turned off (or generators with p = 0) generators will be pv buses, they will maintain voltage (default)")
        .def("get_turnedoff_gen_pv", &GridModel::get_turnedoff_gen_pv, "TODO")
        .def("update_slack_weights", &GridModel::update_slack_weights, "TODO")
        .def("assign_slack_to_most_connected", &GridModel::assign_slack_to_most_connected, "TODO")
        .def("consider_only_main_component", &GridModel::consider_only_main_component, "TODO and TODO DC LINE: one side might be in the connected comp and not the other !")
        
        // names
        .def("set_line_names", &GridModel::set_line_names, "TODO")
        .def("set_dcline_names", &GridModel::set_dcline_names, "TODO")
        .def("set_trafo_names", &GridModel::set_trafo_names, "TODO")
        .def("set_gen_names", &GridModel::set_gen_names, "TODO")
        .def("set_load_names", &GridModel::set_load_names, "TODO")
        .def("set_storage_names", &GridModel::set_storage_names,  "TODO")
        .def("set_sgen_names", &GridModel::set_sgen_names, "TODO")
        .def("set_shunt_names", &GridModel::set_shunt_names, "TODO")

        .def("deactivate_bus", &GridModel::deactivate_bus, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_bus", &GridModel::reactivate_bus, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_powerline", &GridModel::deactivate_powerline, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_powerline", &GridModel::reactivate_powerline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_powerline_or", &GridModel::change_bus_powerline_or, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_powerline_ex", &GridModel::change_bus_powerline_ex, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_powerline_or", &GridModel::get_bus_powerline_or, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_powerline_ex", &GridModel::get_bus_powerline_ex, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_trafo", &GridModel::deactivate_trafo, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_trafo", &GridModel::reactivate_trafo, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_trafo_hv", &GridModel::change_bus_trafo_hv, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_trafo_lv", &GridModel::change_bus_trafo_lv, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_trafo_hv", &GridModel::get_bus_trafo_hv, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_trafo_lv", &GridModel::get_bus_trafo_lv, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_load", &GridModel::deactivate_load, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_load", &GridModel::reactivate_load, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_load", &GridModel::change_bus_load, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_load", &GridModel::get_bus_load, DocGridModel::_internal_do_not_use.c_str())
        .def("change_p_load", &GridModel::change_p_load, DocGridModel::_internal_do_not_use.c_str())
        .def("change_q_load", &GridModel::change_q_load, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_gen", &GridModel::deactivate_gen, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_gen", &GridModel::reactivate_gen, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_gen", &GridModel::change_bus_gen, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_gen", &GridModel::get_bus_gen, DocGridModel::_internal_do_not_use.c_str())
        .def("change_p_gen", &GridModel::change_p_gen, DocGridModel::_internal_do_not_use.c_str())
        .def("change_v_gen", &GridModel::change_v_gen, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_shunt", &GridModel::deactivate_shunt, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_shunt", &GridModel::reactivate_shunt, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_shunt", &GridModel::change_bus_shunt, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_shunt", &GridModel::get_bus_shunt, DocGridModel::_internal_do_not_use.c_str())
        .def("change_p_shunt", &GridModel::change_p_shunt, DocGridModel::_internal_do_not_use.c_str())
        .def("change_q_shunt", &GridModel::change_q_shunt, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_sgen", &GridModel::deactivate_sgen, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_sgen", &GridModel::reactivate_sgen, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_sgen", &GridModel::change_bus_sgen, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_sgen", &GridModel::get_bus_sgen, DocGridModel::_internal_do_not_use.c_str())
        .def("change_p_sgen", &GridModel::change_p_sgen, DocGridModel::_internal_do_not_use.c_str())
        .def("change_q_sgen", &GridModel::change_q_sgen, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_storage", &GridModel::deactivate_storage, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_storage", &GridModel::reactivate_storage, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_storage", &GridModel::change_bus_storage, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_storage", &GridModel::get_bus_storage, DocGridModel::_internal_do_not_use.c_str())
        .def("change_p_storage", &GridModel::change_p_storage, DocGridModel::_internal_do_not_use.c_str())
        .def("change_q_storage", &GridModel::change_q_storage, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_dcline", &GridModel::deactivate_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_dcline", &GridModel::reactivate_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_p_dcline", &GridModel::change_p_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_v_or_dcline", &GridModel::change_v_or_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_v_ex_dcline", &GridModel::change_v_ex_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_dcline_or", &GridModel::change_bus_dcline_or, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_dcline_ex", &GridModel::change_bus_dcline_ex, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_dcline_or", &GridModel::get_bus_dcline_or, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_dcline_ex", &GridModel::get_bus_dcline_ex, DocGridModel::_internal_do_not_use.c_str())

        // get back the results
        .def("get_V", &GridModel::get_V, DocGridModel::get_V.c_str())
        .def("get_Va", &GridModel::get_Va, DocGridModel::get_Va.c_str())
        .def("get_Vm", &GridModel::get_Vm, DocGridModel::get_Vm.c_str())
        .def("get_J", &GridModel::get_J_python, DocGridModel::get_J_python.c_str())
        .def("id_me_to_ac_solver", &GridModel::id_me_to_ac_solver, DocGridModel::id_me_to_ac_solver.c_str())
        .def("id_ac_solver_to_me", &GridModel::id_ac_solver_to_me, DocGridModel::id_ac_solver_to_me.c_str())
        .def("id_me_to_dc_solver", &GridModel::id_me_to_dc_solver, DocGridModel::id_me_to_dc_solver.c_str())
        .def("id_dc_solver_to_me", &GridModel::id_dc_solver_to_me, DocGridModel::id_dc_solver_to_me.c_str())
        .def("total_bus", &GridModel::total_bus, DocGridModel::total_bus.c_str())
        .def("nb_bus", &GridModel::nb_bus, DocGridModel::nb_bus.c_str())
        .def("get_pv", &GridModel::get_pv, DocGridModel::get_pv.c_str())
        .def("get_pq", &GridModel::get_pq, DocGridModel::get_pq.c_str())
        .def("get_slack_ids", &GridModel::get_slack_ids, DocGridModel::get_slack_ids.c_str())
        .def("get_slack_weights", &GridModel::get_slack_weights, DocGridModel::get_slack_weights.c_str())
        .def("get_Ybus", &GridModel::get_Ybus, DocGridModel::get_Ybus.c_str())
        .def("get_dcYbus", &GridModel::get_dcYbus, DocGridModel::get_dcYbus.c_str())
        .def("get_Sbus", &GridModel::get_Sbus, DocGridModel::get_Sbus.c_str())
        .def("get_dcSbus", &GridModel::get_dcSbus, DocGridModel::_internal_do_not_use.c_str())

        .def("check_solution", &GridModel::check_solution, DocGridModel::check_solution.c_str())

        // TODO optimize that for speed, results are copied apparently
        .def("get_loads_res", &GridModel::get_loads_res, DocGridModel::_internal_do_not_use.c_str())
        .def("get_loads_status", &GridModel::get_loads_status, DocGridModel::_internal_do_not_use.c_str())
        .def("get_shunts_res", &GridModel::get_shunts_res, DocGridModel::_internal_do_not_use.c_str())
        .def("get_shunts_status", &GridModel::get_shunts_status, DocGridModel::_internal_do_not_use.c_str())
        .def("get_gen_res", &GridModel::get_gen_res, DocGridModel::_internal_do_not_use.c_str())
        .def("get_gen_status", &GridModel::get_gen_status, DocGridModel::_internal_do_not_use.c_str())
        .def("get_lineor_res", &GridModel::get_lineor_res, DocGridModel::_internal_do_not_use.c_str())
        .def("get_lineex_res", &GridModel::get_lineex_res, DocGridModel::_internal_do_not_use.c_str())
        .def("get_lines_status", &GridModel::get_lines_status, DocGridModel::_internal_do_not_use.c_str())
        .def("get_trafohv_res", &GridModel::get_trafohv_res, DocGridModel::_internal_do_not_use.c_str())
        .def("get_trafolv_res", &GridModel::get_trafolv_res, DocGridModel::_internal_do_not_use.c_str())
        .def("get_trafo_status", &GridModel::get_trafo_status, DocGridModel::_internal_do_not_use.c_str())
        .def("get_storages_res", &GridModel::get_storages_res, DocGridModel::_internal_do_not_use.c_str())
        .def("get_storages_status", &GridModel::get_storages_status, DocGridModel::_internal_do_not_use.c_str())
        .def("get_sgens_res", &GridModel::get_sgens_res, DocGridModel::_internal_do_not_use.c_str())
        .def("get_sgens_status", &GridModel::get_sgens_status, DocGridModel::_internal_do_not_use.c_str())

        .def("get_gen_theta", &GridModel::get_gen_theta, DocGridModel::_internal_do_not_use.c_str())
        .def("get_load_theta", &GridModel::get_load_theta, DocGridModel::_internal_do_not_use.c_str())
        .def("get_shunt_theta", &GridModel::get_shunt_theta, DocGridModel::_internal_do_not_use.c_str())
        .def("get_storage_theta", &GridModel::get_storage_theta, DocGridModel::_internal_do_not_use.c_str())
        .def("get_lineor_theta", &GridModel::get_lineor_theta, DocGridModel::_internal_do_not_use.c_str())
        .def("get_lineex_theta", &GridModel::get_lineex_theta, DocGridModel::_internal_do_not_use.c_str())
        .def("get_trafohv_theta", &GridModel::get_trafohv_theta, DocGridModel::_internal_do_not_use.c_str())
        .def("get_trafolv_theta", &GridModel::get_trafolv_theta, DocGridModel::_internal_do_not_use.c_str())

        .def("get_all_shunt_buses", &GridModel::get_all_shunt_buses, DocGridModel::_internal_do_not_use.c_str())
        .def("get_loads_res_full", &GridModel::get_loads_res_full, DocGridModel::_internal_do_not_use.c_str())
        .def("get_shunts_res_full", &GridModel::get_shunts_res_full, DocGridModel::_internal_do_not_use.c_str())
        .def("get_gen_res_full", &GridModel::get_gen_res_full, DocGridModel::_internal_do_not_use.c_str())
        .def("get_lineor_res_full", &GridModel::get_lineor_res_full, DocGridModel::_internal_do_not_use.c_str())
        .def("get_lineex_res_full", &GridModel::get_lineex_res_full, DocGridModel::_internal_do_not_use.c_str())
        .def("get_trafohv_res_full", &GridModel::get_trafohv_res_full, DocGridModel::_internal_do_not_use.c_str())
        .def("get_trafolv_res_full", &GridModel::get_trafolv_res_full, DocGridModel::_internal_do_not_use.c_str())
        .def("get_storages_res_full", &GridModel::get_storages_res_full, DocGridModel::_internal_do_not_use.c_str())
        .def("get_sgens_res_full", &GridModel::get_sgens_res_full, DocGridModel::_internal_do_not_use.c_str())
        .def("get_dclineor_res_full", &GridModel::get_dclineor_res_full, DocGridModel::_internal_do_not_use.c_str())
        .def("get_dclineex_res_full", &GridModel::get_dclineex_res_full, DocGridModel::_internal_do_not_use.c_str())
        
        // do something with the grid
        // .def("init_Ybus", &DataModel::init_Ybus) // temporary
        .def("deactivate_result_computation", &GridModel::deactivate_result_computation, DocGridModel::deactivate_result_computation.c_str())
        .def("reactivate_result_computation", &GridModel::reactivate_result_computation, DocGridModel::reactivate_result_computation.c_str())
        .def("dc_pf", &GridModel::dc_pf, DocGridModel::dc_pf.c_str())
        .def("ac_pf", &GridModel::ac_pf, DocGridModel::ac_pf.c_str())
        .def("unset_changes", &GridModel::unset_changes, DocGridModel::_internal_do_not_use.c_str())
        .def("tell_recompute_ybus", &GridModel::tell_recompute_ybus, DocGridModel::_internal_do_not_use.c_str())
        .def("tell_recompute_sbus", &GridModel::tell_recompute_sbus, DocGridModel::_internal_do_not_use.c_str())
        .def("tell_solver_need_reset", &GridModel::tell_solver_need_reset, DocGridModel::_internal_do_not_use.c_str())
        .def("tell_ybus_change_sparsity_pattern", &GridModel::tell_ybus_change_sparsity_pattern, DocGridModel::_internal_do_not_use.c_str())
        .def("get_solver_control", &GridModel::get_solver_control, "TODO")
        .def("compute_newton", &GridModel::ac_pf, DocGridModel::ac_pf.c_str())
        .def("get_ptdf", &GridModel::get_ptdf, DocGridModel::_internal_do_not_use.c_str()) // TODO PTDF
        .def("get_Bf", &GridModel::get_Bf, DocGridModel::_internal_do_not_use.c_str()) // TODO PTDF

         // apply action faster (optimized for grid2op representation)
         // it is not recommended to use it outside of grid2Op.
        .def("update_bus_status", &GridModel::update_bus_status, DocGridModel::_internal_do_not_use.c_str())
        .def("update_gens_p", &GridModel::update_gens_p, DocGridModel::_internal_do_not_use.c_str())
        .def("update_sgens_p", &GridModel::update_sgens_p, DocGridModel::_internal_do_not_use.c_str())
        .def("update_gens_v", &GridModel::update_gens_v, DocGridModel::_internal_do_not_use.c_str())
        .def("update_loads_p", &GridModel::update_loads_p, DocGridModel::_internal_do_not_use.c_str())
        .def("update_loads_q", &GridModel::update_loads_q, DocGridModel::_internal_do_not_use.c_str())
        .def("update_topo", &GridModel::update_topo, DocGridModel::_internal_do_not_use.c_str())
        .def("update_storages_p", &GridModel::update_storages_p, DocGridModel::_internal_do_not_use.c_str())

        // auxiliary functions
        .def("set_n_sub", &GridModel::set_n_sub, DocGridModel::_internal_do_not_use.c_str())
        .def("set_max_nb_bus_per_sub", &GridModel::set_max_nb_bus_per_sub, DocGridModel::_internal_do_not_use.c_str())
        .def("set_load_pos_topo_vect", &GridModel::set_load_pos_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_gen_pos_topo_vect", &GridModel::set_gen_pos_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_line_or_pos_topo_vect", &GridModel::set_line_or_pos_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_line_ex_pos_topo_vect", &GridModel::set_line_ex_pos_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_trafo_hv_pos_topo_vect", &GridModel::set_trafo_hv_pos_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_trafo_lv_pos_topo_vect", &GridModel::set_trafo_lv_pos_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_storage_pos_topo_vect", &GridModel::set_storage_pos_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_load_to_subid", &GridModel::set_load_to_subid, DocGridModel::_internal_do_not_use.c_str())
        .def("set_gen_to_subid", &GridModel::set_gen_to_subid, DocGridModel::_internal_do_not_use.c_str())
        .def("set_line_or_to_subid", &GridModel::set_line_or_to_subid, DocGridModel::_internal_do_not_use.c_str())
        .def("set_line_ex_to_subid", &GridModel::set_line_ex_to_subid, DocGridModel::_internal_do_not_use.c_str())
        .def("set_trafo_hv_to_subid", &GridModel::set_trafo_hv_to_subid, DocGridModel::_internal_do_not_use.c_str())
        .def("set_trafo_lv_to_subid", &GridModel::set_trafo_lv_to_subid, DocGridModel::_internal_do_not_use.c_str())
        .def("set_storage_to_subid", &GridModel::set_storage_to_subid, DocGridModel::_internal_do_not_use.c_str())

        // debug function (might disappear without further notice)
        .def("debug_get_Bp_python", &GridModel::debug_get_Bp_python, DocGridModel::_internal_do_not_use.c_str())
        .def("debug_get_Bpp_python", &GridModel::debug_get_Bpp_python, DocGridModel::_internal_do_not_use.c_str())
        ;

    py::class_<TimeSeries>(m, "TimeSeriesCPP", DocComputers::Computers.c_str())
        .def(py::init<const GridModel &>())

        // solver control
        .def("change_solver", &TimeSeries::change_solver, DocGridModel::change_solver.c_str())
        .def("available_solvers", &TimeSeries::available_solvers, DocGridModel::available_solvers.c_str())
        .def("get_solver_type", &TimeSeries::get_solver_type, DocGridModel::get_solver_type.c_str())

        // timers
        .def("total_time", &TimeSeries::total_time, DocComputers::total_time.c_str())
        .def("solver_time", &TimeSeries::solver_time, DocComputers::solver_time.c_str())
        .def("preprocessing_time", &TimeSeries::preprocessing_time, DocComputers::preprocessing_time.c_str())
        .def("amps_computation_time", &TimeSeries::amps_computation_time, DocComputers::amps_computation_time.c_str())
        .def("nb_solved", &TimeSeries::nb_solved, DocComputers::nb_solved.c_str())

        // status
        .def("get_status", &TimeSeries::get_status, DocComputers::get_status.c_str())
        .def("clear", &TimeSeries::clear, DocComputers::clear.c_str())
        .def("close", &TimeSeries::clear, DocComputers::clear.c_str())

        // perform the computations
        .def("compute_Vs", &TimeSeries::compute_Vs, py::call_guard<py::gil_scoped_release>(), DocComputers::compute_Vs.c_str())
        .def("compute_flows", &TimeSeries::compute_flows, py::call_guard<py::gil_scoped_release>(), DocComputers::compute_flows.c_str())
        .def("compute_power_flows", &TimeSeries::compute_power_flows, DocComputers::compute_power_flows.c_str())  // need to be done after "compute_Vs"  and "compute_flows"
        
        // results (for now only flow (at each -line origin- or voltages -at each buses)
        .def("get_flows", &TimeSeries::get_flows, DocComputers::get_flows.c_str())  // need to be done after "compute_Vs"  and "compute_flows"
        .def("get_power_flows", &TimeSeries::get_power_flows, DocComputers::get_power_flows.c_str())  // need to be done after "compute_Vs"  and "compute_flows"
        .def("get_voltages", &TimeSeries::get_voltages, DocComputers::get_voltages.c_str())  // need to be done after "compute_Vs" 
        .def("get_sbuses", &TimeSeries::get_sbuses, DocComputers::get_sbuses.c_str())  // need to be done after "compute_Vs" 
        ;

    py::class_<ContingencyAnalysis>(m, "ContingencyAnalysisCPP", DocSecurityAnalysis::SecurityAnalysis.c_str())
        .def(py::init<const GridModel &>())
        // solver control
        .def("change_solver", &ContingencyAnalysis::change_solver, DocGridModel::change_solver.c_str())
        .def("available_solvers", &ContingencyAnalysis::available_solvers, DocGridModel::available_solvers.c_str())
        .def("get_solver_type", &ContingencyAnalysis::get_solver_type, DocGridModel::get_solver_type.c_str())

        // add some defaults
        .def("add_all_n1", &ContingencyAnalysis::add_all_n1, DocSecurityAnalysis::add_all_n1.c_str())
        .def("add_n1", &ContingencyAnalysis::add_n1, DocSecurityAnalysis::add_n1.c_str())
        .def("add_nk", &ContingencyAnalysis::add_nk, DocSecurityAnalysis::add_nk.c_str())
        .def("add_multiple_n1", &ContingencyAnalysis::add_multiple_n1, DocSecurityAnalysis::add_multiple_n1.c_str())

        // remove some defaults (TODO)
        .def("reset", &ContingencyAnalysis::clear, DocSecurityAnalysis::clear.c_str())
        .def("clear", &ContingencyAnalysis::clear, DocSecurityAnalysis::clear.c_str())
        .def("close", &ContingencyAnalysis::clear, DocComputers::clear.c_str())
        .def("remove_n1", &ContingencyAnalysis::remove_n1, DocSecurityAnalysis::remove_n1.c_str())
        .def("remove_nk", &ContingencyAnalysis::remove_nk, DocSecurityAnalysis::remove_nk.c_str())
        .def("remove_multiple_n1", &ContingencyAnalysis::remove_multiple_n1, DocSecurityAnalysis::remove_multiple_n1.c_str())
        
        // inspect the class
        .def("my_defaults", &ContingencyAnalysis::my_defaults_vect, DocSecurityAnalysis::my_defaults_vect.c_str())

        // perform the computation
        .def("compute", &ContingencyAnalysis::compute, py::call_guard<py::gil_scoped_release>(), DocSecurityAnalysis::compute.c_str())
        .def("compute_flows", &ContingencyAnalysis::compute_flows, py::call_guard<py::gil_scoped_release>(), DocSecurityAnalysis::compute_flows.c_str())
        .def("compute_power_flows", &ContingencyAnalysis::compute_power_flows, DocSecurityAnalysis::compute_power_flows.c_str())

        // results (for now only flow (at each -line origin- or voltages -at each buses)
        .def("get_flows", &ContingencyAnalysis::get_flows, DocSecurityAnalysis::get_flows.c_str())
        .def("get_voltages", &ContingencyAnalysis::get_voltages, DocSecurityAnalysis::get_voltages.c_str())
        .def("get_power_flows", &ContingencyAnalysis::get_power_flows, DocSecurityAnalysis::get_power_flows.c_str())

        // timers
        .def("total_time", &ContingencyAnalysis::total_time, DocComputers::total_time.c_str())
        .def("solver_time", &ContingencyAnalysis::solver_time, DocComputers::solver_time.c_str())
        .def("preprocessing_time", &ContingencyAnalysis::preprocessing_time, DocSecurityAnalysis::preprocessing_time.c_str())
        .def("amps_computation_time", &ContingencyAnalysis::amps_computation_time, DocComputers::amps_computation_time.c_str())
        .def("modif_Ybus_time", &ContingencyAnalysis::modif_Ybus_time, DocSecurityAnalysis::modif_Ybus_time.c_str())
        .def("nb_solved", &ContingencyAnalysis::nb_solved, DocComputers::nb_solved.c_str())
        ;
}
