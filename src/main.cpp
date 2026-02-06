// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "ChooseSolver.hpp"
#include "DataConverter.hpp"
#include "GridModel.hpp"

#include "batch_algorithm/TimeSeries.hpp"
#include "batch_algorithm/ContingencyAnalysis.hpp"

#include "help_fun_msg.hpp"

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
        .def("get_timers_ptdf_lodf", &ChooseSolver::get_timers_ptdf_lodf, "TODO")
        .def("get_fdpf_xb_lu", &ChooseSolver::get_fdpf_xb_lu, py::return_value_policy::reference, DocGridModel::_internal_do_not_use.c_str())  // TODO this for all solver !
        .def("get_fdpf_bx_lu", &ChooseSolver::get_fdpf_bx_lu, py::return_value_policy::reference, DocGridModel::_internal_do_not_use.c_str());

    // iterator for generators
    py::class_<GeneratorContainer>(m, "GeneratorContainer", DocIterator::GeneratorContainer.c_str())
        .def("__len__", [](const GeneratorContainer & data) { return data.nb(); })
        .def("__getitem__", [](const GeneratorContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const GeneratorContainer & data)  {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def("get_bus_id", &GeneratorContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def(py::pickle(
                [](const GeneratorContainer &sub) { // __getstate__
                    // Return a tuple that fully encodes the state of the object
                    return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, sub.get_state());
                },
                [](py::tuple py_state) { // __setstate__
                    if (py_state.size() != 4){
                        std::cout << "GeneratorContainer.__setstate__ : state size " << py_state.size() << std::endl;
                        throw std::runtime_error("Invalid state size when loading GeneratorContainer");
                        }
                    // Create a new C++ instance
                    GeneratorContainer res = GeneratorContainer();
                    // TODO check the size of the input tuple!

                    // now set the status
                    std::string major = py_state[0].cast<std::string>();
                    if (major != VERSION_MAJOR){
                        throw std::runtime_error("Invalid state size when loading GeneratorContainer: wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string minor = py_state[1].cast<std::string>();
                    if (minor != VERSION_MEDIUM){
                        throw std::runtime_error("Invalid state size when loading GeneratorContainer: wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string patch = py_state[2].cast<std::string>();
                    if (patch != VERSION_MINOR){
                        throw std::runtime_error("Invalid state size when loading GeneratorContainer: wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
                    }
                    GeneratorContainer::StateRes state = py_state[3].cast<GeneratorContainer::StateRes>();
                    res.set_state(state);
                    return res;
        })) 
        ; 

    py::class_<GenInfo>(m, "GenInfo", DocIterator::GenInfo.c_str())
        .def_readonly("id", &GenInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &GenInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &GenInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &GenInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &GenInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &GenInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("is_slack", &GenInfo::is_slack, DocIterator::is_slack.c_str())
        .def_readonly("slack_weight", &GenInfo::slack_weight, DocIterator::slack_weight.c_str())
        .def_readonly("voltage_regulator_on", &GenInfo::voltage_regulator_on, "TODO")
        .def_readonly("target_p_mw", &GenInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_vm_pu", &GenInfo::target_vm_pu, DocIterator::target_vm_pu.c_str())
        .def_readonly("target_q_mvar", &GenInfo::target_q_mvar, "TODO")
        .def_readonly("min_q_mvar", &GenInfo::min_q_mvar, DocIterator::min_q_mvar.c_str())
        .def_readonly("max_q_mvar", &GenInfo::max_q_mvar, DocIterator::max_q_mvar.c_str())
        .def_readonly("has_res", &GenInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &GenInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &GenInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &GenInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &GenInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        // pypowsybl
        .def_readonly("voltage_level_id", &GenInfo::sub_id, DocIterator::sub_id.c_str());

    // iterator for sgens
    py::class_<SGenContainer>(m, "SGenContainer", DocIterator::SGenContainer.c_str())
        .def("__len__", [](const SGenContainer & data) { return data.nb(); })
        .def("__getitem__", [](const SGenContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const SGenContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def("get_bus_id", &SGenContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def(py::pickle(
                [](const SGenContainer &sub) { // __getstate__
                    // Return a tuple that fully encodes the state of the object
                    return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, sub.get_state());
                },
                [](py::tuple py_state) { // __setstate__
                    if (py_state.size() != 4){
                        std::cout << "SGenContainer.__setstate__ : state size " << py_state.size() << std::endl;
                        throw std::runtime_error("Invalid state size when loading SGenContainer");
                        }
                    // Create a new C++ instance
                    SGenContainer res = SGenContainer();
                    // TODO check the size of the input tuple!

                    // now set the status
                    std::string major = py_state[0].cast<std::string>();
                    if (major != VERSION_MAJOR){
                        throw std::runtime_error("Invalid state size when loading SGenContainer: wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string minor = py_state[1].cast<std::string>();
                    if (minor != VERSION_MEDIUM){
                        throw std::runtime_error("Invalid state size when loading SGenContainer: wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string patch = py_state[2].cast<std::string>();
                    if (patch != VERSION_MINOR){
                        throw std::runtime_error("Invalid state size when loading SGenContainer: wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
                    }
                    SGenContainer::StateRes state = py_state[3].cast<SGenContainer::StateRes>();
                    res.set_state(state);
                    return res;
        })) 
        ; 

    py::class_<SGenInfo>(m, "SGenInfo", DocIterator::SGenInfo.c_str())
        .def_readonly("id", &SGenInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &SGenInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &SGenInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &SGenInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &SGenInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &SGenInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("min_q_mvar", &SGenInfo::min_q_mvar, DocIterator::min_q_mvar.c_str())
        .def_readonly("max_q_mvar", &SGenInfo::max_q_mvar, DocIterator::max_q_mvar.c_str())
        .def_readonly("min_p_mw", &SGenInfo::min_p_mw, DocIterator::min_p_mw.c_str())
        .def_readonly("max_p_mw", &SGenInfo::max_p_mw, DocIterator::max_p_mw.c_str())
        .def_readonly("target_p_mw", &SGenInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &SGenInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &SGenInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &SGenInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &SGenInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &SGenInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &SGenInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        // pypowsybl
        .def_readonly("voltage_level_id", &SGenInfo::sub_id, DocIterator::sub_id.c_str())
        ;

    // iterator for loads (and storage units)
    py::class_<LoadContainer>(m, "LoadContainer", DocIterator::LoadContainer.c_str())
        .def("__len__", [](const LoadContainer & data) { return data.nb(); })
        .def("__getitem__", [](const LoadContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const LoadContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def("get_bus_id", &LoadContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def(py::pickle(
                [](const LoadContainer &sub) { // __getstate__
                    // Return a tuple that fully encodes the state of the object
                    return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, sub.get_state());
                },
                [](py::tuple py_state) { // __setstate__
                    if (py_state.size() != 4){
                        std::cout << "LoadContainer.__setstate__ : state size " << py_state.size() << std::endl;
                        throw std::runtime_error("Invalid state size when loading LoadContainer");
                        }
                    // Create a new C++ instance
                    LoadContainer res = LoadContainer();
                    // TODO check the size of the input tuple!

                    // now set the status
                    std::string major = py_state[0].cast<std::string>();
                    if (major != VERSION_MAJOR){
                        throw std::runtime_error("Invalid state size when loading LoadContainer: wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string minor = py_state[1].cast<std::string>();
                    if (minor != VERSION_MEDIUM){
                        throw std::runtime_error("Invalid state size when loading LoadContainer: wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string patch = py_state[2].cast<std::string>();
                    if (patch != VERSION_MINOR){
                        throw std::runtime_error("Invalid state size when loading LoadContainer: wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
                    }
                    LoadContainer::StateRes state = py_state[3].cast<LoadContainer::StateRes>();
                    res.set_state(state);
                    return res;
        })) 
        ; 

    py::class_<LoadInfo>(m, "LoadInfo", DocIterator::LoadInfo.c_str())
        .def_readonly("id", &LoadInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &LoadInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &LoadInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &LoadInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &LoadInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &LoadInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("target_p_mw", &LoadInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &LoadInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &LoadInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &LoadInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &LoadInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &LoadInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &LoadInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        // pypowsybl
        .def_readonly("voltage_level_id", &LoadInfo::sub_id, DocIterator::sub_id.c_str())
        ;

    // iterator for shunts
    py::class_<ShuntContainer>(m, "ShuntContainer", DocIterator::ShuntContainer.c_str())
        .def("__len__", [](const ShuntContainer & data) { return data.nb(); })
        .def("__getitem__", [](const ShuntContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const ShuntContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def("get_bus_id", &ShuntContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def(py::pickle(
                [](const ShuntContainer &sub) { // __getstate__
                    // Return a tuple that fully encodes the state of the object
                    return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, sub.get_state());
                },
                [](py::tuple py_state) { // __setstate__
                    if (py_state.size() != 4){
                        std::cout << "ShuntContainer.__setstate__ : state size " << py_state.size() << std::endl;
                        throw std::runtime_error("Invalid state size when loading ShuntContainer");
                        }
                    // Create a new C++ instance
                    ShuntContainer res = ShuntContainer();
                    // TODO check the size of the input tuple!

                    // now set the status
                    std::string major = py_state[0].cast<std::string>();
                    if (major != VERSION_MAJOR){
                        throw std::runtime_error("Invalid state size when loading ShuntContainer: wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string minor = py_state[1].cast<std::string>();
                    if (minor != VERSION_MEDIUM){
                        throw std::runtime_error("Invalid state size when loading ShuntContainer: wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string patch = py_state[2].cast<std::string>();
                    if (patch != VERSION_MINOR){
                        throw std::runtime_error("Invalid state size when loading ShuntContainer: wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
                    }
                    ShuntContainer::StateRes state = py_state[3].cast<ShuntContainer::StateRes>();
                    res.set_state(state);
                    return res;
        }))      
        ; 

    py::class_<ShuntInfo>(m, "ShuntInfo", DocIterator::ShuntInfo.c_str())
        .def_readonly("id", &ShuntInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &ShuntInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &ShuntInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &ShuntInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &ShuntInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &ShuntInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("target_p_mw", &ShuntInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &ShuntInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &ShuntInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &ShuntInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &ShuntInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &ShuntInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &ShuntInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        // pypowsybl
        .def_readonly("voltage_level_id", &ShuntInfo::sub_id, DocIterator::sub_id.c_str())
        ;

    // iterator for trafos
    py::class_<TrafoContainer>(m, "TrafoContainer", DocIterator::TrafoContainer.c_str())
        .def("__len__", [](const TrafoContainer & data) { return data.nb(); })
        .def("__getitem__", [](const TrafoContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const TrafoContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def_property_readonly("ignore_tap_side_for_shift", &TrafoContainer::ignore_tap_side_for_shift, 
            R"mydelimiter(
            Whether ignore the tap side is ignored when using the 
            'shift' attribute (should be True for pandapower, 
            where it is ignored and False otherwise).)mydelimiter")
        .def("get_bus_id_side_1", &TrafoContainer::get_bus_id_side_1_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def("get_bus_id_side_2", &TrafoContainer::get_bus_id_side_2_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def(py::pickle(
                [](const TrafoContainer &sub) { // __getstate__
                    // Return a tuple that fully encodes the state of the object
                    return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, sub.get_state());
                },
                [](py::tuple py_state) { // __setstate__
                    if (py_state.size() != 4){
                        std::cout << "TrafoContainer.__setstate__ : state size " << py_state.size() << std::endl;
                        throw std::runtime_error("Invalid state size when loading TrafoContainer");
                        }
                    // Create a new C++ instance
                    TrafoContainer res = TrafoContainer();
                    // TODO check the size of the input tuple!

                    // now set the status
                    std::string major = py_state[0].cast<std::string>();
                    if (major != VERSION_MAJOR){
                        throw std::runtime_error("Invalid state size when loading TrafoContainer: wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string minor = py_state[1].cast<std::string>();
                    if (minor != VERSION_MEDIUM){
                        throw std::runtime_error("Invalid state size when loading TrafoContainer: wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string patch = py_state[2].cast<std::string>();
                    if (patch != VERSION_MINOR){
                        throw std::runtime_error("Invalid state size when loading TrafoContainer: wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
                    }
                    TrafoContainer::StateRes state = py_state[3].cast<TrafoContainer::StateRes>();
                    res.set_state(state);
                    return res;
        }))       
        ; 

    py::class_<TrafoInfo>(m, "TrafoInfo", DocIterator::TrafoInfo.c_str())
        .def_readonly("id", &TrafoInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &TrafoInfo::name, DocIterator::name.c_str())
        .def_readonly("sub1_id", &TrafoInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("sub2_id", &TrafoInfo::sub_2_id, DocIterator::sub_id.c_str())
        .def_readonly("pos1_topo_vect", &TrafoInfo::pos_1_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("pos2_topo_vect", &TrafoInfo::pos_2_topo_vect, DocIterator::pos_topo_vect.c_str())

        .def_readonly("connected_global", &TrafoInfo::connected_global, DocIterator::connected.c_str())
        .def_readonly("connected1", &TrafoInfo::connected_1, DocIterator::connected.c_str())
        .def_readonly("connected2", &TrafoInfo::connected_2, DocIterator::connected.c_str())

        .def_readonly("bus1_id", &TrafoInfo::bus_1_id, DocIterator::bus_hv_id.c_str())
        .def_readonly("bus2_id", &TrafoInfo::bus_2_id, DocIterator::bus_lv_id.c_str())

        .def_readonly("r_pu", &TrafoInfo::r_pu, DocIterator::r_pu.c_str())
        .def_readonly("x_pu", &TrafoInfo::x_pu, DocIterator::x_pu.c_str())
        .def_readonly("h1_pu", &TrafoInfo::h1_pu, DocIterator::h_pu.c_str())
        .def_readonly("h2_pu", &TrafoInfo::h2_pu, DocIterator::h_pu.c_str())
        .def_readonly("is_tap_side_1", &TrafoInfo::is_tap_side1, DocIterator::is_tap_hv_side.c_str())
        .def_readonly("ratio", &TrafoInfo::ratio, DocIterator::ratio.c_str())
        .def_readonly("shift_rad", &TrafoInfo::shift_rad, DocIterator::shift_rad.c_str())
        .def_readonly("has_res", &TrafoInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p1_mw", &TrafoInfo::res_p1_mw, DocIterator::res_p_hv_mw.c_str())
        .def_readonly("res_q1_mvar", &TrafoInfo::res_q1_mvar, DocIterator::res_q_hv_mvar.c_str())
        .def_readonly("res_v1_kv", &TrafoInfo::res_v1_kv, DocIterator::res_v_hv_kv.c_str())
        .def_readonly("res_a1_ka", &TrafoInfo::res_a1_ka, DocIterator::res_a_hv_ka.c_str())
        .def_readonly("res_p2_mw", &TrafoInfo::res_p2_mw, DocIterator::res_p_lv_mw.c_str())
        .def_readonly("res_q2_mvar", &TrafoInfo::res_q2_mvar, DocIterator::res_q_lv_mvar.c_str())
        .def_readonly("res_v2_kv", &TrafoInfo::res_v2_kv, DocIterator::res_v_lv_kv.c_str())
        .def_readonly("res_a2_ka", &TrafoInfo::res_a2_ka, DocIterator::res_a_lv_ka.c_str())
        .def_readonly("res_theta1_deg", &TrafoInfo::res_theta1_deg, DocIterator::res_theta_hv_deg.c_str())
        .def_readonly("res_theta2_deg", &TrafoInfo::res_theta2_deg, DocIterator::res_theta_lv_deg.c_str())

        .def_readonly("yac_11", &TrafoInfo::yac_11, "TODO doc")
        .def_readonly("yac_12", &TrafoInfo::yac_12, "TODO doc")
        .def_readonly("yac_21", &TrafoInfo::yac_21, "TODO doc")
        .def_readonly("yac_22", &TrafoInfo::yac_22, "TODO doc")
        .def_readonly("ydc_11", &TrafoInfo::ydc_11, "TODO doc")
        .def_readonly("ydc_12", &TrafoInfo::ydc_12, "TODO doc")
        .def_readonly("ydc_21", &TrafoInfo::ydc_21, "TODO doc")
        .def_readonly("ydc_22", &TrafoInfo::ydc_22, "TODO doc")

        // pypowsybl
        .def_readonly("voltage_level1_id", &TrafoInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("voltage_level2_id", &TrafoInfo::sub_2_id, DocIterator::sub_id.c_str())
        ;

    // iterator for trafos
    py::class_<LineContainer>(m, "LineContainer", DocIterator::LineContainer.c_str())
        .def("__len__", [](const LineContainer & data) { return data.nb(); })
        .def("__getitem__", [](const LineContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const LineContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def("get_bus_id_side_1", &LineContainer::get_bus_id_side_1_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def("get_bus_id_side_2", &LineContainer::get_bus_id_side_2_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def(py::pickle(
                [](const LineContainer &sub) { // __getstate__
                    // Return a tuple that fully encodes the state of the object
                    return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, sub.get_state());
                },
                [](py::tuple py_state) { // __setstate__
                    if (py_state.size() != 4){
                        std::cout << "LineContainer.__setstate__ : state size " << py_state.size() << std::endl;
                        throw std::runtime_error("Invalid state size when loading LineContainer");
                        }
                    // Create a new C++ instance
                    LineContainer res = LineContainer();
                    // TODO check the size of the input tuple!

                    // now set the status
                    std::string major = py_state[0].cast<std::string>();
                    if (major != VERSION_MAJOR){
                        throw std::runtime_error("Invalid state size when loading LineContainer: wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string minor = py_state[1].cast<std::string>();
                    if (minor != VERSION_MEDIUM){
                        throw std::runtime_error("Invalid state size when loading LineContainer: wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string patch = py_state[2].cast<std::string>();
                    if (patch != VERSION_MINOR){
                        throw std::runtime_error("Invalid state size when loading LineContainer: wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
                    }
                    LineContainer::StateRes state = py_state[3].cast<LineContainer::StateRes>();
                    res.set_state(state);
                    return res;
        }))  
        ;

    py::class_<LineInfo>(m, "LineInfo", DocIterator::LineInfo.c_str())
        .def_readonly("id", &LineInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &LineInfo::name, DocIterator::name.c_str())
        .def_readonly("sub1_id", &LineInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("sub2_id", &LineInfo::sub_2_id, DocIterator::sub_id.c_str())
        .def_readonly("pos1_topo_vect", &LineInfo::pos_1_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("pos2_topo_vect", &LineInfo::pos_2_topo_vect, DocIterator::pos_topo_vect.c_str())

        .def_readonly("connected_global", &LineInfo::connected_global, DocIterator::connected.c_str())
        .def_readonly("connected1", &LineInfo::connected_1, DocIterator::connected.c_str())
        .def_readonly("connected2", &LineInfo::connected_2, DocIterator::connected.c_str())

        .def_readonly("bus1_id", &LineInfo::bus_1_id, DocIterator::bus_or_id.c_str())
        .def_readonly("bus2_id", &LineInfo::bus_2_id, DocIterator::bus_ex_id.c_str())

        .def_readonly("r_pu", &LineInfo::r_pu, DocIterator::r_pu.c_str())
        .def_readonly("x_pu", &LineInfo::x_pu, DocIterator::x_pu.c_str())
        .def_readonly("h1_pu", &LineInfo::h1_pu, DocIterator::h_pu.c_str())
        .def_readonly("h2_pu", &LineInfo::h2_pu, DocIterator::h_pu.c_str())
        .def_readonly("has_res", &LineInfo::has_res, DocIterator::has_res.c_str())

        .def_readonly("res_p1_mw", &LineInfo::res_p1_mw, DocIterator::res_p_or_mw.c_str())
        .def_readonly("res_q1_mvar", &LineInfo::res_q1_mvar, DocIterator::res_q_or_mvar.c_str())
        .def_readonly("res_v1_kv", &LineInfo::res_v1_kv, DocIterator::res_v_or_kv.c_str())
        .def_readonly("res_a1_ka", &LineInfo::res_a1_ka, DocIterator::res_a_or_ka.c_str())
        .def_readonly("res_p2_mw", &LineInfo::res_p2_mw, DocIterator::res_p_ex_mw.c_str())
        .def_readonly("res_q2_mvar", &LineInfo::res_q2_mvar, DocIterator::res_q_ex_mvar.c_str())
        .def_readonly("res_v2_kv", &LineInfo::res_v2_kv, DocIterator::res_v_ex_kv.c_str())
        .def_readonly("res_a2_ka", &LineInfo::res_a2_ka, DocIterator::res_a_ex_ka.c_str())
        .def_readonly("res_theta1_deg", &LineInfo::res_theta1_deg, DocIterator::res_theta_or_deg.c_str())
        .def_readonly("res_theta2_deg", &LineInfo::res_theta2_deg, DocIterator::res_theta_ex_deg.c_str())

        .def_readonly("yac_11", &LineInfo::yac_11, "TODO doc")
        .def_readonly("yac_12", &LineInfo::yac_12, "TODO doc")
        .def_readonly("yac_21", &LineInfo::yac_21, "TODO doc")
        .def_readonly("yac_22", &LineInfo::yac_22, "TODO doc")
        .def_readonly("ydc_11", &LineInfo::ydc_11, "TODO doc")
        .def_readonly("ydc_12", &LineInfo::ydc_12, "TODO doc")
        .def_readonly("ydc_21", &LineInfo::ydc_21, "TODO doc")
        .def_readonly("ydc_22", &LineInfo::ydc_22, "TODO doc")

        // pypowsybl
        .def_readonly("voltage_level1_id", &LineInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("voltage_level2_id", &LineInfo::sub_2_id, DocIterator::sub_id.c_str())
        ;

    // iterator for dc lines
    py::class_<DCLineContainer>(m, "DCLineContainer", DocIterator::DCLineContainer.c_str())
        .def("__len__", [](const DCLineContainer & data) { return data.nb(); })
        .def("__getitem__", [](const DCLineContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const DCLineContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        .def("get_bus_id_side_1", &DCLineContainer::get_bus_id_side_1_numpy)
        .def("get_bus_id_side_2", &DCLineContainer::get_bus_id_side_2_numpy)
        .def(py::pickle(
                [](const DCLineContainer &sub) { // __getstate__
                    // Return a tuple that fully encodes the state of the object
                    return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, sub.get_state());
                },
                [](py::tuple py_state) { // __setstate__
                    if (py_state.size() != 4){
                        std::cout << "DCLineContainer.__setstate__ : state size " << py_state.size() << std::endl;
                        throw std::runtime_error("Invalid state size when loading DCLineContainer");
                        }
                    // Create a new C++ instance
                    DCLineContainer res = DCLineContainer();
                    // TODO check the size of the input tuple!

                    // now set the status
                    std::string major = py_state[0].cast<std::string>();
                    if (major != VERSION_MAJOR){
                        throw std::runtime_error("Invalid state size when loading DCLineContainer: wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string minor = py_state[1].cast<std::string>();
                    if (minor != VERSION_MEDIUM){
                        throw std::runtime_error("Invalid state size when loading DCLineContainer: wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string patch = py_state[2].cast<std::string>();
                    if (patch != VERSION_MINOR){
                        throw std::runtime_error("Invalid state size when loading DCLineContainer: wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
                    }
                    DCLineContainer::StateRes state = py_state[3].cast<DCLineContainer::StateRes>();
                    res.set_state(state);
                    return res;
        }))  
        ;

    py::class_<DCLineInfo>(m, "DCLineInfo", DocIterator::DCLineInfo.c_str())
        .def_readonly("id", &DCLineInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &DCLineInfo::name, DocIterator::name.c_str())
        .def_readonly("sub1_id", &DCLineInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("sub2_id", &DCLineInfo::sub_2_id, DocIterator::sub_id.c_str())
        .def_readonly("pos1_topo_vect", &DCLineInfo::pos_1_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("pos2_topo_vect", &DCLineInfo::pos_2_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected_global", &DCLineInfo::connected_global, DocIterator::connected.c_str())
        .def_readonly("connected1", &DCLineInfo::connected_1, DocIterator::connected.c_str())
        .def_readonly("connected2", &DCLineInfo::connected_2, DocIterator::connected.c_str())
        .def_readonly("bus1_id", &DCLineInfo::bus_1_id, DocIterator::bus_or_id.c_str())
        .def_readonly("bus2_id", &DCLineInfo::bus_2_id, DocIterator::bus_ex_id.c_str())
        .def_readonly("target_p1_mw", &DCLineInfo::target_p_1_mw, DocIterator::target_p_or_mw.c_str())
        .def_readonly("p2_mw", &DCLineInfo::p_2_mw, DocIterator::target_p_or_mw.c_str())
        .def_readonly("target_vm1_pu", &DCLineInfo::target_vm_1_pu, DocIterator::target_vm_or_pu.c_str())
        .def_readonly("target_vm2_pu", &DCLineInfo::target_vm_2_pu, DocIterator::target_vm_ex_pu.c_str())
        .def_readonly("loss_pct", &DCLineInfo::loss_pct, DocIterator::loss_pct.c_str())
        .def_readonly("loss_mw", &DCLineInfo::loss_mw, DocIterator::loss_mw.c_str())
        .def_readonly("gen1", &DCLineInfo::gen_side_1, DocIterator::gen_or.c_str())
        .def_readonly("gen2", &DCLineInfo::gen_side_2, DocIterator::gen_ex.c_str())
        .def_readonly("has_res", &DCLineInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p1_mw", &DCLineInfo::res_p1_mw, DocIterator::res_p_or_mw_dcline.c_str())
        .def_readonly("res_p2_mw", &DCLineInfo::res_p2_mw, DocIterator::res_p_ex_mw_dcline.c_str())
        .def_readonly("res_q1_mvar", &DCLineInfo::res_q1_mvar, DocIterator::res_q_or_mvar_dcline.c_str())
        .def_readonly("res_q2_mvar", &DCLineInfo::res_q2_mvar, DocIterator::res_q_ex_mvar_dcline.c_str())
        .def_readonly("res_v1_kv", &DCLineInfo::res_v1_kv, DocIterator::res_v_or_kv_dcline.c_str())
        .def_readonly("res_v2_kv", &DCLineInfo::res_v2_kv, DocIterator::res_v_ex_kv_dcline.c_str())
        .def_readonly("res_theta1_deg", &DCLineInfo::res_theta1_deg, DocIterator::res_theta_or_deg_dcline.c_str())
        .def_readonly("res_theta2_deg", &DCLineInfo::res_theta2_deg, DocIterator::res_theta_ex_deg_dcline.c_str())

        // pypowsybl
        .def_readonly("voltage_level1_id", &DCLineInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("voltage_level2_id", &DCLineInfo::sub_2_id, DocIterator::sub_id.c_str())
        ;

    py::class_<SubstationContainer>(m, "SubstationContainer", "TODO")
        .def("__len__", [](const SubstationContainer & data) { return data.nb(); })
        .def("__getitem__", [](const SubstationContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const SubstationContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */        // pickle
        .def(py::pickle(
                [](const SubstationContainer &sub) { // __getstate__
                    // Return a tuple that fully encodes the state of the object
                    return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, sub.get_state());
                },
                [](py::tuple py_state) { // __setstate__
                    if (py_state.size() != 4){
                        std::cout << "SubstationContainer.__setstate__ : state size " << py_state.size() << std::endl;
                        throw std::runtime_error("Invalid state size when loading SubstationContainer");
                        }
                    // Create a new C++ instance
                    SubstationContainer res = SubstationContainer();
                    // TODO check the size of the input tuple!

                    // now set the status
                    std::string major = py_state[0].cast<std::string>();
                    if (major != VERSION_MAJOR){
                        throw std::runtime_error("Invalid state size when loading SubstationContainer: wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string minor = py_state[1].cast<std::string>();
                    if (minor != VERSION_MEDIUM){
                        throw std::runtime_error("Invalid state size when loading SubstationContainer: wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string patch = py_state[2].cast<std::string>();
                    if (patch != VERSION_MINOR){
                        throw std::runtime_error("Invalid state size when loading SubstationContainer: wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
                    }
                    SubstationContainer::StateRes state = py_state[3].cast<SubstationContainer::StateRes>();
                    res.set_state(state);
                    return res;
        })) 
        ;

    py::class_<SubstationInfo>(m, "SubstationInfo", "TODO")
        .def_readonly("id", &SubstationInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &SubstationInfo::name, DocIterator::name.c_str())
        .def_readonly("nb_max_busbars", &SubstationInfo::nb_max_busbars, DocIterator::name.c_str())
        .def_readonly("vn_kv", &SubstationInfo::vn_kv, DocIterator::name.c_str())
        ;

    // converters
    py::class_<PandaPowerConverter>(m, "PandaPowerConverter")
        .def(py::init<>())
        .def("set_f_hz", &PandaPowerConverter::set_f_hz)
        .def("set_sn_mva", &PandaPowerConverter::set_sn_mva)
        .def("get_line_param_legacy", &PandaPowerConverter::get_line_param_legacy)
        .def("get_line_param", &PandaPowerConverter::get_line_param)
        .def("get_trafo_param_pp3", &PandaPowerConverter::get_trafo_param_pp3)
        .def("get_trafo_param_pp2", &PandaPowerConverter::get_trafo_param_pp2);

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
        .def("has_one_el_changed_bus", &SolverControl::has_one_el_changed_bus, "TODO")
        ;

    py::class_<GridModel>(m, "GridModel", DocGridModel::GridModel.c_str())
        .def(py::init<>())
        .def("copy", &GridModel::copy, "TODO", py::return_value_policy::take_ownership)
        .def_property("_ls_to_orig",
                      &GridModel::get_ls_to_orig,
                      &GridModel::set_ls_to_orig,
                      R"mydelimiter(
_ls_to_orig: has the size of the number of possible buses in lightsim2grid 
(*ie* `n_sub_ * max_nb_bus_per_sub_` ) and gives the id of the corresponding
bus in the original grid (pandapower or pypowsybl).

If a "-1" is present, then this bus does not exist in the original grid, 
it is only present in the lightsim2grid gridmodel.
)mydelimiter")
        .def_property("_orig_to_ls",
                      &GridModel::get_orig_to_ls,
                      &GridModel::set_orig_to_ls,
                      R"mydelimiter(
Opposite to _ls_to_orig. The vector _orig_to_ls has the size of the number
of buses in the original grid (pandapower or pypowsybl) and tells 
to which bus of lightsim2grid it corresponds. It should be a >= integer
between 0 and `n_sub_ * max_nb_bus_per_sub_`

)mydelimiter"
                    )
        .def_property("_max_nb_bus_per_sub",
                      &GridModel::get_max_nb_bus_per_sub,
                      &GridModel::set_max_nb_bus_per_sub,
                      "do not modify it after loading !")
        .def_property_readonly("timer_last_ac_pf", &GridModel::timer_last_ac_pf, "TODO")
        .def_property_readonly("timer_last_dc_pf", &GridModel::timer_last_dc_pf, "TODO")

        // pickle
        .def(py::pickle(
                [](const GridModel &sub) { // __getstate__
                    // Return a tuple that fully encodes the state of the object
                    return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, sub.get_state());
                },
                [](py::tuple py_state) { // __setstate__
                    if (py_state.size() != 4){
                        std::cout << "GridModel.__setstate__ : state size " << py_state.size() << std::endl;
                        throw std::runtime_error("Invalid state size when loading GridModel");
                        }
                    // Create a new C++ instance
                    GridModel res = GridModel();
                    // TODO check the size of the input tuple!

                    // now set the status
                    std::string major = py_state[0].cast<std::string>();
                    if (major != VERSION_MAJOR){
                        throw std::runtime_error("Invalid state size when loading GridModel: wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string minor = py_state[1].cast<std::string>();
                    if (minor != VERSION_MEDIUM){
                        throw std::runtime_error("Invalid state size when loading GridModel: wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
                    }
                    std::string patch = py_state[2].cast<std::string>();
                    if (patch != VERSION_MINOR){
                        throw std::runtime_error("Invalid state size when loading GridModel: wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
                    }
                    GridModel::StateRes state = py_state[3].cast<GridModel::StateRes>();
                    res.set_state(state);
                    return res;
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
        .def("init_trafo_pandapower", &GridModel::init_trafo_pandapower, DocGridModel::_internal_do_not_use.c_str())  // same 
        .def("init_trafo", &GridModel::init_trafo, DocGridModel::_internal_do_not_use.c_str())  // same 
        .def("init_generators", &GridModel::init_generators, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_generators_full", &GridModel::init_generators_full, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_loads", &GridModel::init_loads, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_storages", &GridModel::init_storages, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_sgens", &GridModel::init_sgens, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_dclines", &GridModel::init_dclines, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("add_gen_slackbus", &GridModel::add_gen_slackbus, DocGridModel::_internal_do_not_use.c_str()) // same
        .def("remove_gen_slackbus", &GridModel::remove_gen_slackbus, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("get_bus_vn_kv", &GridModel::get_bus_vn_kv, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_bus_status", &GridModel::get_bus_status, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        // inspect the grid
        .def("get_substations", &GridModel::get_substations, "TODO", py::return_value_policy::reference)
        .def("get_lines", &GridModel::get_lines, DocGridModel::get_lines.c_str(), py::return_value_policy::reference)
        .def("get_dclines", &GridModel::get_dclines, DocGridModel::get_dclines.c_str(), py::return_value_policy::reference)
        .def("get_trafos", &GridModel::get_trafos, DocGridModel::get_trafos.c_str(), py::return_value_policy::reference)
        .def("get_generators", &GridModel::get_generators, DocGridModel::get_generators.c_str(), py::return_value_policy::reference)
        .def("get_static_generators", &GridModel::get_static_generators, DocGridModel::get_static_generators.c_str(), py::return_value_policy::reference)
        .def("get_shunts", &GridModel::get_shunts, DocGridModel::get_shunts.c_str(), py::return_value_policy::reference)
        .def("get_storages", &GridModel::get_storages, DocGridModel::get_storages.c_str(), py::return_value_policy::reference)
        .def("get_loads", &GridModel::get_loads, DocGridModel::get_loads.c_str(), py::return_value_policy::reference)

        // pypowsybl compat names
        .def("get_voltage_levels", &GridModel::get_substations, "TODO", py::return_value_policy::reference)
        .def("get_2_windings_transformers", &GridModel::get_trafos, DocGridModel::get_trafos.c_str(), py::return_value_policy::reference)
        .def("get_shunt_compensators", &GridModel::get_shunts, DocGridModel::get_shunts.c_str(), py::return_value_policy::reference)

        // modify the grid
        .def("turnedoff_no_pv", &GridModel::turnedoff_no_pv, "Turned off (or generators with p = 0) generators will not be pv buses, they will not maintain voltage")
        .def("turnedoff_pv", &GridModel::turnedoff_pv, "Turned off (or generators with p = 0) generators will be pv buses, they will maintain voltage (default)")
        .def("get_turnedoff_gen_pv", &GridModel::get_turnedoff_gen_pv, "TODO")
        .def("update_slack_weights", &GridModel::update_slack_weights, "TODO")
        .def("update_slack_weights_by_id", &GridModel::update_slack_weights_by_id, "TODO")
        .def("assign_slack_to_most_connected", &GridModel::assign_slack_to_most_connected, "TODO")
        .def("consider_only_main_component", &GridModel::consider_only_main_component, "TODO and TODO DC LINE: one side might be in the connected comp and not the other !")
        .def("set_ignore_status_global", &GridModel::set_ignore_status_global, "Ignore the 'global_status' flags for powerlines and trafo (set to true if you want to control independantly each side of powerlines and trafo). Default: false.")
        .def("set_synch_status_both_side", &GridModel::set_synch_status_both_side, "Synch the status of each side of the powerlines and trafo. It means that if you disconnect one side of a powerline / trafo, the other side will also be disconnected. Default: true.")
        .def("get_ignore_status_global", &GridModel::get_ignore_status_global, "TODO doc")
        .def("get_synch_status_both_side", &GridModel::get_synch_status_both_side, "TODO doc")

        // names
        .def("set_line_names", &GridModel::set_line_names, "TODO")
        .def("set_dcline_names", &GridModel::set_dcline_names, "TODO")
        .def("set_trafo_names", &GridModel::set_trafo_names, "TODO")
        .def("set_gen_names", &GridModel::set_gen_names, "TODO")
        .def("set_load_names", &GridModel::set_load_names, "TODO")
        .def("set_storage_names", &GridModel::set_storage_names,  "TODO")
        .def("set_sgen_names", &GridModel::set_sgen_names, "TODO")
        .def("set_shunt_names", &GridModel::set_shunt_names, "TODO")
        .def("set_substation_names", &GridModel::set_substation_names, DocGridModel::_internal_do_not_use.c_str())
        .def("get_substation_names", &GridModel::get_substation_names, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_bus", &GridModel::deactivate_bus_python, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_bus", &GridModel::reactivate_bus_python, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_powerline", &GridModel::deactivate_powerline, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_powerline", &GridModel::reactivate_powerline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus1_powerline", &GridModel::change_bus1_powerline_python, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus2_powerline", &GridModel::change_bus2_powerline_python, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus1_powerline", &GridModel::get_bus1_powerline, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_bus2_powerline", &GridModel::get_bus2_powerline, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        .def("deactivate_trafo", &GridModel::deactivate_trafo, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_trafo", &GridModel::reactivate_trafo, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus1_trafo", &GridModel::change_bus1_trafo_python, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus2_trafo", &GridModel::change_bus2_trafo_python, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus1_trafo", &GridModel::get_bus1_trafo, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_bus2_trafo", &GridModel::get_bus2_trafo, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_ratio_trafo", &GridModel::change_ratio_trafo, "TODO")
        .def("change_shift_trafo", &GridModel::change_shift_trafo, 
            R"mydelimiter(
            TODO Change the phase shift ratio for a given transformer. 

            .. warning::
                It should be expressed in rad (not in deg). 
            
            If the flag
            `ignore_tap_side_for_shift` (*eg* gridmodel.get_trafos().ignore_tap_side_for_shift) 
            is set to False (should be default), then the ratio should be given 
            at the side of the tap (side1 or side2). If this
            flag is True (*eg* the grid comes from pandapower) then the phase
            shift ratio should be given in in the side1 (hv side in pandapower).
            )mydelimiter")
        .def("change_shift_trafo_deg", &GridModel::change_shift_trafo_deg, 
            "Same as :ref:`change_shift_trafo` but phase shift is expressed in degree and NOT in rad.")
        .def("deactivate_load", &GridModel::deactivate_load, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_load", &GridModel::reactivate_load, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_load", &GridModel::change_bus_load_python, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_load", &GridModel::get_bus_load, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_p_load", &GridModel::change_p_load, DocGridModel::_internal_do_not_use.c_str())
        .def("change_q_load", &GridModel::change_q_load, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_gen", &GridModel::deactivate_gen, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_gen", &GridModel::reactivate_gen, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_gen", &GridModel::change_bus_gen_python, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_gen", &GridModel::get_bus_gen, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_p_gen", &GridModel::change_p_gen, DocGridModel::_internal_do_not_use.c_str())
        .def("change_v_gen", &GridModel::change_v_gen, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_shunt", &GridModel::deactivate_shunt, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_shunt", &GridModel::reactivate_shunt, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_shunt", &GridModel::change_bus_shunt_python, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_shunt", &GridModel::get_bus_shunt, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_p_shunt", &GridModel::change_p_shunt, DocGridModel::_internal_do_not_use.c_str())
        .def("change_q_shunt", &GridModel::change_q_shunt, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_sgen", &GridModel::deactivate_sgen, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_sgen", &GridModel::reactivate_sgen, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_sgen", &GridModel::change_bus_sgen_python, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_sgen", &GridModel::get_bus_sgen, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_p_sgen", &GridModel::change_p_sgen, DocGridModel::_internal_do_not_use.c_str())
        .def("change_q_sgen", &GridModel::change_q_sgen, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_storage", &GridModel::deactivate_storage, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_storage", &GridModel::reactivate_storage, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus_storage", &GridModel::change_bus_storage_python, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus_storage", &GridModel::get_bus_storage, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_p_storage", &GridModel::change_p_storage, DocGridModel::_internal_do_not_use.c_str())
        .def("change_q_storage", &GridModel::change_q_storage, DocGridModel::_internal_do_not_use.c_str())

        .def("deactivate_dcline", &GridModel::deactivate_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("reactivate_dcline", &GridModel::reactivate_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_p_dcline", &GridModel::change_p_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_v1_dcline", &GridModel::change_v1_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_v2_dcline", &GridModel::change_v2_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus1_dcline", &GridModel::change_bus1_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("change_bus2_dcline", &GridModel::change_bus2_dcline, DocGridModel::_internal_do_not_use.c_str())
        .def("get_bus1_dcline", &GridModel::get_bus1_dcline, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_bus2_dcline", &GridModel::get_bus2_dcline, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        // get back the results
        // todo _solver
        .def("get_V", &GridModel::get_V, DocGridModel::get_V.c_str())  // this is copied, so no particular return policy is needed
        .def("get_Va", &GridModel::get_Va, DocGridModel::get_Va.c_str())  // this is copied, so no particular return policy is needed
        .def("get_Vm", &GridModel::get_Vm, DocGridModel::get_Vm.c_str())  // this is copied, so no particular return policy is needed
        .def("get_V_solver", &GridModel::get_V_solver, DocGridModel::get_V_solver.c_str(), py::return_value_policy::reference)
        .def("get_Va_solver", &GridModel::get_Va_solver, DocGridModel::get_Va_solver.c_str(), py::return_value_policy::reference)
        .def("get_Vm_solver", &GridModel::get_Vm_solver, DocGridModel::get_Vm_solver.c_str(), py::return_value_policy::reference)
        .def("get_J_solver", &GridModel::get_J_python_solver, DocGridModel::get_J_python_solver.c_str(), py::return_value_policy::reference)

        .def("id_me_to_ac_solver", &GridModel::id_ac_solver_to_me_numpy, DocGridModel::id_me_to_ac_solver.c_str(), py::return_value_policy::reference)
        .def("id_ac_solver_to_me", &GridModel::id_ac_solver_to_me_numpy, DocGridModel::id_ac_solver_to_me.c_str(), py::return_value_policy::reference)
        .def("id_me_to_dc_solver", &GridModel::id_me_to_dc_solver_numpy, DocGridModel::id_me_to_dc_solver.c_str(), py::return_value_policy::reference)
        .def("id_dc_solver_to_me", &GridModel::id_dc_solver_to_me_numpy, DocGridModel::id_dc_solver_to_me.c_str(), py::return_value_policy::reference)
        .def("total_bus", &GridModel::total_bus, DocGridModel::total_bus.c_str())
        .def("nb_connected_bus", &GridModel::nb_connected_bus, DocGridModel::nb_connected_bus.c_str())

        .def("get_pv", &GridModel::get_pv_numpy, DocGridModel::get_pv.c_str(), py::return_value_policy::reference)
        .def("get_pq", &GridModel::get_pq_numpy, DocGridModel::get_pq.c_str(), py::return_value_policy::reference)
        .def("get_slack_ids", &GridModel::get_slack_ids_numpy, DocGridModel::get_slack_ids.c_str(), py::return_value_policy::reference)
        .def("get_slack_ids_dc", &GridModel::get_slack_ids_dc_numpy, DocGridModel::get_slack_ids_dc.c_str(), py::return_value_policy::reference)
        .def("get_slack_weights", &GridModel::get_slack_weights, DocGridModel::get_slack_weights.c_str(), py::return_value_policy::reference)
        .def("get_pv_solver", &GridModel::get_pv_solver_numpy, DocGridModel::get_pv_solver.c_str(), py::return_value_policy::reference)
        .def("get_pq_solver", &GridModel::get_pq_solver_numpy, DocGridModel::get_pq_solver.c_str(), py::return_value_policy::reference)
        .def("get_slack_ids_solver", &GridModel::get_slack_ids_solver_numpy, DocGridModel::get_slack_ids_solver.c_str(), py::return_value_policy::reference)
        .def("get_slack_ids_dc_solver", &GridModel::get_slack_ids_dc_solver_numpy, DocGridModel::get_slack_ids_dc_solver.c_str(), py::return_value_policy::reference)
        .def("get_slack_weights_solver", &GridModel::get_slack_weights_solver, DocGridModel::get_slack_weights_solver.c_str(), py::return_value_policy::reference)

        .def("get_Ybus", &GridModel::get_Ybus, DocGridModel::get_Ybus.c_str())  // this is copied, so no particular return policy is needed
        .def("get_dcYbus", &GridModel::get_dcYbus, DocGridModel::get_dcYbus.c_str())  // this is copied, so no particular return policy is needed
        .def("get_Sbus", &GridModel::get_Sbus, DocGridModel::get_Sbus.c_str())  // this is copied, so no particular return policy is needed
        .def("get_dcSbus", &GridModel::get_dcSbus, DocGridModel::get_dcSbus.c_str())  // this is copied, so no particular return policy is needed
        .def("get_Ybus_solver", &GridModel::get_Ybus_solver, DocGridModel::get_Ybus_solver.c_str(), py::return_value_policy::reference)
        .def("get_dcYbus_solver", &GridModel::get_dcYbus_solver, DocGridModel::get_dcYbus_solver.c_str(), py::return_value_policy::reference)
        .def("get_Sbus_solver", &GridModel::get_Sbus_solver, DocGridModel::get_Sbus_solver.c_str(), py::return_value_policy::reference)
        .def("get_dcSbus_solver", &GridModel::get_dcSbus_solver, DocGridModel::get_dcSbus_solver.c_str(), py::return_value_policy::reference)

        .def("check_solution", &GridModel::check_solution, DocGridModel::check_solution.c_str())

        // TODO optimize that for speed, results are copied apparently
        .def("get_loads_res", &GridModel::get_loads_res, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_loads_status", &GridModel::get_loads_status, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_shunts_res", &GridModel::get_shunts_res, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_shunts_status", &GridModel::get_shunts_status, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_gen_res", &GridModel::get_gen_res, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_gen_status", &GridModel::get_gen_status, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_res1", &GridModel::get_line_res1, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_res2", &GridModel::get_line_res2, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_lines_status", &GridModel::get_lines_status, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res1", &GridModel::get_trafo_res1, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res2", &GridModel::get_trafo_res2, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_status", &GridModel::get_trafo_status, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_storages_res", &GridModel::get_storages_res, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_storages_status", &GridModel::get_storages_status, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_sgens_res", &GridModel::get_sgens_res, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_sgens_status", &GridModel::get_sgens_status, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        .def("get_gen_theta", &GridModel::get_gen_theta, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_load_theta", &GridModel::get_load_theta, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_shunt_theta", &GridModel::get_shunt_theta, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_storage_theta", &GridModel::get_storage_theta, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_theta1", &GridModel::get_line_theta1, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_theta2", &GridModel::get_line_theta2, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_theta1", &GridModel::get_trafo_theta1, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_theta2", &GridModel::get_trafo_theta2, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        .def("get_all_shunt_buses", &GridModel::get_all_shunt_buses_numpy, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_loads_res_full", &GridModel::get_loads_res_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_shunts_res_full", &GridModel::get_shunts_res_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_gen_res_full", &GridModel::get_gen_res_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_res1_full", &GridModel::get_line_res1_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_res2_full", &GridModel::get_line_res2_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res1_full", &GridModel::get_trafo_res1_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res2_full", &GridModel::get_trafo_res2_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_storages_res_full", &GridModel::get_storages_res_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_sgens_res_full", &GridModel::get_sgens_res_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_dcline_res1_full", &GridModel::get_dcline_res1_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_dcline_res2_full", &GridModel::get_dcline_res2_full, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        .def("get_shunt_target_p", &GridModel::get_shunt_target_p, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_load_target_p", &GridModel::get_load_target_p, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_gen_target_p", &GridModel::get_gen_target_p, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_sgen_target_p", &GridModel::get_sgen_target_p, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_storage_target_p", &GridModel::get_storage_target_p, DocGridModel::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        
        // do something with the grid
        .def("deactivate_result_computation", &GridModel::deactivate_result_computation, DocGridModel::deactivate_result_computation.c_str())
        .def("reactivate_result_computation", &GridModel::reactivate_result_computation, DocGridModel::reactivate_result_computation.c_str())
        .def("dc_pf", &GridModel::dc_pf, DocGridModel::dc_pf.c_str())
        .def("ac_pf", &GridModel::ac_pf, DocGridModel::ac_pf.c_str())
        .def("unset_changes", &GridModel::unset_changes, DocGridModel::_internal_do_not_use.c_str())
        .def("tell_recompute_ybus", &GridModel::tell_recompute_ybus, DocGridModel::_internal_do_not_use.c_str())
        .def("tell_recompute_sbus", &GridModel::tell_recompute_sbus, DocGridModel::_internal_do_not_use.c_str())
        .def("tell_solver_need_reset", &GridModel::tell_solver_need_reset, DocGridModel::_internal_do_not_use.c_str())
        .def("tell_ybus_change_sparsity_pattern", &GridModel::tell_ybus_change_sparsity_pattern, DocGridModel::_internal_do_not_use.c_str())
        .def("get_solver_control", &GridModel::get_solver_control, "TODO", py::return_value_policy::reference)
        .def("compute_newton", &GridModel::ac_pf, DocGridModel::ac_pf.c_str())
        .def("get_ptdf", &GridModel::get_ptdf, DocGridModel::get_ptdf.c_str(), py::return_value_policy::reference) 
        .def("get_ptdf_solver", &GridModel::get_ptdf_solver, DocGridModel::get_ptdf_solver.c_str(), py::return_value_policy::reference)
        .def("get_lodf", &GridModel::get_lodf, DocGridModel::get_lodf.c_str(), py::return_value_policy::reference)
        .def("get_Bf", &GridModel::get_Bf, DocGridModel::get_Bf.c_str(), py::return_value_policy::reference)
        .def("get_Bf_solver", &GridModel::get_Bf_solver, DocGridModel::get_Bf_solver.c_str(), py::return_value_policy::reference)

        // apply action faster (optimized for grid2op representation)
        // it is not recommended to use it outside of grid2Op.
        .def("update_gens_p", &GridModel::update_gens_p, DocGridModel::_internal_do_not_use.c_str())
        .def("update_sgens_p", &GridModel::update_sgens_p, DocGridModel::_internal_do_not_use.c_str())
        .def("update_gens_v", &GridModel::update_gens_v, DocGridModel::_internal_do_not_use.c_str())
        .def("update_loads_p", &GridModel::update_loads_p, DocGridModel::_internal_do_not_use.c_str())
        .def("update_loads_q", &GridModel::update_loads_q, DocGridModel::_internal_do_not_use.c_str())
        .def("update_topo", &GridModel::update_topo, DocGridModel::_internal_do_not_use.c_str())
        .def("update_storages_p", &GridModel::update_storages_p, DocGridModel::_internal_do_not_use.c_str())

        // auxiliary functions
        .def("set_n_sub", &GridModel::set_n_sub, DocGridModel::_internal_do_not_use.c_str())
        .def("get_n_sub", &GridModel::get_n_sub, DocGridModel::_internal_do_not_use.c_str())
        .def("set_max_nb_bus_per_sub", &GridModel::set_max_nb_bus_per_sub, DocGridModel::_internal_do_not_use.c_str())
        .def("set_load_pos_topo_vect", &GridModel::set_load_pos_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_gen_pos_topo_vect", &GridModel::set_gen_pos_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_line_pos1_topo_vect", &GridModel::set_line_pos1_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_line_pos2_topo_vect", &GridModel::set_line_pos2_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_trafo_pos1_topo_vect", &GridModel::set_trafo_pos1_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_trafo_pos2_topo_vect", &GridModel::set_trafo_pos2_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_storage_pos_topo_vect", &GridModel::set_storage_pos_topo_vect, DocGridModel::_internal_do_not_use.c_str())
        .def("set_load_to_subid", &GridModel::set_load_to_subid, DocGridModel::_internal_do_not_use.c_str())
        .def("set_gen_to_subid", &GridModel::set_gen_to_subid, DocGridModel::_internal_do_not_use.c_str())
        .def("set_shunt_to_subid", &GridModel::set_shunt_to_subid, DocGridModel::_internal_do_not_use.c_str())
        .def("set_line_to_sub1_id", &GridModel::set_line_to_sub1_id, DocGridModel::_internal_do_not_use.c_str())
        .def("set_line_to_sub2_id", &GridModel::set_line_to_sub2_id, DocGridModel::_internal_do_not_use.c_str())
        .def("set_trafo_to_sub1_id", &GridModel::set_trafo_to_sub1_id, DocGridModel::_internal_do_not_use.c_str())
        .def("set_trafo_to_sub2_id", &GridModel::set_trafo_to_sub2_id, DocGridModel::_internal_do_not_use.c_str())
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
        .def("compute_flows", &TimeSeries::compute_flows, DocComputers::compute_flows.c_str())
        .def("compute_power_flows", &TimeSeries::compute_power_flows, DocComputers::compute_power_flows.c_str())  // need to be done after "compute_Vs"  and "compute_flows"
        
        // results (for now only flow (at each -line origin- or voltages -at each buses)
        // see https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html#returning-values-to-python
        .def("get_flows", &TimeSeries::get_flows, DocComputers::get_flows.c_str(), py::return_value_policy::reference_internal)  // need to be done after "compute_Vs"  and "compute_flows"
        .def("get_power_flows", &TimeSeries::get_power_flows, DocComputers::get_power_flows.c_str(), py::return_value_policy::reference_internal)  // need to be done after "compute_Vs"  and "compute_flows"
        .def("get_voltages", &TimeSeries::get_voltages, DocComputers::get_voltages.c_str(), py::return_value_policy::reference_internal)  // need to be done after "compute_Vs" 
        .def("get_sbuses", &TimeSeries::get_sbuses, DocComputers::get_sbuses.c_str(), py::return_value_policy::reference_internal)  // need to be done after "compute_Vs" 
        ;

    py::class_<ContingencyAnalysis>(m, "ContingencyAnalysisCPP", DocSecurityAnalysis::SecurityAnalysis.c_str())
        .def(py::init<const GridModel &>())
        .def_property("init_from_n_powerflow",
                      &ContingencyAnalysis::get_init_from_n_powerflow,
                      &ContingencyAnalysis::set_init_from_n_powerflow,
                      R"mydelim(Whether to initialize the complex voltages of "
                      "each contingencies with the results of a n-powerflow "
                      "(*ie* a powerflow without any line disconnection) or not. "
                      "Default: false, meaning each simulation is initialized "
                      "with the given input vector)mydelim")

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
        .def("clear_results_only", &ContingencyAnalysis::clear_results_only, DocSecurityAnalysis::clear.c_str())
        .def("close", &ContingencyAnalysis::clear, DocComputers::clear.c_str())
        .def("remove_n1", &ContingencyAnalysis::remove_n1, DocSecurityAnalysis::remove_n1.c_str())
        .def("remove_nk", &ContingencyAnalysis::remove_nk, DocSecurityAnalysis::remove_nk.c_str())
        .def("remove_multiple_n1", &ContingencyAnalysis::remove_multiple_n1, DocSecurityAnalysis::remove_multiple_n1.c_str())
        
        // inspect the class
        .def("my_defaults", &ContingencyAnalysis::my_defaults_vect, DocSecurityAnalysis::my_defaults_vect.c_str())
        .def("is_grid_connected_after_contingency", &ContingencyAnalysis::is_grid_connected_after_contingency, DocGridModel::_internal_do_not_use.c_str())  // TODO

        // perform the computation
        .def("compute", &ContingencyAnalysis::compute, py::call_guard<py::gil_scoped_release>(), DocSecurityAnalysis::compute.c_str())
        .def("compute_flows", &ContingencyAnalysis::compute_flows, DocSecurityAnalysis::compute_flows.c_str())
        .def("compute_power_flows", &ContingencyAnalysis::compute_power_flows, DocSecurityAnalysis::compute_power_flows.c_str())

        // results (for now only flow (at each -line origin- or voltages -at each buses)
        // see https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html#returning-values-to-python
        .def("get_flows", &ContingencyAnalysis::get_flows, DocSecurityAnalysis::get_flows.c_str(), py::return_value_policy::reference_internal)
        .def("get_voltages", &ContingencyAnalysis::get_voltages, DocSecurityAnalysis::get_voltages.c_str(), py::return_value_policy::reference_internal)
        .def("get_power_flows", &ContingencyAnalysis::get_power_flows, DocSecurityAnalysis::get_power_flows.c_str(), py::return_value_policy::reference_internal)

        // timers
        .def("total_time", &ContingencyAnalysis::total_time, DocComputers::total_time.c_str())
        .def("solver_time", &ContingencyAnalysis::solver_time, DocComputers::solver_time.c_str())
        .def("preprocessing_time", &ContingencyAnalysis::preprocessing_time, DocSecurityAnalysis::preprocessing_time.c_str())
        .def("amps_computation_time", &ContingencyAnalysis::amps_computation_time, DocComputers::amps_computation_time.c_str())
        .def("modif_Ybus_time", &ContingencyAnalysis::modif_Ybus_time, DocSecurityAnalysis::modif_Ybus_time.c_str())
        .def("nb_solved", &ContingencyAnalysis::nb_solved, DocComputers::nb_solved.c_str())
        ;
}
