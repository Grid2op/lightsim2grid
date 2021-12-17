// Copyright (c) 2020, RTE (https://www.rte-france.com)
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
#include "Computers.h"
#include "SecurityAnalysis.h"

#include "help_fun_msg.h"

namespace py = pybind11;

PYBIND11_MODULE(lightsim2grid_cpp, m)
{

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
    #endif  // KLU_SOLVER_AVAILABLE (or )

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
    #endif  // NICSLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

    py::class_<GaussSeidelSolver>(m, "GaussSeidelSolver", DocSolver::GaussSeidelSolver.c_str())
        .def(py::init<>())
        .def("get_Va", &GaussSeidelSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
        .def("get_Vm", &GaussSeidelSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
        .def("get_V", &GaussSeidelSolver::get_V, DocSolver::get_V.c_str()) 
        .def("get_error", &GaussSeidelSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &GaussSeidelSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
        .def("reset", &GaussSeidelSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
        .def("converged", &GaussSeidelSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
        .def("compute_pf", &GaussSeidelSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
        .def("get_timers", &GaussSeidelSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &GaussSeidelSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization

    py::class_<GaussSeidelSynchSolver>(m, "GaussSeidelSynchSolver", DocSolver::GaussSeidelSynchSolver.c_str())
        .def(py::init<>())
        .def("get_Va", &GaussSeidelSynchSolver::get_Va, DocSolver::get_Va.c_str())  // get the voltage angle vector (vector of double)
        .def("get_Vm", &GaussSeidelSynchSolver::get_Vm, DocSolver::get_Vm.c_str())  // get the voltage magnitude vector (vector of double)
        .def("get_V", &GaussSeidelSynchSolver::get_V, DocSolver::get_V.c_str()) 
        .def("get_error", &GaussSeidelSynchSolver::get_error, DocSolver::get_error.c_str())  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &GaussSeidelSynchSolver::get_nb_iter, DocSolver::get_nb_iter.c_str())  // return the number of iteration performed at the last optimization
        .def("reset", &GaussSeidelSynchSolver::reset, DocSolver::reset.c_str())  // reset the solver to its original state
        .def("converged", &GaussSeidelSynchSolver::converged, DocSolver::converged.c_str())  // whether the solver has converged
        .def("compute_pf", &GaussSeidelSynchSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())  // compute the powerflow
        .def("get_timers", &GaussSeidelSynchSolver::get_timers, DocSolver::get_timers.c_str())  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &GaussSeidelSynchSolver::compute_pf, py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str());  // perform the newton raphson optimization

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
        .def("get_computation_time", &ChooseSolver::get_computation_time, DocSolver::get_computation_time.c_str());

    // iterator for generators
    py::class_<DataGen>(m, "DataGen", DocIterator::DataGen.c_str())
        .def("__len__", [](const DataGen & data) { return data.nb(); })
        .def("__getitem__", [](const DataGen & data, int k){return data[k]; } )
        .def("__iter__", [](const DataGen & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<DataGen::GenInfo>(m, "GenInfo", DocIterator::GenInfo.c_str())
        .def_readonly("id", &DataGen::GenInfo::id, DocIterator::id.c_str())
        .def_readonly("connected", &DataGen::GenInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &DataGen::GenInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("is_slack", &DataGen::GenInfo::is_slack, DocIterator::is_slack.c_str())
        .def_readonly("slack_weight", &DataGen::GenInfo::slack_weight, DocIterator::slack_weight.c_str())
        .def_readonly("target_p_mw", &DataGen::GenInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_vm_pu", &DataGen::GenInfo::target_vm_pu, DocIterator::target_vm_pu.c_str())
        .def_readonly("min_q_mvar", &DataGen::GenInfo::min_q_mvar, DocIterator::min_q_mvar.c_str())
        .def_readonly("max_q_mvar", &DataGen::GenInfo::max_q_mvar, DocIterator::max_q_mvar.c_str())
        .def_readonly("has_res", &DataGen::GenInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &DataGen::GenInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &DataGen::GenInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &DataGen::GenInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &DataGen::GenInfo::res_v_kv, DocIterator::res_v_kv.c_str());

    // iterator for sgens
    py::class_<DataSGen>(m, "DataSGen", DocIterator::DataSGen.c_str())
        .def("__len__", [](const DataSGen & data) { return data.nb(); })
        .def("__getitem__", [](const DataSGen & data, int k){return data[k]; } )
        .def("__iter__", [](const DataSGen & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<DataSGen::SGenInfo>(m, "SGenInfo", DocIterator::SGenInfo.c_str())
        .def_readonly("id", &DataSGen::SGenInfo::id, DocIterator::id.c_str())
        .def_readonly("connected", &DataSGen::SGenInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &DataSGen::SGenInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("min_q_mvar", &DataSGen::SGenInfo::min_q_mvar, DocIterator::min_q_mvar.c_str())
        .def_readonly("max_q_mvar", &DataSGen::SGenInfo::max_q_mvar, DocIterator::max_q_mvar.c_str())
        .def_readonly("min_p_mw", &DataSGen::SGenInfo::min_p_mw, DocIterator::min_p_mw.c_str())
        .def_readonly("max_p_mw", &DataSGen::SGenInfo::max_p_mw, DocIterator::max_p_mw.c_str())
        .def_readonly("target_p_mw", &DataSGen::SGenInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &DataSGen::SGenInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &DataSGen::SGenInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &DataSGen::SGenInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &DataSGen::SGenInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &DataSGen::SGenInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &DataSGen::SGenInfo::res_v_kv, DocIterator::res_v_kv.c_str());

    // iterator for loads (and storage units)
    py::class_<DataLoad>(m, "DataLoad", DocIterator::DataLoad.c_str())
        .def("__len__", [](const DataLoad & data) { return data.nb(); })
        .def("__getitem__", [](const DataLoad & data, int k){return data[k]; } )
        .def("__iter__", [](const DataLoad & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<DataLoad::LoadInfo>(m, "LoadInfo", DocIterator::LoadInfo.c_str())
        .def_readonly("id", &DataLoad::LoadInfo::id, DocIterator::id.c_str())
        .def_readonly("connected", &DataLoad::LoadInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &DataLoad::LoadInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("target_p_mw", &DataLoad::LoadInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &DataLoad::LoadInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &DataLoad::LoadInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &DataLoad::LoadInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &DataLoad::LoadInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &DataLoad::LoadInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &DataLoad::LoadInfo::res_v_kv, DocIterator::res_v_kv.c_str());

    // iterator for shunts
    py::class_<DataShunt>(m, "DataShunt", DocIterator::DataShunt.c_str())
        .def("__len__", [](const DataShunt & data) { return data.nb(); })
        .def("__getitem__", [](const DataShunt & data, int k){return data[k]; } )
        .def("__iter__", [](const DataShunt & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<DataShunt::ShuntInfo>(m, "ShuntInfo", DocIterator::ShuntInfo.c_str())
        .def_readonly("id", &DataShunt::ShuntInfo::id, DocIterator::id.c_str())
        .def_readonly("connected", &DataShunt::ShuntInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &DataShunt::ShuntInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("target_p_mw", &DataShunt::ShuntInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &DataShunt::ShuntInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &DataShunt::ShuntInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &DataShunt::ShuntInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &DataShunt::ShuntInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &DataShunt::ShuntInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &DataShunt::ShuntInfo::res_v_kv, DocIterator::res_v_kv.c_str());

    // iterator for trafos
    py::class_<DataTrafo>(m, "DataTrafo", DocIterator::DataTrafo.c_str())
        .def("__len__", [](const DataTrafo & data) { return data.nb(); })
        .def("__getitem__", [](const DataTrafo & data, int k){return data[k]; } )
        .def("__iter__", [](const DataTrafo & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<DataTrafo::TrafoInfo>(m, "TrafoInfo", DocIterator::TrafoInfo.c_str())
        .def_readonly("id", &DataTrafo::TrafoInfo::id, DocIterator::id.c_str())
        .def_readonly("connected", &DataTrafo::TrafoInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_hv_id", &DataTrafo::TrafoInfo::bus_hv_id, DocIterator::bus_hv_id.c_str())
        .def_readonly("bus_lv_id", &DataTrafo::TrafoInfo::bus_lv_id, DocIterator::bus_lv_id.c_str())
        .def_readonly("r_pu", &DataTrafo::TrafoInfo::r_pu, DocIterator::r_pu.c_str())
        .def_readonly("x_pu", &DataTrafo::TrafoInfo::x_pu, DocIterator::x_pu.c_str())
        .def_readonly("h_pu", &DataTrafo::TrafoInfo::h_pu, DocIterator::h_pu.c_str())
        .def_readonly("is_tap_hv_side", &DataTrafo::TrafoInfo::is_tap_hv_side, DocIterator::is_tap_hv_side.c_str())
        .def_readonly("ratio", &DataTrafo::TrafoInfo::ratio, DocIterator::ratio.c_str())
        .def_readonly("shift_rad", &DataTrafo::TrafoInfo::shift_rad, DocIterator::shift_rad.c_str())
        .def_readonly("has_res", &DataTrafo::TrafoInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_hv_mw", &DataTrafo::TrafoInfo::res_p_hv_mw, DocIterator::res_p_hv_mw.c_str())
        .def_readonly("res_q_hv_mvar", &DataTrafo::TrafoInfo::res_q_hv_mvar, DocIterator::res_q_hv_mvar.c_str())
        .def_readonly("res_v_hv_kv", &DataTrafo::TrafoInfo::res_v_hv_kv, DocIterator::res_v_hv_kv.c_str())
        .def_readonly("res_a_hv_ka", &DataTrafo::TrafoInfo::res_a_hv_ka, DocIterator::res_a_hv_ka.c_str())
        .def_readonly("res_p_lv_mw", &DataTrafo::TrafoInfo::res_p_lv_mw, DocIterator::res_p_lv_mw.c_str())
        .def_readonly("res_q_lv_mvar", &DataTrafo::TrafoInfo::res_q_lv_mvar, DocIterator::res_q_lv_mvar.c_str())
        .def_readonly("res_v_lv_kv", &DataTrafo::TrafoInfo::res_v_lv_kv, DocIterator::res_v_lv_kv.c_str())
        .def_readonly("res_a_lv_ka", &DataTrafo::TrafoInfo::res_a_lv_ka, DocIterator::res_a_lv_ka.c_str())
        .def_readonly("res_theta_hv_deg", &DataTrafo::TrafoInfo::res_theta_hv_deg, DocIterator::res_theta_hv_deg.c_str())
        .def_readonly("res_theta_lv_deg", &DataTrafo::TrafoInfo::res_theta_lv_deg, DocIterator::res_theta_lv_deg.c_str());

    // iterator for trafos
    py::class_<DataLine>(m, "DataLine", DocIterator::DataLine.c_str())
        .def("__len__", [](const DataLine & data) { return data.nb(); })
        .def("__getitem__", [](const DataLine & data, int k){return data[k]; } )
        .def("__iter__", [](const DataLine & data) {
       return py::make_iterator(data.begin(), data.end());
    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    py::class_<DataLine::LineInfo>(m, "LineInfo", DocIterator::LineInfo.c_str())
        .def_readonly("id", &DataLine::LineInfo::id, DocIterator::id.c_str())
        .def_readonly("connected", &DataLine::LineInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_or_id", &DataLine::LineInfo::bus_or_id, DocIterator::bus_or_id.c_str())
        .def_readonly("bus_ex_id", &DataLine::LineInfo::bus_ex_id, DocIterator::bus_ex_id.c_str())
        .def_readonly("r_pu", &DataLine::LineInfo::r_pu, DocIterator::r_pu.c_str())
        .def_readonly("x_pu", &DataLine::LineInfo::x_pu, DocIterator::x_pu.c_str())
        .def_readonly("h_pu", &DataLine::LineInfo::h_pu, DocIterator::h_pu.c_str())
        .def_readonly("has_res", &DataLine::LineInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_or_mw", &DataLine::LineInfo::res_p_or_mw, DocIterator::res_p_or_mw.c_str())
        .def_readonly("res_q_or_mvar", &DataLine::LineInfo::res_q_or_mvar, DocIterator::res_q_or_mvar.c_str())
        .def_readonly("res_v_or_kv", &DataLine::LineInfo::res_v_or_kv, DocIterator::res_v_or_kv.c_str())
        .def_readonly("res_a_or_ka", &DataLine::LineInfo::res_a_or_ka, DocIterator::res_a_or_ka.c_str())
        .def_readonly("res_p_ex_mw", &DataLine::LineInfo::res_p_ex_mw, DocIterator::res_p_ex_mw.c_str())
        .def_readonly("res_q_ex_mvar", &DataLine::LineInfo::res_q_ex_mvar, DocIterator::res_q_ex_mvar.c_str())
        .def_readonly("res_v_ex_kv", &DataLine::LineInfo::res_v_ex_kv, DocIterator::res_v_ex_kv.c_str())
        .def_readonly("res_a_ex_ka", &DataLine::LineInfo::res_a_ex_ka, DocIterator::res_a_ex_ka.c_str())
        .def_readonly("res_theta_or_deg", &DataLine::LineInfo::res_theta_or_deg, DocIterator::res_theta_or_deg.c_str())
        .def_readonly("res_theta_ex_deg", &DataLine::LineInfo::res_theta_ex_deg, DocIterator::res_theta_ex_deg.c_str());

    // converters
    py::class_<PandaPowerConverter>(m, "PandaPowerConverter")
        .def(py::init<>())
        .def("set_f_hz", &PandaPowerConverter::set_f_hz)
        .def("set_sn_mva", &PandaPowerConverter::set_sn_mva)
        .def("get_line_param", &PandaPowerConverter::get_line_param)
        .def("get_trafo_param", &PandaPowerConverter::get_trafo_param);

    py::class_<GridModel>(m, "GridModel", DocGridModel::GridModel.c_str())
        .def(py::init<>())
        .def("copy", &GridModel::copy)

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
        .def("set_init_vm_pu", &GridModel::set_init_vm_pu, DocGridModel::_internal_do_not_use.c_str())  // TODO use python "property" for that
        .def("get_init_vm_pu", &GridModel::get_init_vm_pu, DocGridModel::_internal_do_not_use.c_str())
        .def("set_sn_mva", &GridModel::set_sn_mva, DocGridModel::_internal_do_not_use.c_str())   // TODO use python "property" for that
        .def("get_sn_mva", &GridModel::get_sn_mva, DocGridModel::_internal_do_not_use.c_str())

        // init its elements
        .def("init_powerlines", &GridModel::init_powerlines, DocGridModel::_internal_do_not_use.c_str())  // TODO code the possibility to add / remove a powerline after creation
        .def("init_shunt", &GridModel::init_shunt, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_trafo", &GridModel::init_trafo, DocGridModel::_internal_do_not_use.c_str())  // same 
        .def("init_generators", &GridModel::init_generators, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_loads", &GridModel::init_loads, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_storages", &GridModel::init_storages, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("init_sgens", &GridModel::init_sgens, DocGridModel::_internal_do_not_use.c_str())  // same
        .def("add_gen_slackbus", &GridModel::add_gen_slackbus, DocGridModel::_internal_do_not_use.c_str()) // same
        .def("remove_gen_slackbus", &GridModel::remove_gen_slackbus, DocGridModel::_internal_do_not_use.c_str())  // same

        // inspect the grid
        .def("get_lines", &GridModel::get_lines, DocGridModel::get_lines.c_str())
        .def("get_trafos", &GridModel::get_trafos, DocGridModel::get_trafos.c_str())
        .def("get_generators", &GridModel::get_generators, DocGridModel::get_generators.c_str())
        .def("get_static_generators", &GridModel::get_static_generators, DocGridModel::get_static_generators.c_str())
        .def("get_shunts", &GridModel::get_shunts, DocGridModel::get_shunts.c_str())
        .def("get_storages", &GridModel::get_storages, DocGridModel::get_storages.c_str())
        .def("get_loads", &GridModel::get_loads, DocGridModel::get_loads.c_str())

        // modify the grid
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

        // do something with the grid
        // .def("init_Ybus", &DataModel::init_Ybus) // temporary
        .def("deactivate_result_computation", &GridModel::deactivate_result_computation, DocGridModel::deactivate_result_computation.c_str())
        .def("reactivate_result_computation", &GridModel::reactivate_result_computation, DocGridModel::reactivate_result_computation.c_str())
        .def("dc_pf", &GridModel::dc_pf, DocGridModel::dc_pf.c_str())
        .def("ac_pf", &GridModel::ac_pf, DocGridModel::ac_pf.c_str())
        .def("unset_topo_changed", &GridModel::unset_topo_changed, DocGridModel::_internal_do_not_use.c_str())
        .def("tell_topo_changed", &GridModel::tell_topo_changed, DocGridModel::_internal_do_not_use.c_str())
        .def("compute_newton", &GridModel::ac_pf, DocGridModel::ac_pf.c_str())

         // apply action faster (optimized for grid2op representation)
         // it is not recommended to use it outside of grid2Op.
        .def("update_bus_status", &GridModel::update_bus_status, DocGridModel::_internal_do_not_use.c_str())
        .def("update_gens_p", &GridModel::update_gens_p, DocGridModel::_internal_do_not_use.c_str())
        .def("update_gens_v", &GridModel::update_gens_v, DocGridModel::_internal_do_not_use.c_str())
        .def("update_loads_p", &GridModel::update_loads_p, DocGridModel::_internal_do_not_use.c_str())
        .def("update_loads_q", &GridModel::update_loads_q, DocGridModel::_internal_do_not_use.c_str())
        .def("update_topo", &GridModel::update_topo, DocGridModel::_internal_do_not_use.c_str())
        .def("update_storages_p", &GridModel::update_storages_p, DocGridModel::_internal_do_not_use.c_str())

        // auxiliary functions
        .def("set_n_sub", &GridModel::set_n_sub, DocGridModel::_internal_do_not_use.c_str())
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
        ;

    py::class_<Computers>(m, "Computers", DocComputers::Computers.c_str())
        .def(py::init<const GridModel &>())

        // solver control
        .def("change_solver", &Computers::change_solver, DocGridModel::change_solver.c_str())
        .def("available_solvers", &Computers::available_solvers, DocGridModel::available_solvers.c_str())
        .def("get_solver_type", &Computers::get_solver_type, DocGridModel::get_solver_type.c_str())

        // timers
        .def("total_time", &Computers::total_time, DocComputers::total_time.c_str())
        .def("solver_time", &Computers::solver_time, DocComputers::solver_time.c_str())
        .def("preprocessing_time", &Computers::preprocessing_time, DocComputers::preprocessing_time.c_str())
        .def("amps_computation_time", &Computers::amps_computation_time, DocComputers::amps_computation_time.c_str())
        .def("nb_solved", &Computers::nb_solved, DocComputers::nb_solved.c_str())

        // status
        .def("get_status", &Computers::get_status, DocComputers::get_status.c_str())

        // perform the computations
        .def("compute_Vs", &Computers::compute_Vs, py::call_guard<py::gil_scoped_release>(), DocComputers::compute_Vs.c_str())
        .def("compute_flows", &Computers::compute_flows, py::call_guard<py::gil_scoped_release>(), DocComputers::compute_flows.c_str())
        .def("compute_power_flows", &Computers::compute_power_flows, DocComputers::compute_power_flows.c_str())  // need to be done after "compute_Vs"  and "compute_flows"
        
        // results (for now only flow (at each -line origin- or voltages -at each buses)
        .def("get_flows", &Computers::get_flows, DocComputers::get_flows.c_str())  // need to be done after "compute_Vs"  and "compute_flows"
        .def("get_power_flows", &Computers::get_power_flows, DocComputers::get_power_flows.c_str())  // need to be done after "compute_Vs"  and "compute_flows"
        .def("get_voltages", &Computers::get_voltages, DocComputers::get_voltages.c_str())  // need to be done after "compute_Vs" 
        .def("get_sbuses", &Computers::get_sbuses, DocComputers::get_sbuses.c_str())  // need to be done after "compute_Vs" 
        ;

    py::class_<SecurityAnalysis>(m, "SecurityAnalysisCPP", DocSecurityAnalysis::SecurityAnalysis.c_str())
        .def(py::init<const GridModel &>())
        // solver control
        .def("change_solver", &Computers::change_solver, DocGridModel::change_solver.c_str())
        .def("available_solvers", &Computers::available_solvers, DocGridModel::available_solvers.c_str())
        .def("get_solver_type", &Computers::get_solver_type, DocGridModel::get_solver_type.c_str())

        // add some defaults
        .def("add_all_n1", &SecurityAnalysis::add_all_n1, DocSecurityAnalysis::add_all_n1.c_str())
        .def("add_n1", &SecurityAnalysis::add_n1, DocSecurityAnalysis::add_n1.c_str())
        .def("add_nk", &SecurityAnalysis::add_nk, DocSecurityAnalysis::add_nk.c_str())
        .def("add_multiple_n1", &SecurityAnalysis::add_multiple_n1, DocSecurityAnalysis::add_multiple_n1.c_str())

        // remove some defaults (TODO)
        .def("reset", &SecurityAnalysis::clear, DocSecurityAnalysis::clear.c_str())
        .def("clear", &SecurityAnalysis::clear, DocSecurityAnalysis::clear.c_str())
        .def("remove_n1", &SecurityAnalysis::remove_n1, DocSecurityAnalysis::remove_n1.c_str())
        .def("remove_nk", &SecurityAnalysis::remove_nk, DocSecurityAnalysis::remove_nk.c_str())
        .def("remove_multiple_n1", &SecurityAnalysis::remove_multiple_n1, DocSecurityAnalysis::remove_multiple_n1.c_str())
        
        // inspect the class
        .def("my_defaults", &SecurityAnalysis::my_defaults_vect, DocSecurityAnalysis::my_defaults_vect.c_str())

        // perform the computation
        .def("compute", &SecurityAnalysis::compute, py::call_guard<py::gil_scoped_release>(), DocSecurityAnalysis::compute.c_str())
        .def("compute_flows", &SecurityAnalysis::compute_flows, py::call_guard<py::gil_scoped_release>(), DocSecurityAnalysis::compute_flows.c_str())
        .def("compute_power_flows", &SecurityAnalysis::compute_power_flows, DocSecurityAnalysis::compute_power_flows.c_str())

        // results (for now only flow (at each -line origin- or voltages -at each buses)
        .def("get_flows", &SecurityAnalysis::get_flows, DocSecurityAnalysis::get_flows.c_str())
        .def("get_voltages", &SecurityAnalysis::get_voltages, DocSecurityAnalysis::get_voltages.c_str())
        .def("get_power_flows", &SecurityAnalysis::get_power_flows, DocSecurityAnalysis::get_power_flows.c_str())

        // timers
        .def("total_time", &SecurityAnalysis::total_time, DocComputers::total_time.c_str())
        .def("solver_time", &SecurityAnalysis::solver_time, DocComputers::solver_time.c_str())
        .def("preprocessing_time", &SecurityAnalysis::preprocessing_time, DocSecurityAnalysis::preprocessing_time.c_str())
        .def("amps_computation_time", &SecurityAnalysis::amps_computation_time, DocComputers::amps_computation_time.c_str())
        .def("modif_Ybus_time", &SecurityAnalysis::modif_Ybus_time, DocSecurityAnalysis::modif_Ybus_time.c_str())
        .def("nb_solved", &SecurityAnalysis::nb_solved, DocComputers::nb_solved.c_str())
        ;
}
