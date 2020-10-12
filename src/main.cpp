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

#include "KLUSolver.h"
#include "SparseLUSolver.h"
#include "GaussSeidelSolver.h"
#include "DataConverter.h"
#include "GridModel.h"

namespace py = pybind11;

PYBIND11_MODULE(lightsim2grid_cpp, m)
{

    // solvers
    py::enum_<SolverType>(m, "SolverType")
        .value("SparseLU", SolverType::SparseLU)
        .value("KLU", SolverType::KLU)
        .value("GaussSeidel", SolverType::GaussSeidel)
        .value("DC", SolverType::DC)
        .export_values();

    #ifdef KLU_SOLVER_AVAILABLE
    py::class_<KLUSolver>(m, "KLUSolver")
        .def(py::init<>())
        .def("get_J", &KLUSolver::get_J)  // (get the jacobian matrix, sparse csc matrix)
        .def("get_Va", &KLUSolver::get_Va)  // get the voltage angle vector (vector of double)
        .def("get_Vm", &KLUSolver::get_Vm)  // get the voltage magnitude vector (vector of double)
        .def("get_error", &KLUSolver::get_error)  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &KLUSolver::get_nb_iter)  // return the number of iteration performed at the last optimization
        .def("reset", &KLUSolver::reset)  // reset the solver to its original state
        .def("converged", &KLUSolver::converged)  // whether the solver has converged
        .def("compute_pf", &KLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>())  // perform the newton raphson optimization
        .def("get_timers", &KLUSolver::get_timers)  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &KLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>() );  // perform the newton raphson optimization
    #endif

    py::class_<SparseLUSolver>(m, "SparseLUSolver")
        .def(py::init<>())
        .def("get_J", &SparseLUSolver::get_J)  // (get the jacobian matrix, sparse csc matrix)
        .def("get_Va", &SparseLUSolver::get_Va)  // get the voltage angle vector (vector of double)
        .def("get_Vm", &SparseLUSolver::get_Vm)  // get the voltage magnitude vector (vector of double)
        .def("get_error", &SparseLUSolver::get_error)  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &SparseLUSolver::get_nb_iter)  // return the number of iteration performed at the last optimization
        .def("reset", &SparseLUSolver::reset)  // reset the solver to its original state
        .def("converged", &SparseLUSolver::converged)  // whether the solver has converged
        .def("compute_pf", &SparseLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>())  // perform the newton raphson optimization
        .def("get_timers", &SparseLUSolver::get_timers)  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &SparseLUSolver::compute_pf, py::call_guard<py::gil_scoped_release>() );  // perform the newton raphson optimization

    py::class_<GaussSeidelSolver>(m, "GaussSeidelSolver")
        .def(py::init<>())
        .def("get_Va", &GaussSeidelSolver::get_Va)  // get the voltage angle vector (vector of double)
        .def("get_Vm", &GaussSeidelSolver::get_Vm)  // get the voltage magnitude vector (vector of double)
        .def("get_error", &GaussSeidelSolver::get_error)  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &GaussSeidelSolver::get_nb_iter)  // return the number of iteration performed at the last optimization
        .def("reset", &GaussSeidelSolver::reset)  // reset the solver to its original state
        .def("converged", &GaussSeidelSolver::converged)  // whether the solver has converged
        .def("compute_pf", &GaussSeidelSolver::compute_pf, py::call_guard<py::gil_scoped_release>())  // compute the powerflow
        .def("get_timers", &GaussSeidelSolver::get_timers)  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &GaussSeidelSolver::compute_pf, py::call_guard<py::gil_scoped_release>() );  // perform the newton raphson optimization

    py::class_<DCSolver>(m, "DCSolver")
        .def(py::init<>())
        .def("get_Va", &DCSolver::get_Va)  // get the voltage angle vector (vector of double)
        .def("get_Vm", &DCSolver::get_Vm)  // get the voltage magnitude vector (vector of double)
        .def("get_error", &DCSolver::get_error)  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &DCSolver::get_nb_iter)  // return the number of iteration performed at the last optimization
        .def("reset", &DCSolver::reset)  // reset the solver to its original state
        .def("converged", &DCSolver::converged)  // whether the solver has converged
        .def("compute_pf", &DCSolver::compute_pf, py::call_guard<py::gil_scoped_release>())  // compute the powerflow
        .def("get_timers", &DCSolver::get_timers)  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &DCSolver::compute_pf, py::call_guard<py::gil_scoped_release>() );  // perform the newton raphson optimization


    // converters
    py::class_<PandaPowerConverter>(m, "PandaPowerConverter")
        .def(py::init<>())
        .def("set_f_hz", &PandaPowerConverter::set_f_hz)
        .def("set_sn_mva", &PandaPowerConverter::set_sn_mva)
        .def("get_line_param", &PandaPowerConverter::get_line_param)
        .def("get_trafo_param", &PandaPowerConverter::get_trafo_param);

    py::class_<GridModel>(m, "GridModel")
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
        .def("change_solver", &GridModel::change_solver)  // change the solver to use (KLU - faster or SparseLU - available everywhere)
        .def("available_solvers", &GridModel::available_solvers)  // retrieve the solver available for your installation
        .def("get_computation_time", &GridModel::get_computation_time)  // get the computation time spent in the solver
        .def("get_solver_type", &GridModel::get_solver_type)  // get the type of solver used

        // init the grid
        .def("init_bus", &GridModel::init_bus)
        .def("init_powerlines", &GridModel::init_powerlines)
        .def("init_shunt", &GridModel::init_shunt)
        .def("init_trafo", &GridModel::init_trafo)
        .def("init_generators", &GridModel::init_generators)
        .def("init_loads", &GridModel::init_loads)
        .def("add_gen_slackbus", &GridModel::add_gen_slackbus)

        // modify the grid
        .def("deactivate_bus", &GridModel::deactivate_bus)
        .def("reactivate_bus", &GridModel::reactivate_bus)
        .def("nb_bus", &GridModel::nb_bus)

        .def("deactivate_powerline", &GridModel::deactivate_powerline)
        .def("reactivate_powerline", &GridModel::reactivate_powerline)
        .def("change_bus_powerline_or", &GridModel::change_bus_powerline_or)
        .def("change_bus_powerline_ex", &GridModel::change_bus_powerline_ex)
        .def("get_bus_powerline_or", &GridModel::get_bus_powerline_or)
        .def("get_bus_powerline_ex", &GridModel::get_bus_powerline_ex)

        .def("deactivate_trafo", &GridModel::deactivate_trafo)
        .def("reactivate_trafo", &GridModel::reactivate_trafo)
        .def("change_bus_trafo_hv", &GridModel::change_bus_trafo_hv)
        .def("change_bus_trafo_lv", &GridModel::change_bus_trafo_lv)
        .def("get_bus_trafo_hv", &GridModel::get_bus_trafo_hv)
        .def("get_bus_trafo_lv", &GridModel::get_bus_trafo_lv)

        .def("deactivate_load", &GridModel::deactivate_load)
        .def("reactivate_load", &GridModel::reactivate_load)
        .def("change_bus_load", &GridModel::change_bus_load)
        .def("get_bus_load", &GridModel::get_bus_load)
        .def("change_p_load", &GridModel::change_p_load)
        .def("change_q_load", &GridModel::change_q_load)

        .def("deactivate_gen", &GridModel::deactivate_gen)
        .def("reactivate_gen", &GridModel::reactivate_gen)
        .def("change_bus_gen", &GridModel::change_bus_gen)
        .def("get_bus_gen", &GridModel::get_bus_gen)
        .def("change_p_gen", &GridModel::change_p_gen)
        .def("change_v_gen", &GridModel::change_v_gen)

        .def("deactivate_shunt", &GridModel::deactivate_shunt)
        .def("reactivate_shunt", &GridModel::reactivate_shunt)
        .def("change_bus_shunt", &GridModel::change_bus_shunt)
        .def("get_bus_shunt", &GridModel::get_bus_shunt)
        .def("change_p_shunt", &GridModel::change_p_shunt)
        .def("change_q_shunt", &GridModel::change_q_shunt)

        // get back the results
        .def("get_Va", &GridModel::get_Va)
        .def("get_Vm", &GridModel::get_Vm)

        .def("get_loads_res", &GridModel::get_loads_res)
        .def("get_loads_status", &GridModel::get_loads_status)
        .def("get_shunts_res", &GridModel::get_shunts_res)
        .def("get_shunts_status", &GridModel::get_shunts_status)
        .def("get_gen_res", &GridModel::get_gen_res)
        .def("get_gen_status", &GridModel::get_gen_status)
        .def("get_lineor_res", &GridModel::get_lineor_res)
        .def("get_lineex_res", &GridModel::get_lineex_res)
        .def("get_lines_status", &GridModel::get_lines_status)
        .def("get_trafohv_res", &GridModel::get_trafohv_res)
        .def("get_trafolv_res", &GridModel::get_trafolv_res)
        .def("get_trafo_status", &GridModel::get_trafo_status)

        // do something with the grid
        // .def("init_Ybus", &DataModel::init_Ybus) // temporary
        .def("get_Ybus", &GridModel::get_Ybus)
        .def("get_Sbus", &GridModel::get_Sbus)
        .def("get_pv", &GridModel::get_pv)
        .def("get_pq", &GridModel::get_pq)

        .def("deactivate_result_computation", &GridModel::deactivate_result_computation)
        .def("reactivate_result_computation", &GridModel::reactivate_result_computation)
        .def("dc_pf", &GridModel::dc_pf)
        .def("dc_pf_old", &GridModel::dc_pf_old)
        .def("ac_pf", &GridModel::ac_pf)
        .def("compute_newton", &GridModel::ac_pf)

         // apply action faster (optimized for grid2op representation)
         // it is not recommended to use it outside of grid2Op.
        .def("update_bus_status", &GridModel::update_bus_status)
        .def("update_gens_p", &GridModel::update_gens_p)
        .def("update_gens_v", &GridModel::update_gens_v)
        .def("update_loads_p", &GridModel::update_loads_p)
        .def("update_loads_q", &GridModel::update_loads_q)
        .def("update_topo", &GridModel::update_topo)
        // auxiliary functions
        .def("set_n_sub", &GridModel::set_n_sub)
        .def("set_load_pos_topo_vect", &GridModel::set_load_pos_topo_vect)
        .def("set_gen_pos_topo_vect", &GridModel::set_gen_pos_topo_vect)
        .def("set_line_or_pos_topo_vect", &GridModel::set_line_or_pos_topo_vect)
        .def("set_line_ex_pos_topo_vect", &GridModel::set_line_ex_pos_topo_vect)
        .def("set_trafo_hv_pos_topo_vect", &GridModel::set_trafo_hv_pos_topo_vect)
        .def("set_trafo_lv_pos_topo_vect", &GridModel::set_trafo_lv_pos_topo_vect)
        .def("set_load_to_subid", &GridModel::set_load_to_subid)
        .def("set_gen_to_subid", &GridModel::set_gen_to_subid)
        .def("set_line_or_to_subid", &GridModel::set_line_or_to_subid)
        .def("set_line_ex_to_subid", &GridModel::set_line_ex_to_subid)
        .def("set_trafo_hv_to_subid", &GridModel::set_trafo_hv_to_subid)
        .def("set_trafo_lv_to_subid", &GridModel::set_trafo_lv_to_subid)
        ;

}