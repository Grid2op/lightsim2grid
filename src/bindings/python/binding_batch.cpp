// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "batch_algorithm/TimeSeries.hpp"
#include "batch_algorithm/ContingencyAnalysis.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

void bind_batch(py::module_& m) {
    py::class_<TimeSeries>(m, "TimeSeriesCPP", DocComputers::Computers.c_str())
        .def(py::init<const GridModel &>())
        .def_property("init_from_n_powerflow",
                      &ContingencyAnalysis::get_init_from_n_powerflow,
                      &ContingencyAnalysis::set_init_from_n_powerflow,
                      R"mydelim(Whether to initialize the complex voltages of "
                      "the first time series with the results of a n-powerflow "
                      "(*ie* a powerflow at the start the simulation) or not. "
                      "Default: false)mydelim")

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
        .def("compute_power_flows", &TimeSeries::compute_power_flows, DocComputers::compute_power_flows.c_str())

        // results
        .def("get_flows", &TimeSeries::get_flows, DocComputers::get_flows.c_str(), py::return_value_policy::reference_internal)
        .def("get_power_flows", &TimeSeries::get_power_flows, DocComputers::get_power_flows.c_str(), py::return_value_policy::reference_internal)
        .def("get_voltages", &TimeSeries::get_voltages, DocComputers::get_voltages.c_str(), py::return_value_policy::reference_internal)
        .def("get_sbuses", &TimeSeries::get_sbuses, DocComputers::get_sbuses.c_str(), py::return_value_policy::reference_internal);

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

        // add contingencies
        .def("add_all_n1", &ContingencyAnalysis::add_all_n1, DocSecurityAnalysis::add_all_n1.c_str())
        .def("add_n1", &ContingencyAnalysis::add_n1, DocSecurityAnalysis::add_n1.c_str())
        .def("add_nk", &ContingencyAnalysis::add_nk, DocSecurityAnalysis::add_nk.c_str())
        .def("add_multiple_n1", &ContingencyAnalysis::add_multiple_n1, DocSecurityAnalysis::add_multiple_n1.c_str())

        // remove contingencies
        .def("reset", &ContingencyAnalysis::clear, DocSecurityAnalysis::clear.c_str())
        .def("clear", &ContingencyAnalysis::clear, DocSecurityAnalysis::clear.c_str())
        .def("clear_results_only", &ContingencyAnalysis::clear_results_only, DocSecurityAnalysis::clear.c_str())
        .def("close", &ContingencyAnalysis::clear, DocComputers::clear.c_str())
        .def("remove_n1", &ContingencyAnalysis::remove_n1, DocSecurityAnalysis::remove_n1.c_str())
        .def("remove_nk", &ContingencyAnalysis::remove_nk, DocSecurityAnalysis::remove_nk.c_str())
        .def("remove_multiple_n1", &ContingencyAnalysis::remove_multiple_n1, DocSecurityAnalysis::remove_multiple_n1.c_str())

        // inspect
        .def("my_defaults", &ContingencyAnalysis::my_defaults_vect, DocSecurityAnalysis::my_defaults_vect.c_str())
        .def("is_grid_connected_after_contingency", &ContingencyAnalysis::is_grid_connected_after_contingency, DocGridModel::_internal_do_not_use.c_str())

        // perform computation
        .def("compute", &ContingencyAnalysis::compute, py::call_guard<py::gil_scoped_release>(), DocSecurityAnalysis::compute.c_str())
        .def("compute_flows", &ContingencyAnalysis::compute_flows, DocSecurityAnalysis::compute_flows.c_str())
        .def("compute_power_flows", &ContingencyAnalysis::compute_power_flows, DocSecurityAnalysis::compute_power_flows.c_str())

        // results
        .def("get_flows", &ContingencyAnalysis::get_flows, DocSecurityAnalysis::get_flows.c_str(), py::return_value_policy::reference_internal)
        .def("get_voltages", &ContingencyAnalysis::get_voltages, DocSecurityAnalysis::get_voltages.c_str(), py::return_value_policy::reference_internal)
        .def("get_power_flows", &ContingencyAnalysis::get_power_flows, DocSecurityAnalysis::get_power_flows.c_str(), py::return_value_policy::reference_internal)

        // timers
        .def("total_time", &ContingencyAnalysis::total_time, DocComputers::total_time.c_str())
        .def("solver_time", &ContingencyAnalysis::solver_time, DocComputers::solver_time.c_str())
        .def("preprocessing_time", &ContingencyAnalysis::preprocessing_time, DocSecurityAnalysis::preprocessing_time.c_str())
        .def("amps_computation_time", &ContingencyAnalysis::amps_computation_time, DocComputers::amps_computation_time.c_str())
        .def("modif_Ybus_time", &ContingencyAnalysis::modif_Ybus_time, DocSecurityAnalysis::modif_Ybus_time.c_str())
        .def("nb_solved", &ContingencyAnalysis::nb_solved, DocComputers::nb_solved.c_str());
}
