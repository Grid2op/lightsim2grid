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
#include "batch_algorithm/LimitViolation.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

void bind_batch(py::module_& m) {
    py::enum_<ViolationElementType>(m, "ViolationElementType", DocContingencyAnalysis::ViolationElementType.c_str())
        .value("BUS", ViolationElementType::BUS)
        .value("LINE", ViolationElementType::LINE)
        .value("TRAFO", ViolationElementType::TRAFO)
        .value("GRID", ViolationElementType::GRID,
               "The whole grid / contingency, not a specific element (see LimitViolationType.NOT_SIMULATED "
               "/ LimitViolationType.DIVERGENCE).");

    py::enum_<LimitViolationType>(m, "LimitViolationType", DocContingencyAnalysis::LimitViolationType.c_str())
        .value("LOW_VOLTAGE", LimitViolationType::LOW_VOLTAGE)
        .value("HIGH_VOLTAGE", LimitViolationType::HIGH_VOLTAGE)
        .value("CURRENT", LimitViolationType::CURRENT)
        .value("NOT_SIMULATED", LimitViolationType::NOT_SIMULATED,
               "A pre-check (graph connectivity) skipped this contingency: the solver was never "
               "invoked (element_type is ViolationElementType.GRID).")
        .value("DIVERGENCE", LimitViolationType::DIVERGENCE,
               "The solver was invoked for this contingency but did not converge (element_type is "
               "ViolationElementType.GRID).");

    py::class_<LimitViolation>(m, "LimitViolation", DocContingencyAnalysis::LimitViolation.c_str())
        .def_readonly("element_type", &LimitViolation::element_type, DocContingencyAnalysis::element_type.c_str())
        .def_readonly("element_id", &LimitViolation::element_id, DocContingencyAnalysis::element_id.c_str())
        .def_readonly("side", &LimitViolation::side, DocContingencyAnalysis::side.c_str())
        .def_readonly("violation_type", &LimitViolation::violation_type, DocContingencyAnalysis::violation_type.c_str())
        .def_readonly("value", &LimitViolation::value, DocContingencyAnalysis::value.c_str())
        .def_readonly("limit", &LimitViolation::limit, DocContingencyAnalysis::limit.c_str())
        .def_readonly("name", &LimitViolation::name, DocContingencyAnalysis::violation_name.c_str());

    py::class_<TimeSeries>(m, "TimeSeriesCPP", DocTimeSeries::TimeSeries.c_str())
        .def(py::init<const LSGrid &>())
        .def_property("init_from_n_powerflow",
                      [](const TimeSeries & self){ return self.get_init_from_n_powerflow(); },
                      [](TimeSeries & self, bool val){ self.set_init_from_n_powerflow(val); },
                      "Whether to initialize the complex voltages of "
                      "the first time series with the results of a n-powerflow "
                      "(*ie* a powerflow at the start the simulation) or not. "
                      "Default: false")

        // solver control
        .def("change_algorithm", &TimeSeries::change_algorithm, DocLSGrid::change_algorithm.c_str())
        .def("change_solver", &TimeSeries::change_algorithm, "DEPRECATED: use 'change_algorithm' instead")
        .def("available_default_algorithms", &TimeSeries::available_default_algorithms, DocLSGrid::available_default_algorithms.c_str())
        .def("get_algo_type", &TimeSeries::get_algo_type, DocLSGrid::get_algo_type.c_str())
        .def("get_algo_config", &TimeSeries::get_algo_config,
             "Config (eg ScalingPolicyType / damping parameters) of the internal solver used for "
             "every step. Copied once from the grid model's own get_ac_algo_config() at "
             "construction time, then independent of it; re-apply with set_algo_config() if you "
             "change the grid model's config afterwards, or after change_algorithm().")
        .def("set_algo_config", &TimeSeries::set_algo_config, py::arg("config"),
             "See get_algo_config().")

        // timers
        .def("total_time", &TimeSeries::total_time, DocTimeSeries::total_time.c_str())
        .def("solver_time", &TimeSeries::solver_time, DocTimeSeries::solver_time.c_str())
        .def("preprocessing_time", &TimeSeries::preprocessing_time, DocTimeSeries::preprocessing_time.c_str())
        .def("amps_computation_time", &TimeSeries::amps_computation_time, DocTimeSeries::amps_computation_time.c_str())
        .def("nb_solved", &TimeSeries::nb_solved, DocTimeSeries::nb_solved.c_str())

        // status
        .def("get_status", &TimeSeries::get_status, DocTimeSeries::get_status.c_str())
        .def("clear", &TimeSeries::clear, DocTimeSeries::clear.c_str())
        .def("close", &TimeSeries::clear, DocTimeSeries::clear.c_str())

        // perform the computations
        .def("compute_Vs", &TimeSeries::compute_Vs, py::call_guard<py::gil_scoped_release>(), DocTimeSeries::compute_Vs.c_str())
        .def("compute_flows", &TimeSeries::compute_flows, DocTimeSeries::compute_flows.c_str())
        .def("compute_power_flows", &TimeSeries::compute_power_flows, DocTimeSeries::compute_power_flows.c_str())

        // results
        .def("get_flows", &TimeSeries::get_flows, DocTimeSeries::get_flows.c_str(), py::return_value_policy::reference_internal)
        .def("get_power_flows", &TimeSeries::get_power_flows, DocTimeSeries::get_power_flows.c_str(), py::return_value_policy::reference_internal)
        .def("get_voltages", &TimeSeries::get_voltages, DocTimeSeries::get_voltages.c_str(), py::return_value_policy::reference_internal)
        .def("get_sbuses", &TimeSeries::get_sbuses, DocTimeSeries::get_sbuses.c_str(), py::return_value_policy::reference_internal);

    py::class_<ContingencyAnalysis>(m, "ContingencyAnalysisCPP", DocContingencyAnalysis::ContingencyAnalysis.c_str())
        .def(py::init<const LSGrid &, bool>(), py::arg("grid_model"), py::arg("compute_limit_violations") = false)
        .def_property("compute_limit_violations",
                      [](const ContingencyAnalysis & self){ return self.get_compute_limit_violations(); },
                      [](ContingencyAnalysis & self, bool val){ self.set_compute_limit_violations(val); },
                      "Whether limit violations are computed inline, per contingency, "
                      "during compute() (see converged / get_violations / converged_n / "
                      "get_violations_n). Defaults to ``False``. Computing violations means an extra "
                      "per-element current / voltage check in every contingency's solve, so "
                      "users who only need compute_flows() / get_flows() should leave this off. "
                      "Changing this flag clears any previously-computed results.")
        .def_property("init_from_n_powerflow",
                      [](const ContingencyAnalysis & self){ return self.get_init_from_n_powerflow(); },
                      [](ContingencyAnalysis & self, bool val){ self.set_init_from_n_powerflow(val); },
                      "Whether to initialize the complex voltages of "
                      "each contingencies with the results of a n-powerflow "
                      "(*ie* a powerflow without any line disconnection) or not. "
                      "Default: false, meaning each simulation is initialized "
                      "with the given input vector")
        .def_property("handle_disconnected_grid",
                      [](const ContingencyAnalysis & self){ return self.get_handle_disconnected_grid(); },
                      [](ContingencyAnalysis & self, bool val){ self.set_handle_disconnected_grid(val); },
                      "Whether to simulate the contingencies that split the grid in "
                      "multiple connected components. When False (default) such contingencies "
                      "are skipped (their voltages are left at 0), reproducing the legacy "
                      "behaviour. When True, the largest connected component is solved while "
                      "the buses of the other component(s) are masked (their voltage is "
                      "reported as 0). Supported by the Newton-Raphson family (AC) and the DC "
                      "solver; a non Newton-Raphson AC algorithm is rejected.")
        .def_property("nb_thread",
                      [](const ContingencyAnalysis & self){ return self.get_nb_thread(); },
                      [](ContingencyAnalysis & self, int val){ self.set_nb_thread(val); },
                      "Number of OS threads used to solve the contingencies (default ``1``). "
                      "With nb_thread == 1 the behaviour is identical to the legacy sequential "
                      "computation. With nb_thread > 1 the contingency list is split into "
                      "contiguous ranges, each solved by its own thread (each with its own solver "
                      "and admittance matrix copy), writing to disjoint rows of the result matrix. "
                      "The results do not depend on the number of threads. Values < 1 are "
                      "clamped to 1.")

        // solver control
        .def("change_algorithm", &ContingencyAnalysis::change_algorithm, DocLSGrid::change_algorithm.c_str())
        .def("change_solver", &ContingencyAnalysis::change_algorithm, "DEPRECATED: use 'change_algorithm' instead")
        .def("available_default_algorithms", &ContingencyAnalysis::available_default_algorithms, DocLSGrid::available_default_algorithms.c_str())
        .def("get_algo_type", &ContingencyAnalysis::get_algo_type, DocLSGrid::get_algo_type.c_str())
        .def("get_algo_config", &ContingencyAnalysis::get_algo_config,
             "Config (eg ScalingPolicyType / damping parameters) of the internal solver used for "
             "the pre-contingency ('n') and every per-contingency powerflow. Copied once from the "
             "grid model's own get_ac_algo_config() at construction time, then independent of it; "
             "re-apply with set_algo_config() if you change the grid model's config afterwards, or "
             "after change_algorithm().")
        .def("set_algo_config", &ContingencyAnalysis::set_algo_config, py::arg("config"),
             "See get_algo_config().")

        // add contingencies
        .def("add_all_n1", &ContingencyAnalysis::add_all_n1, DocContingencyAnalysis::add_all_n1.c_str())
        .def("add_n1", &ContingencyAnalysis::add_n1, DocContingencyAnalysis::add_n1.c_str())
        .def("add_nk", &ContingencyAnalysis::add_nk, DocContingencyAnalysis::add_nk.c_str())
        .def("add_multiple_n1", &ContingencyAnalysis::add_multiple_n1, DocContingencyAnalysis::add_multiple_n1.c_str())

        // remove contingencies
        .def("reset", &ContingencyAnalysis::clear, DocContingencyAnalysis::clear.c_str())
        .def("clear", &ContingencyAnalysis::clear, DocContingencyAnalysis::clear.c_str())
        .def("clear_results_only", &ContingencyAnalysis::clear_results_only, DocContingencyAnalysis::clear.c_str())
        .def("close", &ContingencyAnalysis::clear, DocTimeSeries::clear.c_str())
        .def("remove_n1", &ContingencyAnalysis::remove_n1, DocContingencyAnalysis::remove_n1.c_str())
        .def("remove_nk", &ContingencyAnalysis::remove_nk, DocContingencyAnalysis::remove_nk.c_str())
        .def("remove_multiple_n1", &ContingencyAnalysis::remove_multiple_n1, DocContingencyAnalysis::remove_multiple_n1.c_str())

        // inspect
        .def("my_defaults", &ContingencyAnalysis::my_defaults_vect, DocContingencyAnalysis::my_defaults_vect.c_str())
        .def("is_grid_connected_after_contingency", &ContingencyAnalysis::is_grid_connected_after_contingency, DocLSGrid::_internal_do_not_use.c_str())
        .def("pick_reference_slack", &ContingencyAnalysis::pick_reference_slack,
             "Over the registered contingencies, return the slack bus (gridmodel id) "
             "stranded by the fewest of them — feed it to LSGrid.set_reference_slack_bus "
             "before ac_pf so handle_disconnected_grid skips as few contingencies as possible.")

        // perform computation
        .def("compute", &ContingencyAnalysis::compute, py::call_guard<py::gil_scoped_release>(), DocContingencyAnalysis::compute.c_str())
        .def("compute_flows", &ContingencyAnalysis::compute_flows, DocContingencyAnalysis::compute_flows.c_str())
        .def("compute_power_flows", &ContingencyAnalysis::compute_power_flows, DocContingencyAnalysis::compute_power_flows.c_str())

        // results
        .def("get_flows", &ContingencyAnalysis::get_flows, DocContingencyAnalysis::get_flows.c_str(), py::return_value_policy::reference_internal)
        .def("get_voltages", &ContingencyAnalysis::get_voltages, DocContingencyAnalysis::get_voltages.c_str(), py::return_value_policy::reference_internal)
        .def("get_power_flows", &ContingencyAnalysis::get_power_flows, DocContingencyAnalysis::get_power_flows.c_str(), py::return_value_policy::reference_internal)

        // limit violations (only usable if `compute_limit_violations=True`, see above ;
        // raises otherwise). Row order matches `my_defaults()`.
        .def("converged", [](const ContingencyAnalysis & self){
                 const std::vector<char> & c = self.converged();  // internal storage: char, not
                 // bool, so multi-threaded writes to disjoint indices during compute() can never
                 // race (std::vector<bool> is bit-packed and NOT safe for that). Convert to a
                 // fresh std::vector<bool> here (a copy, so no thread-safety concern) purely so
                 // Python gets a clean list[bool] instead of pybind11's char -> 1-char-str cast.
                 return std::vector<bool>(c.begin(), c.end());
             },
             "Per contingency (row order matches my_defaults()): whether it converged / was "
             "actually simulated (False for skipped or diverged contingencies).")
        .def("get_violations", &ContingencyAnalysis::get_violations,
             "Per contingency (row order matches my_defaults()): list of LimitViolation. A "
             "non-converged contingency (converged() is False) has exactly one LimitViolation "
             "here, with element_type ViolationElementType.GRID and violation_type either "
             "LimitViolationType.NOT_SIMULATED (a pre-check skipped it, eg it splits the grid) or "
             "LimitViolationType.DIVERGENCE (the solver ran but did not converge).",
             py::return_value_policy::reference_internal)
        .def("converged_n", &ContingencyAnalysis::converged_n,
             "Whether the pre-contingency ('n') powerflow converged.")
        .def("get_violations_n", &ContingencyAnalysis::get_violations_n,
             "List of LimitViolation for the pre-contingency ('n') case.",
             py::return_value_policy::reference_internal)

        // timers
        .def("total_time", &ContingencyAnalysis::total_time, DocTimeSeries::total_time.c_str())
        .def("solver_time", &ContingencyAnalysis::solver_time, DocTimeSeries::solver_time.c_str())
        .def("preprocessing_time", &ContingencyAnalysis::preprocessing_time, DocContingencyAnalysis::preprocessing_time.c_str())
        .def("amps_computation_time", &ContingencyAnalysis::amps_computation_time, DocTimeSeries::amps_computation_time.c_str())
        .def("modif_Ybus_time", &ContingencyAnalysis::modif_Ybus_time, DocContingencyAnalysis::modif_Ybus_time.c_str())
        .def("thread_init_time", &ContingencyAnalysis::thread_init_time, "TODO")
        .def("solve_time", &ContingencyAnalysis::solve_time, "TODO")
        .def("nb_solved", &ContingencyAnalysis::nb_solved, DocTimeSeries::nb_solved.c_str());
}
