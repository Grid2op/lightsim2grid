// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "batch_algorithm/LimitViolation.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

namespace {

/**
 * The interface shared by TimeSeriesCPP, InjectionSweepCPP and ScenarioSweepCPP: every
 * instantiation of BaseBatchSweep whose SbusPolicy varies (see
 * batch_algorithm/BaseBatchSweep.hpp). Templated on the instantiation rather than
 * copy-pasted so the three Python classes cannot drift apart.
 *
 * Deliberately does NOT bind compute_Vs()/get_sbuses(): those only exist on
 * TimeSeries/InjectionSweep (see bind_legacy_compute_vs below) -- ScenarioSweep never
 * had the old bundled-call API, so binding them unconditionally here would fail to
 * compile for it (the C++ side gates them out via SFINAE).
 */
template<class T>
void bind_batch_sweep_common(py::class_<T> & cls)
{
    cls
        .def(py::init<const LSGrid &>())

        // solver control
        .def("change_algorithm", py::overload_cast<const AlgorithmType&>(&T::change_algorithm), DocLSGrid::change_algorithm.c_str())
        .def("change_algorithm", py::overload_cast<const std::string&>(&T::change_algorithm), DocLSGrid::change_algorithm_by_name.c_str())
        .def("change_solver", py::overload_cast<const AlgorithmType&>(&T::change_algorithm), "DEPRECATED: use 'change_algorithm' instead")
        .def("change_solver", py::overload_cast<const std::string&>(&T::change_algorithm), "DEPRECATED: use 'change_algorithm' instead")
        .def("available_default_algorithms", &T::available_default_algorithms, DocLSGrid::available_default_algorithms.c_str())
        .def("available_algorithm_names", &T::available_algorithm_names, DocLSGrid::available_algorithm_names.c_str())
        .def("get_algo_type", &T::get_algo_type, DocLSGrid::get_algo_type.c_str())
        .def("get_algo_name", &T::get_algo_name,
             "Registry name of the currently selected algorithm. Unlike get_algo_type(), stays "
             "meaningful for plugin solvers (and built-ins with no dedicated AlgorithmType member).")
        .def("get_algo_config", &T::get_algo_config,
             "Config (eg ScalingPolicyType / damping parameters) of the internal solver used for "
             "every step. Copied once from the grid model's own get_ac_algo_config() at "
             "construction time, then independent of it; re-apply with set_algo_config() if you "
             "change the grid model's config afterwards, or after change_algorithm().")
        .def("set_algo_config", &T::set_algo_config, py::arg("config"),
             "See get_algo_config().")

        // timers
        .def("total_time", &T::total_time, DocTimeSeries::total_time.c_str())
        .def("solver_time", &T::solver_time, DocTimeSeries::solver_time.c_str())
        .def("preprocessing_time", &T::preprocessing_time, DocTimeSeries::preprocessing_time.c_str())
        .def("amps_computation_time", &T::amps_computation_time, DocTimeSeries::amps_computation_time.c_str())
        .def("thread_init_time", &T::thread_init_time, DocTimeSeries::thread_init_time.c_str())
        .def("nb_solved", &T::nb_solved, DocTimeSeries::nb_solved.c_str())

        // status
        .def("get_status", &T::get_status, DocTimeSeries::get_status.c_str())
        .def("clear", &T::clear, DocTimeSeries::clear.c_str())
        .def("close", &T::clear, DocTimeSeries::clear.c_str())

        // new setter-based API (replaces the bundled compute_Vs call): build up the
        // per-step injection with as many of these as are relevant, then call
        // compute(). Any axis never set defaults to the grid's own target value,
        // broadcast across every row. All setters share one row-count lock: the
        // first one called fixes the number of simulations, every later one
        // (including set_contingency_lines/trafos on ScenarioSweep) is checked
        // against it immediately.
        // `&T::template modify_gen_p<>` (explicit empty template-argument list, with
        // the `template` disambiguator since T is itself a dependent name here): a
        // member function TEMPLATE's address cannot be taken via plain `&T::method`
        // in a context with no target type to deduce against (which is exactly what
        // pybind11's `.def(name, F&&, ...)` is -- it deduces `F` FROM this
        // expression, so there is nothing to deduce the SFINAE template parameter
        // against). `<>` forces the compiler to use the default template arguments
        // instead.
        .def("modify_gen_p", &T::template modify_gen_p<>, py::arg("gen_p"),
             "Per-step active generator setpoints, shape (n_simul, n_gen). See the class "
             "docstring: locks / checks the number of simulations against any other "
             "modify_* / set_contingency_* call already made on this object.")
        .def("modify_sgen_p", &T::template modify_sgen_p<>, py::arg("sgen_p"),
             "Per-step active static generator setpoints, shape (n_simul, n_sgen). "
             "See modify_gen_p().")
        .def("modify_load_p", &T::template modify_load_p<>, py::arg("load_p"),
             "Per-step active load setpoints, shape (n_simul, n_load). See modify_gen_p().")
        .def("modify_load_q", &T::template modify_load_q<>, py::arg("load_q"),
             "Per-step reactive load setpoints, shape (n_simul, n_load). See modify_gen_p().")
        .def("compute", &T::compute, py::call_guard<py::gil_scoped_release>(),
             py::arg("Vinit"), py::arg("max_iter"), py::arg("tol"),
             "Run the batch: one powerflow per simulation, using whatever was set by "
             "modify_* (and, on ScenarioSweep, set_contingency_lines / "
             "set_contingency_trafos). Raises if nothing was ever set.")

        // results
        .def("compute_flows", &T::compute_flows, DocTimeSeries::compute_flows.c_str())
        .def("compute_power_flows", &T::compute_power_flows, DocTimeSeries::compute_power_flows.c_str())
        .def("get_flows", &T::get_flows, DocTimeSeries::get_flows.c_str(), py::return_value_policy::reference_internal)
        .def("get_power_flows", &T::get_power_flows, DocTimeSeries::get_power_flows.c_str(), py::return_value_policy::reference_internal)
        .def("get_voltages", &T::get_voltages, DocTimeSeries::get_voltages.c_str(), py::return_value_policy::reference_internal)

        // nb_thread is bound for every one of these classes on purpose, even though
        // TimeSeriesCPP rejects any value but 1 (see the warning in its docstring): a
        // user who discovers the attribute gets an error message pointing at
        // InjectionSweepCPP/ScenarioSweepCPP, instead of an AttributeError that
        // explains nothing.
        .def_property("nb_thread",
                      [](const T & self){ return self.get_nb_thread(); },
                      [](T & self, int val){ self.set_nb_thread(val); },
                      DocTimeSeries::nb_thread.c_str());
}

/**
 * compute_Vs()/get_sbuses(): the legacy bundled-call API, kept (deprecated via
 * docstring only, not removed) on TimeSeriesCPP/InjectionSweepCPP for backwards
 * compatibility. Internally a thin wrapper around the 4 modify_* setters + compute().
 * NOT available on ScenarioSweepCPP -- it never had this call, see
 * bind_batch_sweep_common's own docstring above.
 */
template<class T>
void bind_legacy_compute_vs(py::class_<T> & cls)
{
    cls
        .def("compute_Vs", &T::template compute_Vs<>, py::call_guard<py::gil_scoped_release>(),
             (DocTimeSeries::compute_Vs + " DEPRECATED: prefer modify_gen_p / modify_sgen_p / "
              "modify_load_p / modify_load_q + compute() instead.").c_str())
        .def("get_sbuses", &T::template get_sbuses<>, DocTimeSeries::get_sbuses.c_str(), py::return_value_policy::reference_internal);
}

}  // namespace

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

    // TimeSeriesCPP, InjectionSweepCPP and ScenarioSweepCPP are three instantiations
    // of the same C++ template (see batch_algorithm/BaseBatchSweep.hpp): same
    // per-step injection inputs, same results -- they differ in how each step is
    // initialized (TimeSeries/InjectionSweep) and in whether a contingency also
    // varies per step (ScenarioSweep). Everything but the class docstring and
    // `init_from_n_powerflow` is bound by the shared helpers above.
    py::class_<TimeSeries> time_series(m, "TimeSeriesCPP", DocTimeSeries::TimeSeries.c_str());
    bind_batch_sweep_common(time_series);
    bind_legacy_compute_vs(time_series);
    time_series
        .def_property("init_from_n_powerflow",
                      [](const TimeSeries & self){ return self.get_init_from_n_powerflow(); },
                      [](TimeSeries & self, bool val){ self.set_init_from_n_powerflow(val); },
                      DocTimeSeries::init_from_n_powerflow.c_str());

    py::class_<InjectionSweep> injection_sweep(m, "InjectionSweepCPP", DocInjectionSweep::InjectionSweep.c_str());
    bind_batch_sweep_common(injection_sweep);
    bind_legacy_compute_vs(injection_sweep);
    injection_sweep
        .def_property("init_from_n_powerflow",
                      [](const InjectionSweep & self){ return self.get_init_from_n_powerflow(); },
                      [](InjectionSweep & self, bool val){ self.set_init_from_n_powerflow(val); },
                      DocInjectionSweep::init_from_n_powerflow.c_str());

    // ScenarioSweepCPP: the 4th instantiation -- varies both the injection AND a
    // contingency (line/trafo disconnection) per row, independently, row-aligned.
    py::class_<ScenarioSweep> scenario_sweep(m, "ScenarioSweepCPP",
        "Batch powerflow varying both the injection AND a contingency per simulation, "
        "row-aligned: row i of every modify_* input is solved together with row i of "
        "set_contingency_lines / set_contingency_trafos. Build up the batch with "
        "modify_gen_p / modify_sgen_p / modify_load_p / modify_load_q and "
        "set_contingency_lines / set_contingency_trafos (any axis never set defaults "
        "to the grid's own state for every row), then call compute(). Unlike "
        "ContingencyAnalysisCPP's add_n1/add_nk (a set of distinct scenarios applied "
        "to one shared base case), set_contingency_lines/trafos are dense boolean "
        "masks of shape (n_simul, n_lines) / (n_simul, n_trafos) -- True means "
        "'deactivate this branch for this simulation'; the two APIs are deliberately "
        "not unified, they serve different usages.");
    bind_batch_sweep_common(scenario_sweep);
    scenario_sweep
        .def("set_contingency_lines", &ScenarioSweep::set_contingency_lines<>, py::arg("mask"),
             "Per-step powerline contingency mask, shape (n_simul, n_line), dtype bool. "
             "True means 'deactivate this powerline for this simulation'. See the class "
             "docstring: locks / checks the number of simulations, and is a different "
             "API from ContingencyAnalysisCPP's add_n1/add_nk on purpose.")
        .def("set_contingency_trafos", &ScenarioSweep::set_contingency_trafos<>, py::arg("mask"),
             "Per-step trafo contingency mask, shape (n_simul, n_trafo), dtype bool. "
             "See set_contingency_lines().")
        .def_property("init_from_n_powerflow",
                      [](const ScenarioSweep & self){ return self.get_init_from_n_powerflow(); },
                      [](ScenarioSweep & self, bool val){ self.set_init_from_n_powerflow(val); },
                      "Whether to initialize the complex voltages of each simulation with the "
                      "results of a n-powerflow (ie a powerflow with no injection change and no "
                      "contingency) or not. Default: false, meaning each simulation is "
                      "initialized with the given input vector.");

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
        .def_property("violation_threshold",
                      [](const ContingencyAnalysis & self){ return self.get_violation_threshold(); },
                      [](ContingencyAnalysis & self, real_type val){ self.set_violation_threshold(val); },
                      DocContingencyAnalysis::violation_threshold.c_str())
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
        .def("change_algorithm", py::overload_cast<const AlgorithmType&>(&ContingencyAnalysis::change_algorithm), DocLSGrid::change_algorithm.c_str())
        .def("change_algorithm", py::overload_cast<const std::string&>(&ContingencyAnalysis::change_algorithm), DocLSGrid::change_algorithm_by_name.c_str())
        .def("change_solver", py::overload_cast<const AlgorithmType&>(&ContingencyAnalysis::change_algorithm), "DEPRECATED: use 'change_algorithm' instead")
        .def("change_solver", py::overload_cast<const std::string&>(&ContingencyAnalysis::change_algorithm), "DEPRECATED: use 'change_algorithm' instead")
        .def("available_default_algorithms", &ContingencyAnalysis::available_default_algorithms, DocLSGrid::available_default_algorithms.c_str())
        .def("available_algorithm_names", &ContingencyAnalysis::available_algorithm_names, DocLSGrid::available_algorithm_names.c_str())
        .def("get_algo_type", &ContingencyAnalysis::get_algo_type, DocLSGrid::get_algo_type.c_str())
        .def("get_algo_name", &ContingencyAnalysis::get_algo_name,
             "Registry name of the currently selected algorithm. Unlike get_algo_type(), stays "
             "meaningful for plugin solvers (and built-ins with no dedicated AlgorithmType member).")
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
        .def("thread_init_time", &ContingencyAnalysis::thread_init_time, DocTimeSeries::thread_init_time.c_str())
        .def("solve_time", &ContingencyAnalysis::solve_time, "TODO")
        .def("nb_solved", &ContingencyAnalysis::nb_solved, DocTimeSeries::nb_solved.c_str());
}
