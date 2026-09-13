// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "batch_algorithm/ContinuationSweep.hpp"
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
        .def("get_linear_solver_stats", &T::get_linear_solver_stats,
             "Linear-solver counters of the whole compute() (nb_analyze, nb_factorize, "
             "nb_refactorize, ...), summed over the member algorithm and every worker "
             "thread's own.\n\n"
             "A batch keeps the Jacobian's sparsity pattern fixed for the whole run, so "
             "each algorithm analyzes ONCE and every row after its first only "
             "refactorizes. That is 1 analyze per algorithm used: 1 single-threaded, and "
             "with nb_thread > 1 one per worker plus the member one that solved the 'n' "
             "warm-up case. More than that means a row changed the sparsity pattern.\n\n"
             "Use get_linear_solver_stats_per_algo() to tell 'every algorithm analyzed "
             "once' from 'one algorithm analyzed several times' -- the sum alone cannot.")
        .def("get_linear_solver_stats_per_algo", &T::get_linear_solver_stats_per_algo,
             "The same counters kept apart, one entry per algorithm: index 0 is the "
             "member one (which solves the 'n' warm-up case, and the whole batch when "
             "single-threaded), then one per worker thread of the last multi-threaded "
             "compute().")
        .def("nb_converged", &T::nb_converged, DocTimeSeries::nb_converged.c_str())
        .def("converged_mask", [](const T & self){
                 const std::vector<char> & c = self.converged_mask();  // char, not bool: see
                 // ContingencyAnalysis's own converged() binding below for why (disjoint
                 // multi-threaded writes into std::vector<bool> are not safe, it is bit-packed).
                 return std::vector<bool>(c.begin(), c.end());
             },
             DocTimeSeries::converged_mask.c_str())

        // base-case reuse (see BaseBatchSweep::set_reuse_base_case)
        .def_property("reuse_base_case",
                      [](const T & self){ return self.get_reuse_base_case(); },
                      [](T & self, bool val){ self.set_reuse_base_case(val); },
                      "Whether the base case is kept between two compute() calls on this object "
                      "(default: ``True``).\n\n"
                      "Before any row is solved, a batch has a base case to establish: read the "
                      "grid (admittance matrix, bus labelling, pv/pq split), walk the graph to "
                      "settle what each contingency strands, solve one 'n' powerflow, and analyze "
                      "and factorize the Jacobian. None of that depends on the injections, so "
                      "calling compute() again -- with new injections, which is what a loop over "
                      "scenarios does -- repeats all of it for nothing. With this ``True`` it is "
                      "done once and kept.\n\n"
                      "It is dropped, and rebuilt on the next compute(), whenever something it is "
                      "made of changes: clear(), change_algorithm(), algo_config, any contingency "
                      "registration, handle_disconnected_grid, nb_thread, init_from_n_powerflow, "
                      "or a different number of simulations. Every such modifier says so itself -- "
                      "internally each names one of three nested cache levels (the grid, this "
                      "batch's inputs, the results) and dropping one drops the levels below it. "
                      "The grid cannot change underneath -- this object holds its own copy of it, "
                      "taken when it was built, and offers no way to modify it.\n\n"
                      "Set it to ``False`` to make every compute() rebuild everything, as it did "
                      "before this existed: useful to tell a suspected caching problem from a real "
                      "one.")
        .def("invalidate_base_case", &T::invalidate_base_case,
             "Drop the kept base case, so the next compute() builds a fresh one (see "
             "``reuse_base_case``). This object already does this itself whenever it changes "
             "anything the base case is made of, so there is normally no reason to call it.")
        .def("base_case_was_reused", &T::base_case_was_reused,
             "Whether the last compute() kept a base case instead of building one. Mostly "
             "of interest when measuring where a batch's time goes.")
        .def_property("reuse_threads",
                      [](const T & self){ return self.get_reuse_threads(); },
                      [](T & self, bool val){ self.set_reuse_threads(val); },
                      "Whether the worker THREADS of the multi-threaded path are kept alive "
                      "between compute() calls (default: ``True``). Independent of "
                      "``reuse_thread_algos``, which keeps what those threads work with: this "
                      "is only about not paying std::thread construction and destruction on "
                      "every call.")
        .def_property("reuse_thread_algos",
                      [](const T & self){ return self.get_reuse_thread_algos(); },
                      [](T & self, bool val){ self.set_reuse_thread_algos(val); },
                      "Whether the worker algorithms of the multi-threaded path are kept "
                      "between compute() calls (default: ``True``), on top of the base case "
                      "itself.\n\n"
                      "At ``nb_thread == 1`` the rows are run by the member algorithm, so "
                      "``reuse_base_case`` alone already keeps its factorization. At "
                      "``nb_thread > 1`` the rows are run by workers, which used to be "
                      "rebuilt every call -- each paying a fresh Jacobian analysis on its "
                      "first row. This gives a threaded batch the same treatment.\n\n"
                      "Has no effect when ``reuse_base_case`` is ``False``: the workers are "
                      "part of the same batch state.")
        .def("thread_algos_were_reused", &T::thread_algos_were_reused,
             "Whether the last compute() also kept the WORKER algorithms of the "
             "multi-threaded path, instead of building and analyzing one per thread. "
             "Always False for a single-threaded batch: it has no workers (the member "
             "algorithm runs the rows itself, and keeping that is base_case_was_reused).")

        // reverse-mode differentiation (see BatchAdjoint.hpp)
        .def_property("keep_jacobian",
                      [](const T & self){ return self.get_keep_jacobian(); },
                      [](T & self, bool val){ self.set_keep_jacobian(val); },
                      "Whether each row's converged Jacobian is kept during compute(), so that "
                      "solve_JT() can run afterwards. Defaults to ``False``; must be set BEFORE "
                      "compute().\n\n"
                      "This is what makes a batch differentiable: with the Jacobians kept, one "
                      "transposed solve per row turns a loss's sensitivity to the voltages into "
                      "its gradient with respect to every injection of every row (the adjoint "
                      "method / the implicit function theorem).\n\n"
                      "It costs `nb_rows * nnz(J)` floats -- see adjoint_memory_bytes(), and mind "
                      "that on a large grid with many rows this is gigabytes -- plus one extra "
                      "Jacobian evaluation per row (the Newton-Raphson loop stops on the Jacobian "
                      "of its previous iterate, which is not quite the one at the solution). Only "
                      "the AC Newton-Raphson algorithms build a Jacobian at all; compute() raises "
                      "with any other.")
        .def("adjoint_memory_bytes", &T::adjoint_memory_bytes,
             "Bytes the kept Jacobians occupy after a compute() (0 when none were kept). "
             "Grows as nb_rows * nnz(J): worth reading before scaling a batch up.")
        .def("dim_J", &T::dim_J,
             "Dimension of the augmented Newton-Raphson system -- the length of one "
             "cotangent handed to solve_JT, and of one row of its result. 0 until a "
             "compute() that kept the Jacobians.")
        .def("get_theta_col_of_bus", &T::get_theta_col_of_bus,
             "For each grid bus (same numbering as the columns of get_voltages()), the column "
             "of the Jacobian holding its voltage-ANGLE unknown, or -1 where it has none (a "
             "reference slack, a bus outside the solver). Where the angle part of that bus's "
             "cotangent goes in solve_JT's input.")
        .def("get_vm_col_of_bus", &T::get_vm_col_of_bus,
             "Same as get_theta_col_of_bus for the voltage-MAGNITUDE unknown: -1 at a bus whose "
             "magnitude is held fixed (a PV bus, unless a generator contingency may release it).")
        .def("get_p_row_of_bus", &T::get_p_row_of_bus,
             "For each grid bus, the row of the Jacobian holding its ACTIVE power mismatch "
             "equation, or -1 where it has none. Where that bus's active injection gradient is "
             "read out of solve_JT's result.")
        .def("get_q_row_of_bus", &T::get_q_row_of_bus,
             "Same as get_p_row_of_bus for the REACTIVE power mismatch equation.")
        .def("solve_JT", &T::solve_JT, py::arg("xbar"),
             "Solve the adjoint system `J_i^T . lambda_i = xbar_i` of every row, and return "
             "lambda. Requires `keep_jacobian` to have been True during compute().\n\n"
             "`xbar` has shape (nb_rows, k * dim_J): one row per simulation, holding k "
             "cotangents of dim_J coefficients laid end to end. Differentiating a scalar loss "
             "uses k = 1; several directions are only needed for a full Jacobian, and cost one "
             "triangular solve each rather than one factorization each.\n\n"
             "lambda is the gradient with respect to the per-unit bus injection: its real part "
             "sits at get_p_row_of_bus()[bus], its imaginary part at get_q_row_of_bus()[bus]. "
             "Rows that did not converge come back as zeros -- see adjoint_row_ok().")
        .def("adjoint_row_ok", [](const T & self){
                 const std::vector<char> & c = self.adjoint_row_ok();
                 return std::vector<bool>(c.begin(), c.end());
             },
             "Per row: True where the last solve_JT() actually solved that row's adjoint "
             "system. False for a row the batch never solved (it diverged, or a contingency "
             "islanded it), whose lambda is zero.")
        .def("adjoint_solver_stats", &T::adjoint_solver_stats,
             "Linear-solver counters and timings of the last solve_JT(), summed over its "
             "worker threads: how much of the backward pass went into refactorizing each "
             "row's Jacobian versus into the transposed solves themselves.")

        // status
        // `<>`: see the comment above modify_gen_p -- Clang (used for the macOS/arm64
        // wheels) rejects `&T::get_status` outright for a member function template
        // with no target type to deduce against, where GCC happens to accept it via
        // the default template argument. Always spell it out explicitly.
        .def("get_status", &T::template get_status<>, DocTimeSeries::get_status.c_str())
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
        .def("modify_gen_v", &T::template modify_gen_v<>, py::arg("gen_v"),
             "Per-step generator target voltage magnitude, shape (n_simul, n_gen), in pu "
             "(vm_pu), NOT kV. "
             "Unlike modify_gen_p/modify_sgen_p/modify_load_p/modify_load_q, this does NOT "
             "feed the injection (Sbus) -- it only re-seeds |V| at each voltage-regulating "
             "generator's regulated bus before that step's solve. See modify_gen_p() for "
             "the shared row-count-lock behavior.")
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

    // ContinuationSweepCPP: a SIBLING of the four BaseBatchSweep instantiations, not
    // a fifth one -- see batch_algorithm/ContinuationSweep.hpp for why. It shares the
    // base class (and therefore the one-analyze-per-run guarantee, the timers and the
    // result accessors) but decides its own rows, so it has its own binding block
    // rather than going through bind_batch_sweep_common.
    py::class_<ContinuationSweep>(m, "ContinuationSweepCPP",
        "Continuation powerflow: traces the solution curve from the grid's own "
        "injection state (lambda = 0) to a target one (lambda = 1), and stops at the "
        "voltage-collapse 'nose' or at a requested lambda.\n\n"
        "Set the target with set_target_load_p / set_target_load_q / set_target_gen_p "
        "/ set_target_sgen_p (any axis left unset stays at the grid's own values), "
        "then call compute(). The option names follow MATPOWER's cpf.* so a runcpf "
        "user is on familiar ground. Newton-Raphson algorithms only.\n\n"
        "Prefer the `lightsim2grid.continuationPowerflow.ContinuationPowerFlow` "
        "wrapper, which builds the target from a loading factor and per-load / "
        "per-generator steering vectors.")
        .def(py::init<const LSGrid &>())

        // solver control
        .def("change_algorithm", py::overload_cast<const AlgorithmType&>(&ContinuationSweep::change_algorithm), DocLSGrid::change_algorithm.c_str())
        .def("change_algorithm", py::overload_cast<const std::string&>(&ContinuationSweep::change_algorithm), DocLSGrid::change_algorithm_by_name.c_str())
        .def("available_default_algorithms", &ContinuationSweep::available_default_algorithms, DocLSGrid::available_default_algorithms.c_str())
        .def("get_algo_type", &ContinuationSweep::get_algo_type, DocLSGrid::get_algo_type.c_str())
        .def("get_algo_name", &ContinuationSweep::get_algo_name, "Registry name of the selected algorithm.")
        .def("get_algo_config", &ContinuationSweep::get_algo_config, "Config of the internal solver.")
        .def("set_algo_config", &ContinuationSweep::set_algo_config, py::arg("config"), "See get_algo_config().")

        // the target state
        .def("set_target_gen_p", &ContinuationSweep::set_target_gen_p, py::arg("gen_p"),
             "Target active generator setpoints (n_gen,), in MW. Unset means 'unchanged'.")
        .def("set_target_sgen_p", &ContinuationSweep::set_target_sgen_p, py::arg("sgen_p"),
             "Target active static-generator setpoints (n_sgen,), in MW. Unset means 'unchanged'.")
        .def("set_target_load_p", &ContinuationSweep::set_target_load_p, py::arg("load_p"),
             "Target active load setpoints (n_load,), in MW. Unset means 'unchanged'.")
        .def("set_target_load_q", &ContinuationSweep::set_target_load_q, py::arg("load_q"),
             "Target reactive load setpoints (n_load,), in MVAr. Unset means 'unchanged'.")
        .def("clear_target", &ContinuationSweep::clear_target, "Forget every target axis set so far.")

        // options (MATPOWER cpf.* names and defaults)
        .def_property("step", &ContinuationSweep::get_step, &ContinuationSweep::set_step,
                      "Nominal continuation step, as an arc length along the unit tangent "
                      "(MATPOWER cpf.step, default 0.05).")
        .def_property("step_min", &ContinuationSweep::get_step_min, &ContinuationSweep::set_step_min,
                      "Smallest step the corrector-failure retry may shrink to; failing at this "
                      "step ends the curve (MATPOWER cpf.step_min, default 1e-4).")
        .def_property("step_max", &ContinuationSweep::get_step_max, &ContinuationSweep::set_step_max,
                      "Largest step the adaptation may grow to (MATPOWER cpf.step_max, default 0.2).")
        .def_property("adapt_step", &ContinuationSweep::get_adapt_step, &ContinuationSweep::set_adapt_step,
                      "Adapt the step to the predictor's error (MATPOWER cpf.adapt_step, default False).")
        .def_property("adapt_step_damping", &ContinuationSweep::get_adapt_step_damping, &ContinuationSweep::set_adapt_step_damping,
                      "Damping of the step adaptation (MATPOWER cpf.adapt_step_damping, default 0.7).")
        .def_property("adapt_step_tol", &ContinuationSweep::get_adapt_step_tol, &ContinuationSweep::set_adapt_step_tol,
                      "Target predictor error the adaptation aims at (MATPOWER cpf.adapt_step_tol, default 1e-3).")
        .def_property("nose_tol", &ContinuationSweep::get_nose_tol, &ContinuationSweep::set_nose_tol,
                      "The curve is declared at the nose when the tangent's lambda component falls "
                      "below this (MATPOWER cpf.nose_tol, default 1e-5). Note that this "
                      "parameterisation makes that component strictly positive, tending to zero at "
                      "the nose -- it is a threshold, never a sign change.")
        .def_property("stop_at_lam", &ContinuationSweep::get_stop_at_lam, &ContinuationSweep::set_stop_at_lam,
                      "Stop once lambda reaches this value (MATPOWER's numeric cpf.stop_at). "
                      "Non-positive (the default) means 'trace until the nose'.")
        .def_property("max_steps", &ContinuationSweep::get_max_steps, &ContinuationSweep::set_max_steps,
                      "Hard cap on the number of traced points (default 1000).")
        .def_property("exact_tangent", &ContinuationSweep::get_exact_tangent, &ContinuationSweep::set_exact_tangent,
                      "Rebuild and refactorize the Jacobian at each converged point before taking "
                      "its tangent (default False). Off, the tangent uses the factorization the "
                      "corrector left standing, which is one NR iterate behind the converged point.")

        .def("compute", &ContinuationSweep::compute, py::call_guard<py::gil_scoped_release>(),
             py::arg("Vinit"), py::arg("max_iter"), py::arg("tol"),
             "Trace the curve. Raises if the target is identical to the base state, or if "
             "the selected algorithm is not Newton-Raphson based.")

        // results
        .def("get_status", &ContinuationSweep::get_status,
             "1 if the curve reached its requested end (the nose, or stop_at_lam), 0 otherwise.")
        .def("get_msg", &ContinuationSweep::get_msg, "Why the run stopped, in words.")
        .def("nb_points", &ContinuationSweep::nb_points,
             "Number of traced points, the base case included.")
        .def("get_lam", &ContinuationSweep::get_lam, py::return_value_policy::reference_internal,
             "Lambda at each traced point; lam[0] == 0 is the base case, lam == 1 the target.")
        .def("get_tangent_lam", &ContinuationSweep::get_tangent_lam, py::return_value_policy::reference_internal,
             "The tangent's lambda component at each point, in (0, 1]; it tends to 0 at the "
             "nose. The last point has none and reads 0.")
        .def("get_lam_max", &ContinuationSweep::get_lam_max, "Largest lambda reached.")
        .def("nb_retries", &ContinuationSweep::nb_retries,
             "How many times a corrector failed and the step had to be halved.")
        .def("get_direction_solver", &ContinuationSweep::get_direction_solver, py::return_value_policy::reference_internal,
             "The direction actually used (Sbus_target - Sbus_base), in solver bus ordering "
             "and per unit.")
        .def("get_voltages", &ContinuationSweep::get_voltages, DocTimeSeries::get_voltages.c_str(), py::return_value_policy::reference_internal)
        .def("compute_flows", &ContinuationSweep::compute_flows, DocTimeSeries::compute_flows.c_str())
        .def("compute_power_flows", &ContinuationSweep::compute_power_flows, DocTimeSeries::compute_power_flows.c_str())
        .def("get_flows", &ContinuationSweep::get_flows, DocTimeSeries::get_flows.c_str(), py::return_value_policy::reference_internal)
        .def("get_power_flows", &ContinuationSweep::get_power_flows, DocTimeSeries::get_power_flows.c_str(), py::return_value_policy::reference_internal)

        // timers / counters
        .def("total_time", &ContinuationSweep::total_time, DocTimeSeries::total_time.c_str())
        .def("solver_time", &ContinuationSweep::solver_time, DocTimeSeries::solver_time.c_str())
        .def("preprocessing_time", &ContinuationSweep::preprocessing_time, DocTimeSeries::preprocessing_time.c_str())
        .def("nb_solved", &ContinuationSweep::nb_solved, DocTimeSeries::nb_solved.c_str())
        .def("nb_converged", &ContinuationSweep::nb_converged, DocTimeSeries::nb_converged.c_str())
        .def("get_linear_solver_stats", &ContinuationSweep::get_linear_solver_stats,
             "Linear-solver counters for the whole curve. nb_analyze must be 1 however many "
             "points were traced -- that is the entire reason a continuation belongs in the "
             "batch layer. More than 1 means something changed the Jacobian's sparsity.")
        .def("clear", &ContinuationSweep::clear, DocTimeSeries::clear.c_str())
        .def("close", &ContinuationSweep::clear, DocTimeSeries::clear.c_str())

        // bound although any value but 1 is rejected, for the same reason TimeSeriesCPP
        // binds it: a user who finds the attribute gets an explanation instead of an
        // AttributeError. The points of a curve are chained, so there is nothing to split.
        .def("set_nb_thread", &ContinuationSweep::set_nb_thread, py::arg("nb_thread"),
             "Always 1: the points of a continuation are chained (each is predicted from "
             "the previous one's tangent), so the curve cannot be split over threads.")
        .def_property("nb_thread",
                      [](const ContinuationSweep & self){ return self.get_nb_thread(); },
                      [](ContinuationSweep & self, int val){ self.set_nb_thread(val); },
                      DocTimeSeries::nb_thread.c_str());

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
        .def("set_contingency_gens", &ScenarioSweep::set_contingency_gens<>, py::arg("mask"),
             "Per-step generator contingency mask, shape (n_simul, n_gen), dtype bool. "
             "True means 'disconnect this generator for this simulation'.\n\n"
             "Unlike the two branch masks this does not edit Ybus -- a generator has no "
             "admittance. It removes the generator's active power (and, if it does not "
             "regulate voltage, its reactive setpoint) from that step's injection, "
             "re-weights the distributed slack without it, and -- when the LAST "
             "generator regulating its own bus is taken out -- turns that bus from PV "
             "to PQ for the step, so its voltage magnitude is solved for instead of "
             "held at the setpoint. The lost MW is picked up by the slack; use "
             "modify_gen_p to express a redispatch instead.\n\n"
             "The PV/PQ relabelling costs no extra symbolic factorization: every bus "
             "that can flip is given a voltage-magnitude unknown and a reactive "
             "equation once, up front, and each step merely masks the equation of the "
             "buses that are still PV. The whole sweep keeps running on one analysis.\n\n"
             "It is NOT free, though: those reserved unknowns and equations make the "
             "Jacobian bigger for EVERY step, whether or not that step disconnects "
             "anything -- one extra row and column per bus that can flip. A handful of "
             "candidate generators is negligible; masking every generator on the grid "
             "grows the Jacobian's dimension by roughly the number of PV buses, and "
             "every factorization and solve pays for it.\n\n"
             "Only generators regulating their OWN bus are supported. compute() raises "
             "if the mask names a generator that regulates a remote bus, or one whose "
             "bus a control group holds (a remote generator, an SVC or an HVDC "
             "converter station): remote voltage control is not supported yet.")

        // limit violations + "handle disconnected grid": same names/semantics as
        // ContingencyAnalysisCPP (see below), now also available here. Deliberately
        // NO converged()/converged_n() -- a non-converged row's get_violations()
        // entry already carries a GRID-type NOT_SIMULATED/DIVERGENCE sentinel
        // LimitViolation, so a separate convergence flag would be redundant; a
        // diverging pre-batch "n" powerflow is likewise stamped with that same
        // sentinel rather than left looking like "converged, no violations".
        .def_property("compute_limit_violations",
                      [](const ScenarioSweep & self){ return self.get_compute_limit_violations(); },
                      [](ScenarioSweep & self, bool val){ self.set_compute_limit_violations(val); },
                      "Whether limit violations are computed inline, per row, during compute() "
                      "(see get_violations() / get_violations_n()). Defaults to ``False``. "
                      "Computing violations means an extra per-element current / voltage check "
                      "in every row's solve, so users who only need compute_flows() / get_flows() "
                      "should leave this off. Changing this flag clears any previously-computed "
                      "results. Unlike ContingencyAnalysisCPP, there is no converged() / "
                      "converged_n() here: a non-converged row's get_violations() entry already "
                      "carries a GRID-type NOT_SIMULATED / DIVERGENCE LimitViolation, which fully "
                      "encodes that row's status by itself.")
        .def_property("violation_threshold",
                      [](const ScenarioSweep & self){ return self.get_violation_threshold(); },
                      [](ScenarioSweep & self, real_type val){ self.set_violation_threshold(val); },
                      DocContingencyAnalysis::violation_threshold.c_str())
        .def_property("handle_disconnected_grid",
                      [](const ScenarioSweep & self){ return self.get_handle_disconnected_grid(); },
                      [](ScenarioSweep & self, bool val){ self.set_handle_disconnected_grid(val); },
                      "Whether to simulate a row whose contingency splits the grid into multiple "
                      "connected components. When False (default) such a row is skipped (its "
                      "voltages are left at 0), reproducing the legacy behaviour. When True, the "
                      "largest connected component is solved while the buses of the other "
                      "component(s) are masked (their voltage is reported as 0). Supported by the "
                      "Newton-Raphson family (AC) and the DC solver; a non Newton-Raphson AC "
                      "algorithm is rejected.")
        .def("get_violations", &ScenarioSweep::get_violations<>,
             "Per row (same order as every modify_* / set_contingency_* input): list of "
             "LimitViolation. A row that did not converge has exactly one LimitViolation here, "
             "with element_type ViolationElementType.GRID and violation_type either "
             "LimitViolationType.NOT_SIMULATED (a pre-check skipped it, eg it splits the grid "
             "with handle_disconnected_grid off) or LimitViolationType.DIVERGENCE (the solver "
             "ran but did not converge, including a diverging pre-batch \"n\" powerflow, which "
             "stamps every row this way). Requires compute_limit_violations=True.",
             py::return_value_policy::reference_internal)
        .def("get_violations_n", &ScenarioSweep::get_violations_n<>,
             "List of LimitViolation for the pre-batch (\"n\") case (no injection change, no "
             "contingency) shared by every row. Requires compute_limit_violations=True.",
             py::return_value_policy::reference_internal)
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
        // `<>` on every one of these below: they are member function TEMPLATES
        // (SFINAE-gated on YbusPolicy::supports_contingency && !SbusPolicy::supports_vary,
        // see BaseBatchSweep.hpp), and `&ContingencyAnalysis::method` (no explicit
        // template-argument list) is only accepted by GCC, which happens to fall back to
        // the default template argument -- Clang and MSVC both reject it outright ("no
        // matching member function for call to 'def'" / "no matching overloaded function
        // found"), since there is no target type here to deduce the SFINAE parameter
        // against. `<>` forces the default template arguments explicitly, which every
        // compiler accepts.
        .def("add_all_n1", &ContingencyAnalysis::add_all_n1<>, DocContingencyAnalysis::add_all_n1.c_str())
        .def("add_n1", &ContingencyAnalysis::add_n1<>, DocContingencyAnalysis::add_n1.c_str())
        .def("add_nk", &ContingencyAnalysis::add_nk<>, DocContingencyAnalysis::add_nk.c_str())
        .def("add_multiple_n1", &ContingencyAnalysis::add_multiple_n1<>, DocContingencyAnalysis::add_multiple_n1.c_str())

        // remove contingencies
        .def("reset", &ContingencyAnalysis::clear, DocContingencyAnalysis::clear.c_str())
        .def("clear", &ContingencyAnalysis::clear, DocContingencyAnalysis::clear.c_str())
        .def("clear_results_only", &ContingencyAnalysis::clear_results_only<>, DocContingencyAnalysis::clear.c_str())
        .def("close", &ContingencyAnalysis::clear, DocTimeSeries::clear.c_str())
        .def("remove_n1", &ContingencyAnalysis::remove_n1<>, DocContingencyAnalysis::remove_n1.c_str())
        .def("remove_nk", &ContingencyAnalysis::remove_nk<>, DocContingencyAnalysis::remove_nk.c_str())
        .def("remove_multiple_n1", &ContingencyAnalysis::remove_multiple_n1<>, DocContingencyAnalysis::remove_multiple_n1.c_str())

        // inspect
        .def("my_defaults", &ContingencyAnalysis::my_defaults_vect<>, DocContingencyAnalysis::my_defaults_vect.c_str())
        .def("is_grid_connected_after_contingency", &ContingencyAnalysis::is_grid_connected_after_contingency<>, DocLSGrid::_internal_do_not_use.c_str())
        .def("pick_reference_slack", &ContingencyAnalysis::pick_reference_slack<>,
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
        .def("get_violations", &ContingencyAnalysis::get_violations<>,
             "Per contingency (row order matches my_defaults()): list of LimitViolation. A "
             "non-converged contingency (converged() is False) has exactly one LimitViolation "
             "here, with element_type ViolationElementType.GRID and violation_type either "
             "LimitViolationType.NOT_SIMULATED (a pre-check skipped it, eg it splits the grid) or "
             "LimitViolationType.DIVERGENCE (the solver ran but did not converge).",
             py::return_value_policy::reference_internal)
        .def("converged_n", &ContingencyAnalysis::converged_n<>,
             "Whether the pre-contingency ('n') powerflow converged.")
        .def("get_violations_n", &ContingencyAnalysis::get_violations_n<>,
             "List of LimitViolation for the pre-contingency ('n') case.",
             py::return_value_policy::reference_internal)

        // timers
        .def("total_time", &ContingencyAnalysis::total_time, DocTimeSeries::total_time.c_str())
        .def("solver_time", &ContingencyAnalysis::solver_time, DocTimeSeries::solver_time.c_str())
        .def("preprocessing_time", &ContingencyAnalysis::preprocessing_time, DocContingencyAnalysis::preprocessing_time.c_str())
        .def("amps_computation_time", &ContingencyAnalysis::amps_computation_time, DocTimeSeries::amps_computation_time.c_str())
        .def("modif_Ybus_time", &ContingencyAnalysis::modif_Ybus_time<>, DocContingencyAnalysis::modif_Ybus_time.c_str())
        .def("thread_init_time", &ContingencyAnalysis::thread_init_time, DocTimeSeries::thread_init_time.c_str())
        .def("solve_time", &ContingencyAnalysis::solve_time<>, "TODO")
        .def("nb_solved", &ContingencyAnalysis::nb_solved, DocTimeSeries::nb_solved.c_str())
        .def("nb_converged", &ContingencyAnalysis::nb_converged, DocTimeSeries::nb_converged.c_str())
        .def("converged_mask", [](const ContingencyAnalysis & self){
                 const std::vector<char> & c = self.converged_mask();
                 return std::vector<bool>(c.begin(), c.end());
             },
             DocTimeSeries::converged_mask.c_str());
}
