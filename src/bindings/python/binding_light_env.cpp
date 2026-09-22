// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "LSGrid.hpp"
#include "light_env/inj_action.hpp"
#include "light_env/topo_action.hpp"
#include "light_env/protections.hpp"
#include "light_env/light_env.hpp"

using namespace ls2g;

void bind_light_env(py::module_& m) {
    // enum for the topology element
    py::enum_<ElementType>(m, "ElementType", "This enum controls the topology action in the light environment of lightsim2grid")
        .value("load", ElementType::load, "denotes a load")
        .value("gen", ElementType::gen, "denotes a gen")
        .value("line_or", ElementType::line_or, "denotes a line (origin side)")
        .value("line_ex", ElementType::line_ex, "denotes a line (ext side)")
        .value("storage", ElementType::storage, "denotes a storage unit")
        .value("trafo_hv", ElementType::trafo_hv, "denotes a trafo (hv side)")
        .value("trafo_lv", ElementType::trafo_lv, "denotes a trafo (lv side)")
        .value("shunt", ElementType::shunt, "denotes a shunt")
        .value("static_gen", ElementType::static_gen, "denotes a static generator")
        .value("dc_line_or", ElementType::dc_line_or, "denotes a dc powerline (or side)")
        .value("dc_line_ex", ElementType::dc_line_ex, "denotes a dc powerline (ex side)")
        .export_values();

    // topology action
    py::class_<TopoAction>(m, "TopoAction",
                           "Limited topological action used for the light environment, with grid2op semantics: "
                           "`add_element` is grid2op `set_bus` (local busbar id: -1 disconnect, 0 nothing, 1..n_busbar_per_sub), "
                           "`set_line_status` is grid2op `set_line_status` (-1, 0, +1). Line ids follow grid2op numbering "
                           "(powerlines then transformers). `check_validity(grid)` validates the action and must be called "
                           "before it can be applied (LightEnv.init_actions does it).")
        .def(py::init<>())
        .def("add_element", &TopoAction::add_element, py::arg("el_type"), py::arg("el_id"), py::arg("local_bus_id"),
             "grid2op `set_bus` on one element: add_element(ElementType el_type, int el_id, int local_bus_id) with "
             "local_bus_id -1 (disconnect), 0 (nothing) or 1..n_busbar_per_sub.")
        .def("set_line_status", &TopoAction::set_line_status, py::arg("line_id"), py::arg("status"),
             "grid2op `set_line_status` on one line (grid2op numbering): -1 disconnect, 0 nothing, +1 reconnect.")
        .def("check_validity", &TopoAction::check_validity, py::arg("grid"),
             "Check the action against a grid (element exists, busbar exists, no contradiction) and resolve the busbars. "
             "Raises ValueError / IndexError on an invalid action.")
        .def("apply_to_gridmodel", &TopoAction::apply_to_gridmodel, py::arg("grid"),
             "Play the (checked) action on a grid: disconnect / reconnect / move its elements. The grid should be "
             "the one it was checked against, or a copy of it.")
        .def("is_do_nothing", &TopoAction::is_do_nothing, "True if the action does not modify anything")
        .def("has_been_checked", &TopoAction::has_been_checked, "True if check_validity was called since the last modification")
        .def_property_readonly("nb_set_bus", &TopoAction::nb_set_bus, "Number of (resolved) set_bus entries, after check_validity")
        .def_property_readonly("nb_set_line_status", &TopoAction::nb_set_line_status, "Number of (resolved) set_line_status entries, after check_validity")
        ;

    // protections
    py::class_<Protections>(m, "Protections", "Limited protection used for the light environment")
        .def(py::init<>())
        .def_property_readonly("total_time", &Protections::get_total_time, "Total time spent in the 'next_grid_state' function, cumulated over an episode")
        .def_property_readonly("overflow_time", &Protections::get_check_overflow_time, "Total time spent to check the overflow, cumulated over an episode")
        .def_property_readonly("powerflow_time", &Protections::get_powerflow_time, "Total time spent to perform the powerflow (init DC + AC), cumulated over an episode")
        .def_property_readonly("update_rho_time", &Protections::get_update_rho_time, "Total time spent to update the value of rho, cumulated over an episode")
        .def_property_readonly("nb_iter", &Protections::get_nb_iter, "Number of 'iterations' used to perform the 'next_grid_state' (number of powerflow run) for the last call to 'next_grid_state' only")
        .def_property_readonly("rho", &Protections::get_rho, "Last computed 'rho' values.")
        .def_property_readonly("line_disconnected", &Protections::get_line_disconnected, "Information on line / trafo disconnected")
        .def_property_readonly("line_time_step_overflow", &Protections::get_line_time_step_overflow, "For how long each line has been on overflow")
        .def_property_readonly("v_init_dc", &Protections::get_v_init_dc, "Complex voltage used to initialize the last DC powerflow")
        .def_property_readonly("v_init_ac", &Protections::get_v_init_ac, "Complex voltage used to initialize the last AC powerflow (result of a DC powerflow)")
        .def_property_readonly("v_res_ac", &Protections::get_v_res_ac, "Complex voltage at the end of the last powerflow")

        .def("lines_disconnected_this_step", &Protections::lines_disconnected_this_step, "Lines (grid2op numbering) disconnected by the protections during the last step, all iterations together")
        .def("get_thermal_limit_or", &Protections::get_thermal_limit_or,
             "Thermal limit (kA) of each line on its origin side (hv side for a trafo), grid2op numbering (powerlines then transformers)")
        .def("get_thermal_limit_ex", &Protections::get_thermal_limit_ex,
             "Thermal limit (kA) of each line on its extremity side (lv side for a trafo), grid2op numbering (powerlines then transformers)")
        .def("get_max_line_time_step_overflow", &Protections::get_max_line_time_step_overflow,
             "For each line, the number of steps it can stay in overflow: it is disconnected once its overflow counter "
             "exceeds this value (grid2op NB_TIMESTEP_OVERFLOW_ALLOWED)")
        .def("set_thermal_limit_or", &Protections::set_thermal_limit_or, py::arg("thermal_limit_or"),
             "Set the thermal limits (kA) on the origin side, one per line (powerlines then transformers). "
             "rho is the max over both sides of current / thermal limit.")
        .def("set_thermal_limit_ex", &Protections::set_thermal_limit_ex, py::arg("thermal_limit_ex"),
             "Set the thermal limits (kA) on the extremity side, one per line (powerlines then transformers). "
             "rho is the max over both sides of current / thermal limit.")
        .def("set_max_line_time_step_overflow", &Protections::set_max_line_time_step_overflow, py::arg("max_line_time_step_overflow"),
             "Set, for each line (powerlines then transformers), the number of steps it can stay in overflow before "
             "the protections disconnect it (int32 vector)")
        ;

    // observation: every array is a read-only numpy view on the env's memory (no copy), kept
    // alive by the observation, itself kept alive by its env (def_property_readonly is
    // reference_internal)
    py::class_<LightEnvObservation>(m, "LightEnvObservation",
        "Observation of a LightEnv: a read-only view on the current state of its environment, nothing is copied. "
        "Every attribute is a read-only numpy array on the environment's memory, so it follows the environment as "
        "it steps (it is not a snapshot: use np.array(obs.p_or) to keep one). An array is valid until the next reset "
        "of the environment or a step ending the episode by a divergence; once done is True the values are not "
        "meaningful. Lines use grid2op numbering (powerlines then transformers, 'or' being a transformer's hv side), "
        "powers are in MW / MVAr, currents in kA (like the thermal limits of Protections).")
        .def_property_readonly("rho", &LightEnvObservation::get_rho, "For each line, the max over its two sides of current / thermal limit")
        .def_property_readonly("p_or", &LightEnvObservation::get_p_or, "Active power flow (MW) at the origin side of each line")
        .def_property_readonly("q_or", &LightEnvObservation::get_q_or, "Reactive power flow (MVAr) at the origin side of each line")
        .def_property_readonly("a_or", &LightEnvObservation::get_a_or, "Current flow (kA) at the origin side of each line")
        .def_property_readonly("p_ex", &LightEnvObservation::get_p_ex, "Active power flow (MW) at the extremity side of each line")
        .def_property_readonly("q_ex", &LightEnvObservation::get_q_ex, "Reactive power flow (MVAr) at the extremity side of each line")
        .def_property_readonly("a_ex", &LightEnvObservation::get_a_ex, "Current flow (kA) at the extremity side of each line")
        .def_property_readonly("load_p", &LightEnvObservation::get_load_p, "Active power (MW) consumed by each load")
        .def_property_readonly("gen_p", &LightEnvObservation::get_gen_p, "Active power (MW) produced by each generator")
        .def_property_readonly("topo_vect", &LightEnvObservation::get_topo_vect,
                               "grid2op topology vector: the local busbar (1..n_busbar_per_sub, -1 if disconnected) of every "
                               "element, at its grid2op position. Raises if the grid carries no position in the topology "
                               "vector (they are set by LightSimBackend).")
        .def_property_readonly("time_before_cooldown_line", &LightEnvObservation::get_time_before_cooldown_line,
                               "For each line, number of steps before its status can be changed again")
        .def_property_readonly("time_before_cooldown_sub", &LightEnvObservation::get_time_before_cooldown_sub,
                               "For each substation, number of steps before it can be acted on again")
        .def_property_readonly("current_step", &LightEnvObservation::get_current_step, "Step of the environment (0 at reset)")
        ;

    // env
    py::class_<LightEnv>(m, "LightEnv", "Fast implementation of a grid2op env in pure c++ with (very) limited functionality")
        .def(py::init<const LSGrid &>())
        .def_property_readonly("step_time", &LightEnv::get_step_time, "Total time spent in the 'step' function, cumulated over the entire episode")
        .def_property_readonly("reset_time", &LightEnv::get_reset_time, "Total time spent in the 'reset' function, cumulated over the entire episode")
        .def_property_readonly("obs_time", &LightEnv::get_obs_time, "Total time spent to retrieve the observation, cumulated over the entire episode")
        .def_property_readonly("update_gridmodel_time", &LightEnv::get_update_gridmodel_time, "Total time spent to update the gridmodel (new injections), cumulated over an episode")
        .def_property_readonly("max_step", &LightEnv::get_max_step, "Maximum step for this environment (number of rows of the input matrices)")
        .def_property_readonly("current_step", &LightEnv::get_current_step, "Current step count for this environment (number of env.step). current_step is 0 at reset and increased by 1 for each env.step.")
        .def_property_readonly("grid", &LightEnv::get_grid, "State of the underlying powergrid (readonly)")
        .def_property("protections", &LightEnv::get_protections, &LightEnv::assign_protections, "The current protections of the env")
        .def_property("tol", &LightEnv::get_tol, &LightEnv::set_tol, "Tolerance used to compute the powerflows")
        .def_property("max_iter", &LightEnv::get_max_iter, &LightEnv::set_max_iter, "Maximum number of iterations for the powerflows")
        .def_property("nb_timestep_cooldown_sub", &LightEnv::get_nb_timestep_cooldown_sub, &LightEnv::set_nb_timestep_cooldown_sub,
                      "Number of steps a substation cannot be acted on after an action modified it (grid2op NB_TIMESTEP_COOLDOWN_SUB)")
        .def_property("nb_timestep_cooldown_line", &LightEnv::get_nb_timestep_cooldown_line, &LightEnv::set_nb_timestep_cooldown_line,
                      "Number of steps a line cannot be acted on after an action modified its status (grid2op NB_TIMESTEP_COOLDOWN_LINE)")
        .def_property("nb_timestep_reconnection", &LightEnv::get_nb_timestep_reconnection, &LightEnv::set_nb_timestep_reconnection,
                      "Number of steps a line disconnected by the protections cannot be reconnected (grid2op NB_TIMESTEP_RECONNECTION)")
        .def_property_readonly("time_before_cooldown_sub", &LightEnv::get_time_before_cooldown_sub,
                               "For each substation, number of steps before it can be acted on again (grid2op obs.time_before_cooldown_sub)")
        .def_property_readonly("time_before_cooldown_line", &LightEnv::get_time_before_cooldown_line,
                               "For each line (grid2op numbering), number of steps before its status can be changed again (grid2op obs.time_before_cooldown_line)")
        .def_property_readonly("nb_actions", &LightEnv::nb_actions, "Number of actions registered with init_actions")
        .def("init_actions", &LightEnv::init_actions, py::arg("actions"),
             "Register the actions the agent can take: step(i) plays actions[i]. Every action is checked against the "
             "initial grid and a ValueError naming the invalid action is raised (and nothing registered) if one is invalid.")
        .def("get_actions", &LightEnv::get_actions, "The (checked) actions registered with init_actions")
        .def("assign_time_series", &LightEnv::assign_time_series,
             py::arg("load_p"), py::arg("load_q"), py::arg("gen_p"), py::arg("gen_v"), py::arg("storage_p"),
             py::arg("shunt_p"), py::arg("shunt_q"), py::arg("sgen_p"), py::arg("sgen_q"),
             "Give the injections replayed by the episode: one matrix per quantity, one row per step and one column "
             "per element (MW, MVAr, gen_v in pu), NaN meaning 'unchanged'. All must have the same number of rows, "
             "which becomes max_step; the sizes are checked at the next reset. For now only load_p, load_q, gen_p "
             "and gen_v are applied to the grid, the other matrices are only checked.")
        .def("reset",
             [](py::object self){
                 LightEnv & env = self.cast<LightEnv &>();
                 const auto res = env.reset();
                 return py::make_tuple(py::cast(&std::get<0>(res), py::return_value_policy::reference_internal, self),
                                       std::get<1>(res));
             },
             "Start a new episode: restore the initial topology, the protections' counters and the cooldowns, apply "
             "the injections of the first row and run a powerflow. Returns (obs, info), obs being the LightEnvObservation of the env. Must be "
             "called before the first step, and again after assign_time_series or a change of protections.")
        .def("step",
             [](py::object self, int act_id){
                 LightEnv & env = self.cast<LightEnv &>();
                 const auto res = env.step(act_id);
                 return py::make_tuple(py::cast(&std::get<0>(res), py::return_value_policy::reference_internal, self),
                                       std::get<1>(res), std::get<2>(res), std::get<3>(res), std::get<4>(res));
             }, py::arg("act_id"),
             "Play one step with the action of id act_id (see init_actions). Without any action registered, only act_id = 0 (do nothing) is valid. "
             "Returns (obs, reward, done, truncated, info): obs is the LightEnvObservation of the env (the same object at every step, "
             "a view on the current state), the reward is the fraction of the episode survived, info['is_illegal'] is 'true' if the "
             "action was refused because of a cooldown.")
        .def("get_obs", &LightEnv::get_obs, py::return_value_policy::reference_internal,
             "The observation of the env (a LightEnvObservation, a view on its current state)")
        ;
}
