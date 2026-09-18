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
    py::class_<TopoAction>(m, "TopoAction", "Limited topological action used for the light environment")
        .def(py::init<>())
        .def("add_element", &TopoAction::add_element, "allows to add an element that will be affected by this topological action (add_element(ElementType el_type, int el_id, int bus_global_id))")
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

        .def("get_thermal_limit_or", &Protections::get_thermal_limit_or, "TODO")
        .def("get_thermal_limit_ex", &Protections::get_thermal_limit_ex, "TODO")
        .def("get_max_line_time_step_overflow", &Protections::get_max_line_time_step_overflow, "TODO")
        .def("set_thermal_limit_or", &Protections::set_thermal_limit_or, "TODO")
        .def("set_thermal_limit_ex", &Protections::set_thermal_limit_ex, "TODO")
        .def("set_max_line_time_step_overflow", &Protections::set_max_line_time_step_overflow, "TODO")
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
        .def("assign_time_series", &LightEnv::assign_time_series, "TODO")
        .def("reset", &LightEnv::reset, "TODO")
        .def("step", &LightEnv::step, "TODO")
        .def("get_obs", &LightEnv::get_obs, "Last computed observation")
        ;
}
