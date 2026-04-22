// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "pickle_helpers.hpp"
#include "GridModel.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

void bind_gridmodel(py::module_& m) {
    auto gridmodel_cls = py::class_<GridModel>(m, "GridModel", DocGridModel::GridModel.c_str())
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
        .def_property_readonly("timer_last_dc_pf", &GridModel::timer_last_dc_pf, "TODO");
    add_pickle(gridmodel_cls, "GridModel");
    gridmodel_cls
        // solver control
        .def("change_solver", py::overload_cast<const SolverType&>(&GridModel::change_solver), DocGridModel::change_solver.c_str())
        .def("change_solver", py::overload_cast<const std::string&>(&GridModel::change_solver), "Change the AC (or DC) solver by registry name. Accepts built-in names and plugin names registered via load_solver_plugin().")
        .def("available_solvers", &GridModel::available_solvers, DocGridModel::available_solvers.c_str())
        .def("available_solver_names", &GridModel::available_solver_names, "Returns names of all registered solvers, including any loaded plugins.")
        .def("get_computation_time", &GridModel::get_computation_time, DocGridModel::get_computation_time.c_str())
        .def("get_dc_computation_time", &GridModel::get_dc_computation_time, DocGridModel::get_dc_computation_time.c_str())
        .def("get_solver_type", &GridModel::get_solver_type, DocGridModel::get_solver_type.c_str())
        .def("get_dc_solver_type", &GridModel::get_dc_solver_type, DocGridModel::get_dc_solver_type.c_str())
        .def("get_solver", &GridModel::get_solver, py::return_value_policy::reference, DocGridModel::get_solver.c_str())
        .def("get_dc_solver", &GridModel::get_dc_solver, py::return_value_policy::reference, DocGridModel::get_dc_solver.c_str())

        // init the grid
        .def("init_bus", &GridModel::init_bus, DocGridModel::_internal_do_not_use.c_str())
        .def("init_bus_status", &GridModel::init_bus_status, DocGridModel::_internal_do_not_use.c_str())
        .def("set_init_vm_pu", &GridModel::set_init_vm_pu, DocGridModel::_internal_do_not_use.c_str())
        .def("get_init_vm_pu", &GridModel::get_init_vm_pu, DocGridModel::_internal_do_not_use.c_str())
        .def("set_sn_mva", &GridModel::set_sn_mva, DocGridModel::_internal_do_not_use.c_str())
        .def("get_sn_mva", &GridModel::get_sn_mva, DocGridModel::_internal_do_not_use.c_str())

        // init elements
        .def("init_powerlines", &GridModel::init_powerlines, DocGridModel::_internal_do_not_use.c_str())
        .def("init_powerlines_full", &GridModel::init_powerlines_full, DocGridModel::_internal_do_not_use.c_str())
        .def("init_shunt", &GridModel::init_shunt, DocGridModel::_internal_do_not_use.c_str())
        .def("init_trafo_pandapower", &GridModel::init_trafo_pandapower, DocGridModel::_internal_do_not_use.c_str())
        .def("init_trafo", &GridModel::init_trafo, DocGridModel::_internal_do_not_use.c_str())
        .def("init_generators", &GridModel::init_generators, DocGridModel::_internal_do_not_use.c_str())
        .def("init_generators_full", &GridModel::init_generators_full, DocGridModel::_internal_do_not_use.c_str())
        .def("init_loads", &GridModel::init_loads, DocGridModel::_internal_do_not_use.c_str())
        .def("init_storages", &GridModel::init_storages, DocGridModel::_internal_do_not_use.c_str())
        .def("init_sgens", &GridModel::init_sgens, DocGridModel::_internal_do_not_use.c_str())
        .def("init_dclines", &GridModel::init_dclines, DocGridModel::_internal_do_not_use.c_str())
        .def("add_gen_slackbus", &GridModel::add_gen_slackbus, DocGridModel::_internal_do_not_use.c_str())
        .def("remove_gen_slackbus", &GridModel::remove_gen_slackbus, DocGridModel::_internal_do_not_use.c_str())
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
        .def("set_storage_names", &GridModel::set_storage_names, "TODO")
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
        .def("get_V", &GridModel::get_V, DocGridModel::get_V.c_str())
        .def("get_Va", &GridModel::get_Va, DocGridModel::get_Va.c_str())
        .def("get_Vm", &GridModel::get_Vm, DocGridModel::get_Vm.c_str())
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

        .def("get_Ybus", &GridModel::get_Ybus, DocGridModel::get_Ybus.c_str())
        .def("get_dcYbus", &GridModel::get_dcYbus, DocGridModel::get_dcYbus.c_str())
        .def("get_Sbus", &GridModel::get_Sbus, DocGridModel::get_Sbus.c_str())
        .def("get_dcSbus", &GridModel::get_dcSbus, DocGridModel::get_dcSbus.c_str())
        .def("get_Ybus_solver", &GridModel::get_Ybus_solver, DocGridModel::get_Ybus_solver.c_str(), py::return_value_policy::reference)
        .def("get_dcYbus_solver", &GridModel::get_dcYbus_solver, DocGridModel::get_dcYbus_solver.c_str(), py::return_value_policy::reference)
        .def("get_Sbus_solver", &GridModel::get_Sbus_solver, DocGridModel::get_Sbus_solver.c_str(), py::return_value_policy::reference)
        .def("get_dcSbus_solver", &GridModel::get_dcSbus_solver, DocGridModel::get_dcSbus_solver.c_str(), py::return_value_policy::reference)

        .def("check_solution", &GridModel::check_solution, DocGridModel::check_solution.c_str())

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

        // debug functions (might disappear without further notice)
        .def("debug_get_Bp_python", &GridModel::debug_get_Bp_python, DocGridModel::_internal_do_not_use.c_str())
        .def("debug_get_Bpp_python", &GridModel::debug_get_Bpp_python, DocGridModel::_internal_do_not_use.c_str());
}
