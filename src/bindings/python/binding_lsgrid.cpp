// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "pickle_helpers.hpp"
#include "LSGrid.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

void bind_gridmodel(py::module_& m) {
    auto lsgrid_cls = py::class_<LSGrid>(m, "LSGrid", DocLSGrid::LSGrid.c_str())
        .def(py::init<>())
        .def("copy", &LSGrid::copy, "TODO", py::return_value_policy::take_ownership)
        .def_property("_ls_to_orig",
                      &LSGrid::get_ls_to_orig,
                      &LSGrid::set_ls_to_orig,
                      R"mydelimiter(
_ls_to_orig: has the size of the number of possible buses in lightsim2grid
(*ie* `n_sub_ * max_nb_bus_per_sub_` ) and gives the id of the corresponding
bus in the original grid (pandapower or pypowsybl).

If a "-1" is present, then this bus does not exist in the original grid,
it is only present in the lightsim2grid gridmodel.
)mydelimiter")
        .def_property("_orig_to_ls",
                      &LSGrid::get_orig_to_ls,
                      &LSGrid::set_orig_to_ls,
                      R"mydelimiter(
Opposite to _ls_to_orig. The vector _orig_to_ls has the size of the number
of buses in the original grid (pandapower or pypowsybl) and tells
to which bus of lightsim2grid it corresponds. It should be a >= integer
between 0 and `n_sub_ * max_nb_bus_per_sub_`

)mydelimiter"
                    )
        .def_property("_max_nb_bus_per_sub",
                      &LSGrid::get_max_nb_bus_per_sub,
                      &LSGrid::set_max_nb_bus_per_sub,
                      "do not modify it after loading !")
        .def_property_readonly("timer_last_ac_pf", &LSGrid::timer_last_ac_pf, "TODO")
        .def_property_readonly("timer_last_dc_pf", &LSGrid::timer_last_dc_pf, "TODO");
    add_pickle(lsgrid_cls, "LSGrid");
    lsgrid_cls
        // algo config (scaling/refactor policy params)
        .def("get_ac_algo_config", &LSGrid::get_ac_algo_config,
            "Return the AC solver's AlgoConfig (scaling/refactor policy type and parameters).")
        .def("set_ac_algo_config", &LSGrid::set_ac_algo_config, py::arg("config"),
            "Apply an AlgoConfig to the AC solver (restores scaling/refactor policy and parameters).")
        .def("get_dc_algo_config", &LSGrid::get_dc_algo_config,
            "Return the DC solver's AlgoConfig (no-op for non-NR solvers, returns empty config).")
        .def("set_dc_algo_config", &LSGrid::set_dc_algo_config, py::arg("config"),
            "Apply an AlgoConfig to the DC solver.")

        // solver control
        .def("change_algorithm", py::overload_cast<const AlgorithmType&>(&LSGrid::change_algorithm), DocLSGrid::change_algorithm.c_str())
        .def("change_algorithm", py::overload_cast<const std::string&>(&LSGrid::change_algorithm), "Change the AC (or DC) algorithm by registry name. Accepts built-in names and plugin names registered via load_solver_plugin().")
        .def("available_default_algorithms", &LSGrid::available_default_algorithms, DocLSGrid::available_algorithm_names.c_str())
        .def("available_algorithm_names", &LSGrid::available_algorithm_names, "Returns names of all registered algorithms, including any loaded plugins.")
        .def("get_computation_time", &LSGrid::get_computation_time, DocLSGrid::get_computation_time.c_str())
        .def("get_dc_computation_time", &LSGrid::get_dc_computation_time, DocLSGrid::get_dc_computation_time.c_str())
        .def("get_algo_type", &LSGrid::get_algo_type, DocLSGrid::get_algo_type.c_str())
        .def("get_dc_algo_type", &LSGrid::get_dc_algo_type, DocLSGrid::get_dc_algo_type.c_str())
        .def("get_algo", &LSGrid::get_algo, py::return_value_policy::reference, DocLSGrid::get_algo.c_str())
        .def("get_dc_algo", &LSGrid::get_dc_algo, py::return_value_policy::reference, DocLSGrid::get_dc_algo.c_str())
        // deprecated method
        .def("change_solver", py::overload_cast<const AlgorithmType&>(&LSGrid::change_algorithm), "DEPRECATED: use 'change_algorithm' instead")
        .def("change_solver", py::overload_cast<const std::string&>(&LSGrid::change_algorithm), "DEPRECATED: use 'change_algorithm' instead")
        .def("available_solvers", &LSGrid::available_default_algorithms, "DEPRECATED: use 'available_default_algorithms' instead")
        .def("available_solver_names", &LSGrid::available_algorithm_names, "DEPRECATED: use 'available_algorithm_names' instead")
        .def("get_solver_type", &LSGrid::get_algo_type, "DEPRECATED: use 'get_algo_type' instead")
        .def("get_dc_solver_type", &LSGrid::get_dc_algo_type, "DEPRECATED: use 'get_dc_algo_type' instead")
        .def("get_solver", &LSGrid::get_algo, py::return_value_policy::reference, "DEPRECATED: use 'get_algo' instead")
        .def("get_dc_solver", &LSGrid::get_dc_algo, py::return_value_policy::reference, "DEPRECATED: use 'get_dc_algo' instead")

        // init the grid
        .def("init_bus", &LSGrid::init_bus, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_bus_status", &LSGrid::init_bus_status, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_init_vm_pu", &LSGrid::set_init_vm_pu, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_init_vm_pu", &LSGrid::get_init_vm_pu, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_sn_mva", &LSGrid::set_sn_mva, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_sn_mva", &LSGrid::get_sn_mva, DocLSGrid::_internal_do_not_use.c_str())

        // init elements
        .def("init_powerlines", &LSGrid::init_powerlines, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_powerlines_full", &LSGrid::init_powerlines_full, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_shunt", &LSGrid::init_shunt, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_trafo_pandapower", &LSGrid::init_trafo_pandapower, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_trafo", &LSGrid::init_trafo, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_generators", &LSGrid::init_generators, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_generators_full", &LSGrid::init_generators_full, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_loads", &LSGrid::init_loads, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_storages", &LSGrid::init_storages, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_sgens", &LSGrid::init_sgens, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_dclines", &LSGrid::init_dclines, DocLSGrid::_internal_do_not_use.c_str())
        .def("add_gen_slackbus", &LSGrid::add_gen_slackbus, DocLSGrid::_internal_do_not_use.c_str())
        .def("remove_gen_slackbus", &LSGrid::remove_gen_slackbus, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_bus_vn_kv", &LSGrid::get_bus_vn_kv, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_bus_status", &LSGrid::get_bus_status, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        // inspect the grid
        .def("get_substations", &LSGrid::get_substations, "TODO", py::return_value_policy::reference)
        .def("get_lines", &LSGrid::get_lines, DocLSGrid::get_lines.c_str(), py::return_value_policy::reference)
        .def("get_dclines", &LSGrid::get_dclines, DocLSGrid::get_dclines.c_str(), py::return_value_policy::reference)
        .def("get_trafos", &LSGrid::get_trafos, DocLSGrid::get_trafos.c_str(), py::return_value_policy::reference)
        .def("get_generators", &LSGrid::get_generators, DocLSGrid::get_generators.c_str(), py::return_value_policy::reference)
        .def("get_static_generators", &LSGrid::get_static_generators, DocLSGrid::get_static_generators.c_str(), py::return_value_policy::reference)
        .def("get_shunts", &LSGrid::get_shunts, DocLSGrid::get_shunts.c_str(), py::return_value_policy::reference)
        .def("get_storages", &LSGrid::get_storages, DocLSGrid::get_storages.c_str(), py::return_value_policy::reference)
        .def("get_loads", &LSGrid::get_loads, DocLSGrid::get_loads.c_str(), py::return_value_policy::reference)

        // pypowsybl compat names
        .def("get_voltage_levels", &LSGrid::get_substations, "TODO", py::return_value_policy::reference)
        .def("get_2_windings_transformers", &LSGrid::get_trafos, DocLSGrid::get_trafos.c_str(), py::return_value_policy::reference)
        .def("get_shunt_compensators", &LSGrid::get_shunts, DocLSGrid::get_shunts.c_str(), py::return_value_policy::reference)

        // modify the grid
        .def("turnedoff_no_pv", &LSGrid::turnedoff_no_pv, "Turned off (or generators with p = 0) generators will not be pv buses, they will not maintain voltage")
        .def("turnedoff_pv", &LSGrid::turnedoff_pv, "Turned off (or generators with p = 0) generators will be pv buses, they will maintain voltage (default)")
        .def("get_turnedoff_gen_pv", &LSGrid::get_turnedoff_gen_pv, "TODO")
        .def("update_slack_weights", &LSGrid::update_slack_weights, "TODO")
        .def("update_slack_weights_by_id", &LSGrid::update_slack_weights_by_id, "TODO")
        .def("assign_slack_to_most_connected", &LSGrid::assign_slack_to_most_connected, "TODO")
        .def("consider_only_main_component", &LSGrid::consider_only_main_component, "TODO and TODO DC LINE: one side might be in the connected comp and not the other !")
        .def("set_ignore_status_global", &LSGrid::set_ignore_status_global, "Ignore the 'global_status' flags for powerlines and trafo (set to true if you want to control independantly each side of powerlines and trafo). Default: false.")
        .def("set_synch_status_both_side", &LSGrid::set_synch_status_both_side, "Synch the status of each side of the powerlines and trafo. It means that if you disconnect one side of a powerline / trafo, the other side will also be disconnected. Default: true.")
        .def("get_ignore_status_global", &LSGrid::get_ignore_status_global, "TODO doc")
        .def("get_synch_status_both_side", &LSGrid::get_synch_status_both_side, "TODO doc")

        // names
        .def("set_line_names", &LSGrid::set_line_names, "TODO")
        .def("set_dcline_names", &LSGrid::set_dcline_names, "TODO")
        .def("set_trafo_names", &LSGrid::set_trafo_names, "TODO")
        .def("set_gen_names", &LSGrid::set_gen_names, "TODO")
        .def("set_load_names", &LSGrid::set_load_names, "TODO")
        .def("set_storage_names", &LSGrid::set_storage_names, "TODO")
        .def("set_sgen_names", &LSGrid::set_sgen_names, "TODO")
        .def("set_shunt_names", &LSGrid::set_shunt_names, "TODO")
        .def("set_substation_names", &LSGrid::set_substation_names, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_substation_names", &LSGrid::get_substation_names, DocLSGrid::_internal_do_not_use.c_str())

        .def("deactivate_bus", &LSGrid::deactivate_bus_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("reactivate_bus", &LSGrid::reactivate_bus_python, DocLSGrid::_internal_do_not_use.c_str())

        .def("deactivate_powerline", &LSGrid::deactivate_powerline, DocLSGrid::_internal_do_not_use.c_str())
        .def("reactivate_powerline", &LSGrid::reactivate_powerline, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus1_powerline", &LSGrid::change_bus1_powerline_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus2_powerline", &LSGrid::change_bus2_powerline_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_bus1_powerline", &LSGrid::get_bus1_powerline, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_bus2_powerline", &LSGrid::get_bus2_powerline, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        .def("deactivate_trafo", &LSGrid::deactivate_trafo, DocLSGrid::_internal_do_not_use.c_str())
        .def("reactivate_trafo", &LSGrid::reactivate_trafo, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus1_trafo", &LSGrid::change_bus1_trafo_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus2_trafo", &LSGrid::change_bus2_trafo_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_bus1_trafo", &LSGrid::get_bus1_trafo, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_bus2_trafo", &LSGrid::get_bus2_trafo, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_ratio_trafo", &LSGrid::change_ratio_trafo, "TODO")
        .def("change_shift_trafo", &LSGrid::change_shift_trafo,
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
        .def("change_shift_trafo_deg", &LSGrid::change_shift_trafo_deg,
            "Same as :ref:`change_shift_trafo` but phase shift is expressed in degree and NOT in rad.")
        .def("deactivate_load", &LSGrid::deactivate_load, DocLSGrid::_internal_do_not_use.c_str())
        .def("reactivate_load", &LSGrid::reactivate_load, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus_load", &LSGrid::change_bus_load_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_bus_load", &LSGrid::get_bus_load, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_p_load", &LSGrid::change_p_load, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_q_load", &LSGrid::change_q_load, DocLSGrid::_internal_do_not_use.c_str())

        .def("deactivate_gen", &LSGrid::deactivate_gen, DocLSGrid::_internal_do_not_use.c_str())
        .def("reactivate_gen", &LSGrid::reactivate_gen, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus_gen", &LSGrid::change_bus_gen_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_bus_gen", &LSGrid::get_bus_gen, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_p_gen", &LSGrid::change_p_gen, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_v_gen", &LSGrid::change_v_gen, DocLSGrid::_internal_do_not_use.c_str())

        .def("deactivate_shunt", &LSGrid::deactivate_shunt, DocLSGrid::_internal_do_not_use.c_str())
        .def("reactivate_shunt", &LSGrid::reactivate_shunt, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus_shunt", &LSGrid::change_bus_shunt_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_bus_shunt", &LSGrid::get_bus_shunt, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_p_shunt", &LSGrid::change_p_shunt, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_q_shunt", &LSGrid::change_q_shunt, DocLSGrid::_internal_do_not_use.c_str())

        .def("deactivate_sgen", &LSGrid::deactivate_sgen, DocLSGrid::_internal_do_not_use.c_str())
        .def("reactivate_sgen", &LSGrid::reactivate_sgen, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus_sgen", &LSGrid::change_bus_sgen_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_bus_sgen", &LSGrid::get_bus_sgen, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_p_sgen", &LSGrid::change_p_sgen, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_q_sgen", &LSGrid::change_q_sgen, DocLSGrid::_internal_do_not_use.c_str())

        .def("deactivate_storage", &LSGrid::deactivate_storage, DocLSGrid::_internal_do_not_use.c_str())
        .def("reactivate_storage", &LSGrid::reactivate_storage, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus_storage", &LSGrid::change_bus_storage_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_bus_storage", &LSGrid::get_bus_storage, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("change_p_storage", &LSGrid::change_p_storage, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_q_storage", &LSGrid::change_q_storage, DocLSGrid::_internal_do_not_use.c_str())

        .def("deactivate_dcline", &LSGrid::deactivate_dcline, DocLSGrid::_internal_do_not_use.c_str())
        .def("reactivate_dcline", &LSGrid::reactivate_dcline, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_p_dcline", &LSGrid::change_p_dcline, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_v1_dcline", &LSGrid::change_v1_dcline, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_v2_dcline", &LSGrid::change_v2_dcline, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus1_dcline", &LSGrid::change_bus1_dcline, DocLSGrid::_internal_do_not_use.c_str())
        .def("change_bus2_dcline", &LSGrid::change_bus2_dcline, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_bus1_dcline", &LSGrid::get_bus1_dcline, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_bus2_dcline", &LSGrid::get_bus2_dcline, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        // get back the results
        .def("get_V", &LSGrid::get_V, DocLSGrid::get_V.c_str())
        .def("get_Va", &LSGrid::get_Va, DocLSGrid::get_Va.c_str())
        .def("get_Vm", &LSGrid::get_Vm, DocLSGrid::get_Vm.c_str())
        .def("get_V_solver", &LSGrid::get_V_solver, DocLSGrid::get_V_solver.c_str(), py::return_value_policy::reference)
        .def("get_Va_solver", &LSGrid::get_Va_solver, DocLSGrid::get_Va_solver.c_str(), py::return_value_policy::reference)
        .def("get_Vm_solver", &LSGrid::get_Vm_solver, DocLSGrid::get_Vm_solver.c_str(), py::return_value_policy::reference)
        .def("get_J_solver", &LSGrid::get_J_python_solver, DocLSGrid::get_J_python_solver.c_str(), py::return_value_policy::reference)

        .def("id_me_to_ac_solver", &LSGrid::id_ac_solver_to_me_numpy, DocLSGrid::id_me_to_ac_solver.c_str(), py::return_value_policy::reference)
        .def("id_ac_solver_to_me", &LSGrid::id_ac_solver_to_me_numpy, DocLSGrid::id_ac_solver_to_me.c_str(), py::return_value_policy::reference)
        .def("id_me_to_dc_solver", &LSGrid::id_me_to_dc_solver_numpy, DocLSGrid::id_me_to_dc_solver.c_str(), py::return_value_policy::reference)
        .def("id_dc_solver_to_me", &LSGrid::id_dc_solver_to_me_numpy, DocLSGrid::id_dc_solver_to_me.c_str(), py::return_value_policy::reference)
        .def("total_bus", &LSGrid::total_bus, DocLSGrid::total_bus.c_str())
        .def("nb_connected_bus", &LSGrid::nb_connected_bus, DocLSGrid::nb_connected_bus.c_str())

        .def("get_pv", &LSGrid::get_pv_numpy, DocLSGrid::get_pv.c_str(), py::return_value_policy::reference)
        .def("get_pq", &LSGrid::get_pq_numpy, DocLSGrid::get_pq.c_str(), py::return_value_policy::reference)
        .def("get_slack_ids", &LSGrid::get_slack_ids_numpy, DocLSGrid::get_slack_ids.c_str(), py::return_value_policy::reference)
        .def("get_slack_ids_dc", &LSGrid::get_slack_ids_dc_numpy, DocLSGrid::get_slack_ids_dc.c_str(), py::return_value_policy::reference)
        .def("get_slack_weights", &LSGrid::get_slack_weights, DocLSGrid::get_slack_weights.c_str(), py::return_value_policy::reference)
        .def("get_pv_solver", &LSGrid::get_pv_solver_numpy, DocLSGrid::get_pv_solver.c_str(), py::return_value_policy::reference)
        .def("get_pq_solver", &LSGrid::get_pq_solver_numpy, DocLSGrid::get_pq_solver.c_str(), py::return_value_policy::reference)
        .def("get_slack_ids_solver", &LSGrid::get_slack_ids_solver_numpy, DocLSGrid::get_slack_ids_solver.c_str(), py::return_value_policy::reference)
        .def("get_slack_ids_dc_solver", &LSGrid::get_slack_ids_dc_solver_numpy, DocLSGrid::get_slack_ids_dc_solver.c_str(), py::return_value_policy::reference)
        .def("get_slack_weights_solver", &LSGrid::get_slack_weights_solver, DocLSGrid::get_slack_weights_solver.c_str(), py::return_value_policy::reference)

        .def("get_Ybus", &LSGrid::get_Ybus, DocLSGrid::get_Ybus.c_str())
        .def("get_dcYbus", &LSGrid::get_dcYbus, DocLSGrid::get_dcYbus.c_str())
        .def("get_Sbus", &LSGrid::get_Sbus, DocLSGrid::get_Sbus.c_str())
        .def("get_dcSbus", &LSGrid::get_dcSbus, DocLSGrid::get_dcSbus.c_str())
        .def("get_Ybus_solver", &LSGrid::get_Ybus_solver, DocLSGrid::get_Ybus_solver.c_str(), py::return_value_policy::reference)
        .def("get_dcYbus_solver", &LSGrid::get_dcYbus_solver, DocLSGrid::get_dcYbus_solver.c_str(), py::return_value_policy::reference)
        .def("get_Sbus_solver", &LSGrid::get_Sbus_solver, DocLSGrid::get_Sbus_solver.c_str(), py::return_value_policy::reference)
        .def("get_dcSbus_solver", &LSGrid::get_dcSbus_solver, DocLSGrid::get_dcSbus_solver.c_str(), py::return_value_policy::reference)

        .def("check_solution", &LSGrid::check_solution, DocLSGrid::check_solution.c_str())

        .def("get_loads_res", &LSGrid::get_loads_res, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_loads_status", &LSGrid::get_loads_status, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_shunts_res", &LSGrid::get_shunts_res, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_shunts_status", &LSGrid::get_shunts_status, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_gen_res", &LSGrid::get_gen_res, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_gen_status", &LSGrid::get_gen_status, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_res1", &LSGrid::get_line_res1, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_res2", &LSGrid::get_line_res2, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_lines_status", &LSGrid::get_lines_status, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res1", &LSGrid::get_trafo_res1, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res2", &LSGrid::get_trafo_res2, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_status", &LSGrid::get_trafo_status, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_storages_res", &LSGrid::get_storages_res, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_storages_status", &LSGrid::get_storages_status, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_sgens_res", &LSGrid::get_sgens_res, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_sgens_status", &LSGrid::get_sgens_status, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        .def("get_gen_theta", &LSGrid::get_gen_theta, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_load_theta", &LSGrid::get_load_theta, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_shunt_theta", &LSGrid::get_shunt_theta, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_storage_theta", &LSGrid::get_storage_theta, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_theta1", &LSGrid::get_line_theta1, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_theta2", &LSGrid::get_line_theta2, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_theta1", &LSGrid::get_trafo_theta1, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_theta2", &LSGrid::get_trafo_theta2, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        .def("get_all_shunt_buses", &LSGrid::get_all_shunt_buses_numpy, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_loads_res_full", &LSGrid::get_loads_res_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_shunts_res_full", &LSGrid::get_shunts_res_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_gen_res_full", &LSGrid::get_gen_res_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_res1_full", &LSGrid::get_line_res1_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_line_res2_full", &LSGrid::get_line_res2_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res1_full", &LSGrid::get_trafo_res1_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res2_full", &LSGrid::get_trafo_res2_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_storages_res_full", &LSGrid::get_storages_res_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_sgens_res_full", &LSGrid::get_sgens_res_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_dcline_res1_full", &LSGrid::get_dcline_res1_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_dcline_res2_full", &LSGrid::get_dcline_res2_full, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        .def("get_shunt_target_p", &LSGrid::get_shunt_target_p, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_load_target_p", &LSGrid::get_load_target_p, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_gen_target_p", &LSGrid::get_gen_target_p, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_sgen_target_p", &LSGrid::get_sgen_target_p, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)
        .def("get_storage_target_p", &LSGrid::get_storage_target_p, DocLSGrid::_internal_do_not_use.c_str(), py::return_value_policy::reference)

        // do something with the grid
        .def("deactivate_result_computation", &LSGrid::deactivate_result_computation, DocLSGrid::deactivate_result_computation.c_str())
        .def("reactivate_result_computation", &LSGrid::reactivate_result_computation, DocLSGrid::reactivate_result_computation.c_str())
        .def("dc_pf", &LSGrid::dc_pf, DocLSGrid::dc_pf.c_str())
        .def("ac_pf", &LSGrid::ac_pf, DocLSGrid::ac_pf.c_str())
        .def("unset_changes", &LSGrid::unset_changes, DocLSGrid::_internal_do_not_use.c_str())
        .def("tell_recompute_ybus", &LSGrid::tell_recompute_ybus, DocLSGrid::_internal_do_not_use.c_str())
        .def("tell_recompute_sbus", &LSGrid::tell_recompute_sbus, DocLSGrid::_internal_do_not_use.c_str())
        .def("tell_solver_need_reset", &LSGrid::tell_solver_need_reset, DocLSGrid::_internal_do_not_use.c_str())
        .def("tell_ybus_change_sparsity_pattern", &LSGrid::tell_ybus_change_sparsity_pattern, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_algo_controler", &LSGrid::get_algo_controler, "TODO", py::return_value_policy::reference)
        .def("get_solver_control",  &LSGrid::get_algo_controler, "DEPRECATED use 'get_algo_controler'", py::return_value_policy::reference)
        .def("compute_newton", &LSGrid::ac_pf, DocLSGrid::ac_pf.c_str())
        .def("get_ptdf", &LSGrid::get_ptdf, DocLSGrid::get_ptdf.c_str(), py::return_value_policy::reference)
        .def("get_ptdf_solver", &LSGrid::get_ptdf_solver, DocLSGrid::get_ptdf_solver.c_str(), py::return_value_policy::reference)
        .def("get_lodf", &LSGrid::get_lodf, DocLSGrid::get_lodf.c_str(), py::return_value_policy::reference)
        .def("get_Bf", &LSGrid::get_Bf, DocLSGrid::get_Bf.c_str(), py::return_value_policy::reference)
        .def("get_Bf_solver", &LSGrid::get_Bf_solver, DocLSGrid::get_Bf_solver.c_str(), py::return_value_policy::reference)

        // apply action faster (optimized for grid2op representation)
        .def("update_gens_p", &LSGrid::update_gens_p, DocLSGrid::_internal_do_not_use.c_str())
        .def("update_sgens_p", &LSGrid::update_sgens_p, DocLSGrid::_internal_do_not_use.c_str())
        .def("update_gens_v", &LSGrid::update_gens_v, DocLSGrid::_internal_do_not_use.c_str())
        .def("update_loads_p", &LSGrid::update_loads_p, DocLSGrid::_internal_do_not_use.c_str())
        .def("update_loads_q", &LSGrid::update_loads_q, DocLSGrid::_internal_do_not_use.c_str())
        .def("update_topo", &LSGrid::update_topo, DocLSGrid::_internal_do_not_use.c_str())
        .def("update_storages_p", &LSGrid::update_storages_p, DocLSGrid::_internal_do_not_use.c_str())

        // auxiliary functions
        .def("set_n_sub", &LSGrid::set_n_sub, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_n_sub", &LSGrid::get_n_sub, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_max_nb_bus_per_sub", &LSGrid::set_max_nb_bus_per_sub, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_load_pos_topo_vect", &LSGrid::set_load_pos_topo_vect, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_gen_pos_topo_vect", &LSGrid::set_gen_pos_topo_vect, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_line_pos1_topo_vect", &LSGrid::set_line_pos1_topo_vect, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_line_pos2_topo_vect", &LSGrid::set_line_pos2_topo_vect, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_trafo_pos1_topo_vect", &LSGrid::set_trafo_pos1_topo_vect, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_trafo_pos2_topo_vect", &LSGrid::set_trafo_pos2_topo_vect, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_storage_pos_topo_vect", &LSGrid::set_storage_pos_topo_vect, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_load_to_subid", &LSGrid::set_load_to_subid, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_gen_to_subid", &LSGrid::set_gen_to_subid, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_shunt_to_subid", &LSGrid::set_shunt_to_subid, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_line_to_sub1_id", &LSGrid::set_line_to_sub1_id, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_line_to_sub2_id", &LSGrid::set_line_to_sub2_id, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_trafo_to_sub1_id", &LSGrid::set_trafo_to_sub1_id, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_trafo_to_sub2_id", &LSGrid::set_trafo_to_sub2_id, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_storage_to_subid", &LSGrid::set_storage_to_subid, DocLSGrid::_internal_do_not_use.c_str())

        // debug functions (might disappear without further notice)
        .def("debug_get_Bp_python", &LSGrid::debug_get_Bp_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("debug_get_Bpp_python", &LSGrid::debug_get_Bpp_python, DocLSGrid::_internal_do_not_use.c_str());
}
