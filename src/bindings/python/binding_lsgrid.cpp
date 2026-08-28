// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "pickle_helpers.hpp"
#include "binary_helpers.hpp"
#include "LSGrid.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

void bind_gridmodel(py::module_& m) {
    auto lsgrid_cls = py::class_<LSGrid>(m, "LSGrid", DocLSGrid::LSGrid.c_str())
        .def(py::init<>())
        .def("copy", &LSGrid::copy, DocLSGrid::copy.c_str(), py::return_value_policy::take_ownership)
        .def_property("_ls_to_orig",
                      &LSGrid::get_ls_to_orig,
                      &LSGrid::set_ls_to_orig,
                      DocLSGrid::_ls_to_orig.c_str())
        .def_property("_orig_to_ls",
                      &LSGrid::get_orig_to_ls,
                      &LSGrid::set_orig_to_ls,
                      DocLSGrid::_orig_to_ls.c_str())
        .def_property("_max_nb_bus_per_sub",
                      &LSGrid::get_max_nb_bus_per_sub,
                      &LSGrid::set_max_nb_bus_per_sub,
                      DocLSGrid::_max_nb_bus_per_sub.c_str())
        .def_property("_init_kwargs",
                      &LSGrid::get_init_kwargs,
                      &LSGrid::set_init_kwargs,
                      DocLSGrid::_init_kwargs.c_str())
        .def_property("_bus_fusion_rep",
                      &LSGrid::get_bus_fusion_rep,
                      &LSGrid::set_bus_fusion_rep,
                      DocLSGrid::_bus_fusion_rep.c_str())
        .def_property_readonly("timer_last_ac_pf", &LSGrid::timer_last_ac_pf, DocLSGrid::timer_last_ac_pf.c_str())
        .def_property_readonly("timer_last_dc_pf", &LSGrid::timer_last_dc_pf, DocLSGrid::timer_last_dc_pf.c_str());
    add_pickle(lsgrid_cls, "LSGrid");
    add_binary_serialization(lsgrid_cls);
    // Companion to load_binary(): loads the grid *data* without re-selecting the
    // solver it was saved with. Lets a grid saved with a solver that is not
    // available here (a plugin that has not been loaded, or a built-in needing an
    // optional backend this build lacks) still be loaded -- the grid keeps the
    // default solvers, and you pick one yourself with change_algorithm().
    lsgrid_cls.def_static("load_binary_without_algorithm", [](const std::string& path) {
        return ls2g::load_binary_generic_with<LSGrid>(
            path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, /*restore_algorithm=*/false);
    }, py::arg("path"),
       "Load a grid saved with save_binary(), WITHOUT restoring the AC / DC solver it "
       "was saved with (nor that solver's configuration): the grid keeps the default "
       "solvers and you select one yourself with change_algorithm(). Use this when "
       "load_binary() reports that the saved solver is unavailable here -- typically a "
       "solver plugin that has not been loaded in this process. Every other check "
       "(binary format, corruption, grid consistency) is applied exactly as in "
       "load_binary().");
    lsgrid_cls
        // whole-grid consistency validation
        .def("check_grid", &LSGrid::check_grid, DocLSGrid::check_grid.c_str())

        // algo config (scaling/refactor policy params)
        .def("get_ac_algo_config", &LSGrid::get_ac_algo_config, DocLSGrid::get_ac_algo_config.c_str())
        .def("set_ac_algo_config", &LSGrid::set_ac_algo_config, py::arg("config"), DocLSGrid::set_ac_algo_config.c_str())
        .def("get_dc_algo_config", &LSGrid::get_dc_algo_config, DocLSGrid::get_dc_algo_config.c_str())
        .def("set_dc_algo_config", &LSGrid::set_dc_algo_config, py::arg("config"), DocLSGrid::set_dc_algo_config.c_str())

        // solver control
        .def("change_algorithm", py::overload_cast<const AlgorithmType&>(&LSGrid::change_algorithm), DocLSGrid::change_algorithm.c_str())
        .def("change_algorithm", py::overload_cast<const std::string&>(&LSGrid::change_algorithm), DocLSGrid::change_algorithm_by_name.c_str())
        .def("available_default_algorithms", &LSGrid::available_default_algorithms, DocLSGrid::available_default_algorithms.c_str())
        .def("available_algorithm_names", &LSGrid::available_algorithm_names, DocLSGrid::available_algorithm_names.c_str())
        .def("get_computation_time", &LSGrid::get_computation_time, DocLSGrid::get_computation_time.c_str())
        .def("get_dc_computation_time", &LSGrid::get_dc_computation_time, DocLSGrid::get_dc_computation_time.c_str())
        .def("get_algo_type", &LSGrid::get_algo_type, DocLSGrid::get_algo_type.c_str())
        .def("get_dc_algo_type", &LSGrid::get_dc_algo_type, DocLSGrid::get_dc_algo_type.c_str())
        .def("get_algo", &LSGrid::get_algo, py::return_value_policy::reference_internal, DocLSGrid::get_algo.c_str())
        .def("get_dc_algo", &LSGrid::get_dc_algo, py::return_value_policy::reference_internal, DocLSGrid::get_dc_algo.c_str())
        // deprecated method
        .def("change_solver", py::overload_cast<const AlgorithmType&>(&LSGrid::change_algorithm), "DEPRECATED: use 'change_algorithm' instead")
        .def("change_solver", py::overload_cast<const std::string&>(&LSGrid::change_algorithm), "DEPRECATED: use 'change_algorithm' instead")
        .def("available_solvers", &LSGrid::available_default_algorithms, "DEPRECATED: use 'available_default_algorithms' instead")
        .def("available_solver_names", &LSGrid::available_algorithm_names, "DEPRECATED: use 'available_algorithm_names' instead")
        .def("get_solver_type", &LSGrid::get_algo_type, "DEPRECATED: use 'get_algo_type' instead")
        .def("get_dc_solver_type", &LSGrid::get_dc_algo_type, "DEPRECATED: use 'get_dc_algo_type' instead")
        .def("get_solver", &LSGrid::get_algo, py::return_value_policy::reference_internal, "DEPRECATED: use 'get_algo' instead")
        .def("get_dc_solver", &LSGrid::get_dc_algo, py::return_value_policy::reference_internal, "DEPRECATED: use 'get_dc_algo' instead")

        // init the grid
        .def("init_bus", &LSGrid::init_bus, DocLSGrid::_internal_do_not_use.c_str())
        .def("init_bus_status", &LSGrid::init_bus_status, DocLSGrid::_internal_do_not_use.c_str())
        .def("set_init_vm_pu", &LSGrid::set_init_vm_pu, DocLSGrid::set_init_vm_pu.c_str())
        .def("get_init_vm_pu", &LSGrid::get_init_vm_pu, DocLSGrid::get_init_vm_pu.c_str())
        .def("set_sn_mva", &LSGrid::set_sn_mva, DocLSGrid::set_sn_mva.c_str())
        .def("get_sn_mva", &LSGrid::get_sn_mva, DocLSGrid::get_sn_mva.c_str())

        // init elements
        .def("init_powerlines", &LSGrid::init_powerlines, DocLSGrid::init_powerlines.c_str())
        .def("init_powerlines_full", &LSGrid::init_powerlines_full, DocLSGrid::init_powerlines_full.c_str())
        .def("init_shunt", &LSGrid::init_shunt, DocLSGrid::init_shunt.c_str())
        .def("init_trafo_pandapower", &LSGrid::init_trafo_pandapower, DocLSGrid::init_trafo_pandapower.c_str())
        .def("init_trafo", &LSGrid::init_trafo, DocLSGrid::init_trafo.c_str())
        .def("init_generators", &LSGrid::init_generators, DocLSGrid::init_generators.c_str())
        .def("init_generators_full", &LSGrid::init_generators_full, DocLSGrid::init_generators_full.c_str())
        .def("init_loads", &LSGrid::init_loads, DocLSGrid::init_loads.c_str())
        .def("init_storages", &LSGrid::init_storages, DocLSGrid::init_storages.c_str())
        .def("init_sgens", &LSGrid::init_sgens, DocLSGrid::init_sgens.c_str())
        .def("init_dclines", &LSGrid::init_dclines, DocLSGrid::init_dclines.c_str())
        .def("init_hvdc_lines", &LSGrid::init_hvdc_lines, DocLSGrid::init_hvdc_lines.c_str())
        .def("init_svcs", &LSGrid::init_svcs, DocLSGrid::init_svcs.c_str())
        .def("add_gen_slackbus", &LSGrid::add_gen_slackbus, DocLSGrid::add_gen_slackbus.c_str())
        .def("remove_gen_slackbus", &LSGrid::remove_gen_slackbus, DocLSGrid::remove_gen_slackbus.c_str())
        .def("get_bus_vn_kv", &LSGrid::get_bus_vn_kv, DocLSGrid::get_bus_vn_kv.c_str(), py::return_value_policy::reference_internal)
        .def("get_bus_status", &LSGrid::get_bus_status, DocLSGrid::get_bus_status.c_str(), py::return_value_policy::reference)
        .def("set_bus_voltage_limits", &LSGrid::set_bus_voltage_limits, DocLSGrid::set_bus_voltage_limits.c_str())
        .def("get_bus_vmin_kv", &LSGrid::get_bus_vmin_kv, DocLSGrid::get_bus_vmin_kv.c_str(), py::return_value_policy::reference_internal)
        .def("get_bus_vmax_kv", &LSGrid::get_bus_vmax_kv, DocLSGrid::get_bus_vmax_kv.c_str(), py::return_value_policy::reference_internal)

        // inspect the grid
        .def("get_substations", &LSGrid::get_substations, DocLSGrid::get_substations.c_str(), py::return_value_policy::reference_internal)
        .def("get_lines", &LSGrid::get_lines, DocLSGrid::get_lines.c_str(), py::return_value_policy::reference_internal)
        .def("get_dclines", &LSGrid::get_dclines, DocLSGrid::get_dclines.c_str(), py::return_value_policy::reference_internal)
        .def("get_trafos", &LSGrid::get_trafos, DocLSGrid::get_trafos.c_str(), py::return_value_policy::reference_internal)
        .def("get_generators", &LSGrid::get_generators, DocLSGrid::get_generators.c_str(), py::return_value_policy::reference_internal)
        .def("get_static_generators", &LSGrid::get_static_generators, DocLSGrid::get_static_generators.c_str(), py::return_value_policy::reference_internal)
        .def("get_svcs", &LSGrid::get_svcs, DocLSGrid::get_svcs.c_str(), py::return_value_policy::reference_internal)
        .def("get_shunts", &LSGrid::get_shunts, DocLSGrid::get_shunts.c_str(), py::return_value_policy::reference_internal)
        .def("get_storages", &LSGrid::get_storages, DocLSGrid::get_storages.c_str(), py::return_value_policy::reference_internal)
        .def("get_loads", &LSGrid::get_loads, DocLSGrid::get_loads.c_str(), py::return_value_policy::reference_internal)

        // pypowsybl compat names
        .def("get_voltage_levels", &LSGrid::get_substations, DocLSGrid::get_substations.c_str(), py::return_value_policy::reference_internal)
        .def("get_2_windings_transformers", &LSGrid::get_trafos, DocLSGrid::get_trafos.c_str(), py::return_value_policy::reference_internal)
        .def("get_shunt_compensators", &LSGrid::get_shunts, DocLSGrid::get_shunts.c_str(), py::return_value_policy::reference_internal)

        // modify the grid
        .def("turnedoff_no_pv", &LSGrid::turnedoff_no_pv, DocLSGrid::turnedoff_no_pv.c_str())
        .def("turnedoff_pv", &LSGrid::turnedoff_pv, DocLSGrid::turnedoff_pv.c_str())
        .def("get_turnedoff_gen_pv", &LSGrid::get_turnedoff_gen_pv, DocLSGrid::get_turnedoff_gen_pv.c_str())
        .def("update_slack_weights", &LSGrid::update_slack_weights, DocLSGrid::update_slack_weights.c_str())
        .def("update_slack_weights_by_id", &LSGrid::update_slack_weights_by_id, DocLSGrid::update_slack_weights_by_id.c_str())
        .def("assign_slack_to_most_connected", &LSGrid::assign_slack_to_most_connected, DocLSGrid::assign_slack_to_most_connected.c_str())
        .def("set_reference_slack_bus", &LSGrid::set_reference_slack_bus, DocLSGrid::set_reference_slack_bus.c_str())
        .def("get_reference_slack_bus", &LSGrid::get_reference_slack_bus, DocLSGrid::get_reference_slack_bus.c_str())
        .def("consider_only_main_component", &LSGrid::consider_only_main_component, DocLSGrid::consider_only_main_component.c_str())
        .def("set_ignore_status_global", &LSGrid::set_ignore_status_global, DocLSGrid::set_ignore_status_global.c_str())
        .def("set_synch_status_both_side", &LSGrid::set_synch_status_both_side, DocLSGrid::set_synch_status_both_side.c_str())
        .def("get_ignore_status_global", &LSGrid::get_ignore_status_global, DocLSGrid::get_ignore_status_global.c_str())
        .def("get_synch_status_both_side", &LSGrid::get_synch_status_both_side, DocLSGrid::get_synch_status_both_side.c_str())

        // names
        .def("set_line_names", &LSGrid::set_line_names, DocLSGrid::set_line_names.c_str())
        .def("get_line_names", &LSGrid::get_line_names, DocLSGrid::get_line_names.c_str())
        .def("set_dcline_names", &LSGrid::set_dcline_names, DocLSGrid::set_dcline_names.c_str())
        .def("set_trafo_names", &LSGrid::set_trafo_names, DocLSGrid::set_trafo_names.c_str())
        .def("get_trafo_names", &LSGrid::get_trafo_names, DocLSGrid::get_trafo_names.c_str())
        .def("set_line_current_limit_side1", &LSGrid::set_line_current_limit_side1, DocLSGrid::set_line_current_limit_side1.c_str())
        .def("set_line_current_limit_side2", &LSGrid::set_line_current_limit_side2, DocLSGrid::set_line_current_limit_side2.c_str())
        .def("set_trafo_current_limit_side1", &LSGrid::set_trafo_current_limit_side1, DocLSGrid::set_trafo_current_limit_side1.c_str())
        .def("set_trafo_current_limit_side2", &LSGrid::set_trafo_current_limit_side2, DocLSGrid::set_trafo_current_limit_side2.c_str())
        .def("set_gen_names", &LSGrid::set_gen_names, DocLSGrid::set_gen_names.c_str())
        .def("set_load_names", &LSGrid::set_load_names, DocLSGrid::set_load_names.c_str())
        .def("set_storage_names", &LSGrid::set_storage_names, DocLSGrid::set_storage_names.c_str())
        .def("set_sgen_names", &LSGrid::set_sgen_names, DocLSGrid::set_sgen_names.c_str())
        .def("set_shunt_names", &LSGrid::set_shunt_names, DocLSGrid::set_shunt_names.c_str())
        .def("set_substation_names", &LSGrid::set_substation_names, DocLSGrid::set_substation_names.c_str())
        .def("get_substation_names", &LSGrid::get_substation_names, DocLSGrid::get_substation_names.c_str())

        .def("deactivate_bus", &LSGrid::deactivate_bus_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("reactivate_bus", &LSGrid::reactivate_bus_python, DocLSGrid::_internal_do_not_use.c_str())

        .def("deactivate_powerline", &LSGrid::deactivate_powerline, DocLSGrid::deactivate_powerline.c_str())
        .def("reactivate_powerline", &LSGrid::reactivate_powerline, DocLSGrid::reactivate_powerline.c_str())
        .def("deactivate_powerline_side1", &LSGrid::deactivate_powerline_side1, DocLSGrid::deactivate_powerline_side1.c_str())
        .def("deactivate_powerline_side2", &LSGrid::deactivate_powerline_side2, DocLSGrid::deactivate_powerline_side2.c_str())
        .def("reactivate_powerline_side1", &LSGrid::reactivate_powerline_side1, DocLSGrid::reactivate_powerline_side1.c_str())
        .def("reactivate_powerline_side2", &LSGrid::reactivate_powerline_side2, DocLSGrid::reactivate_powerline_side2.c_str())
        .def("change_bus1_powerline", &LSGrid::change_bus1_powerline_python, DocLSGrid::change_bus1_powerline.c_str())
        .def("change_bus2_powerline", &LSGrid::change_bus2_powerline_python, DocLSGrid::change_bus2_powerline.c_str())
        .def("get_bus1_powerline", &LSGrid::get_bus1_powerline, DocLSGrid::get_bus1_powerline.c_str(), py::return_value_policy::reference)
        .def("get_bus2_powerline", &LSGrid::get_bus2_powerline, DocLSGrid::get_bus2_powerline.c_str(), py::return_value_policy::reference)

        .def("deactivate_trafo", &LSGrid::deactivate_trafo, DocLSGrid::deactivate_trafo.c_str())
        .def("reactivate_trafo", &LSGrid::reactivate_trafo, DocLSGrid::reactivate_trafo.c_str())
        .def("deactivate_trafo_side1", &LSGrid::deactivate_trafo_side1, DocLSGrid::deactivate_trafo_side1.c_str())
        .def("deactivate_trafo_side2", &LSGrid::deactivate_trafo_side2, DocLSGrid::deactivate_trafo_side2.c_str())
        .def("reactivate_trafo_side1", &LSGrid::reactivate_trafo_side1, DocLSGrid::reactivate_trafo_side1.c_str())
        .def("reactivate_trafo_side2", &LSGrid::reactivate_trafo_side2, DocLSGrid::reactivate_trafo_side2.c_str())
        .def("change_bus1_trafo", &LSGrid::change_bus1_trafo_python, DocLSGrid::change_bus1_trafo.c_str())
        .def("change_bus2_trafo", &LSGrid::change_bus2_trafo_python, DocLSGrid::change_bus2_trafo.c_str())
        .def("get_bus1_trafo", &LSGrid::get_bus1_trafo, DocLSGrid::get_bus1_trafo.c_str(), py::return_value_policy::reference)
        .def("get_bus2_trafo", &LSGrid::get_bus2_trafo, DocLSGrid::get_bus2_trafo.c_str(), py::return_value_policy::reference)
        .def("change_ratio_trafo", &LSGrid::change_ratio_trafo, DocLSGrid::change_ratio_trafo.c_str())
        .def("change_shift_trafo", &LSGrid::change_shift_trafo, DocLSGrid::change_shift_trafo.c_str())
        .def("change_shift_trafo_deg", &LSGrid::change_shift_trafo_deg, DocLSGrid::change_shift_trafo_deg.c_str())
        .def("set_trafo_shift_dependent_rx", &LSGrid::set_trafo_shift_dependent_rx,
            py::arg("enable"), py::arg("alpha_rad"), py::arg("rx_corr_pct"),
            DocLSGrid::set_trafo_shift_dependent_rx.c_str())
        .def("deactivate_load", &LSGrid::deactivate_load, DocLSGrid::deactivate_load.c_str())
        .def("reactivate_load", &LSGrid::reactivate_load, DocLSGrid::reactivate_load.c_str())
        .def("change_bus_load", &LSGrid::change_bus_load_python, DocLSGrid::change_bus_load.c_str())
        .def("get_bus_load", &LSGrid::get_bus_load, DocLSGrid::get_bus_load.c_str(), py::return_value_policy::reference)
        .def("change_p_load", &LSGrid::change_p_load, DocLSGrid::change_p_load.c_str())
        .def("change_q_load", &LSGrid::change_q_load, DocLSGrid::change_q_load.c_str())

        .def("deactivate_gen", &LSGrid::deactivate_gen, DocLSGrid::deactivate_gen.c_str())
        .def("reactivate_gen", &LSGrid::reactivate_gen, DocLSGrid::reactivate_gen.c_str())
        .def("change_bus_gen", &LSGrid::change_bus_gen_python, DocLSGrid::change_bus_gen.c_str())
        .def("get_bus_gen", &LSGrid::get_bus_gen, DocLSGrid::get_bus_gen.c_str(), py::return_value_policy::reference)
        .def("change_p_gen", &LSGrid::change_p_gen, DocLSGrid::change_p_gen.c_str())
        .def("change_v_gen", &LSGrid::change_v_gen, DocLSGrid::change_v_gen.c_str())
        .def("set_gen_regulated_bus", &LSGrid::set_gen_regulated_bus, DocLSGrid::set_gen_regulated_bus.c_str())
        .def("deactivate_svc", &LSGrid::deactivate_svc, DocLSGrid::deactivate_svc.c_str())
        .def("reactivate_svc", &LSGrid::reactivate_svc, DocLSGrid::reactivate_svc.c_str())
        .def("change_bus_svc", &LSGrid::change_bus_svc_python, DocLSGrid::change_bus_svc.c_str())
        .def("get_bus_svc", &LSGrid::get_bus_svc, DocLSGrid::get_bus_svc.c_str())
        .def("set_svc_names", &LSGrid::set_svc_names, DocLSGrid::set_svc_names.c_str())

        .def("deactivate_shunt", &LSGrid::deactivate_shunt, DocLSGrid::deactivate_shunt.c_str())
        .def("reactivate_shunt", &LSGrid::reactivate_shunt, DocLSGrid::reactivate_shunt.c_str())
        .def("change_bus_shunt", &LSGrid::change_bus_shunt_python, DocLSGrid::change_bus_shunt.c_str())
        .def("get_bus_shunt", &LSGrid::get_bus_shunt, DocLSGrid::get_bus_shunt.c_str(), py::return_value_policy::reference)
        .def("change_p_shunt", &LSGrid::change_p_shunt, DocLSGrid::change_p_shunt.c_str())
        .def("change_q_shunt", &LSGrid::change_q_shunt, DocLSGrid::change_q_shunt.c_str())

        .def("deactivate_sgen", &LSGrid::deactivate_sgen, DocLSGrid::deactivate_sgen.c_str())
        .def("reactivate_sgen", &LSGrid::reactivate_sgen, DocLSGrid::reactivate_sgen.c_str())
        .def("change_bus_sgen", &LSGrid::change_bus_sgen_python, DocLSGrid::change_bus_sgen.c_str())
        .def("get_bus_sgen", &LSGrid::get_bus_sgen, DocLSGrid::get_bus_sgen.c_str(), py::return_value_policy::reference)
        .def("change_p_sgen", &LSGrid::change_p_sgen, DocLSGrid::change_p_sgen.c_str())
        .def("change_q_sgen", &LSGrid::change_q_sgen, DocLSGrid::change_q_sgen.c_str())

        .def("deactivate_storage", &LSGrid::deactivate_storage, DocLSGrid::deactivate_storage.c_str())
        .def("reactivate_storage", &LSGrid::reactivate_storage, DocLSGrid::reactivate_storage.c_str())
        .def("change_bus_storage", &LSGrid::change_bus_storage_python, DocLSGrid::change_bus_storage.c_str())
        .def("get_bus_storage", &LSGrid::get_bus_storage, DocLSGrid::get_bus_storage.c_str(), py::return_value_policy::reference)
        .def("change_p_storage", &LSGrid::change_p_storage, DocLSGrid::change_p_storage.c_str())
        .def("change_q_storage", &LSGrid::change_q_storage, DocLSGrid::change_q_storage.c_str())

        .def("deactivate_dcline", &LSGrid::deactivate_dcline, DocLSGrid::deactivate_dcline.c_str())
        .def("deactivate_dcline_side1", &LSGrid::deactivate_dcline_side1, DocLSGrid::deactivate_dcline_side1.c_str())
        .def("deactivate_dcline_side2", &LSGrid::deactivate_dcline_side2, DocLSGrid::deactivate_dcline_side2.c_str())
        .def("reactivate_dcline", &LSGrid::reactivate_dcline, DocLSGrid::reactivate_dcline.c_str())
        .def("change_p_dcline", &LSGrid::change_p_dcline, DocLSGrid::change_p_dcline.c_str())
        .def("change_v1_dcline", &LSGrid::change_v1_dcline, DocLSGrid::change_v1_dcline.c_str())
        .def("change_v2_dcline", &LSGrid::change_v2_dcline, DocLSGrid::change_v2_dcline.c_str())
        .def("change_bus1_dcline", &LSGrid::change_bus1_dcline, DocLSGrid::change_bus1_dcline.c_str())
        .def("change_bus2_dcline", &LSGrid::change_bus2_dcline, DocLSGrid::change_bus2_dcline.c_str())
        .def("get_bus1_dcline", &LSGrid::get_bus1_dcline, DocLSGrid::get_bus1_dcline.c_str(), py::return_value_policy::reference)
        .def("get_bus2_dcline", &LSGrid::get_bus2_dcline, DocLSGrid::get_bus2_dcline.c_str(), py::return_value_policy::reference)
        .def("set_status_droop_hvdc", &LSGrid::set_status_droop_hvdc, DocLSGrid::set_status_droop_hvdc.c_str())
        .def("get_status_droop_hvdc", &LSGrid::get_status_droop_hvdc, DocLSGrid::get_status_droop_hvdc.c_str())
        .def("get_status_droop_hvdc_vect", &LSGrid::get_status_droop_hvdc_vect, DocLSGrid::get_status_droop_hvdc_vect.c_str())

        // get back the results
        .def("get_V", &LSGrid::get_V, DocLSGrid::get_V.c_str())
        .def("get_Va", &LSGrid::get_Va, DocLSGrid::get_Va.c_str())
        .def("get_Vm", &LSGrid::get_Vm, DocLSGrid::get_Vm.c_str())
        .def("get_V_solver", &LSGrid::get_V_solver, DocLSGrid::get_V_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_Va_solver", &LSGrid::get_Va_solver, DocLSGrid::get_Va_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_Vm_solver", &LSGrid::get_Vm_solver, DocLSGrid::get_Vm_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_J_solver", &LSGrid::get_J_python_solver, DocLSGrid::get_J_python_solver.c_str(), py::return_value_policy::reference)

        .def("id_me_to_ac_solver", &LSGrid::id_me_to_ac_solver_numpy, DocLSGrid::id_me_to_ac_solver.c_str(), py::return_value_policy::reference)
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
        // per-family variants: the AC and the DC solver each keep their own split
        .def("get_ac_pv_solver", &LSGrid::get_ac_pv_solver_numpy, DocLSGrid::get_ac_pv_solver.c_str(), py::return_value_policy::reference)
        .def("get_dc_pv_solver", &LSGrid::get_dc_pv_solver_numpy, DocLSGrid::get_dc_pv_solver.c_str(), py::return_value_policy::reference)
        .def("get_ac_pq_solver", &LSGrid::get_ac_pq_solver_numpy, DocLSGrid::get_ac_pq_solver.c_str(), py::return_value_policy::reference)
        .def("get_dc_pq_solver", &LSGrid::get_dc_pq_solver_numpy, DocLSGrid::get_dc_pq_solver.c_str(), py::return_value_policy::reference)
        .def("get_slack_ids_solver", &LSGrid::get_slack_ids_solver_numpy, DocLSGrid::get_slack_ids_solver.c_str(), py::return_value_policy::reference)
        .def("get_slack_ids_dc_solver", &LSGrid::get_slack_ids_dc_solver_numpy, DocLSGrid::get_slack_ids_dc_solver.c_str(), py::return_value_policy::reference)
        .def("get_slack_weights_solver", &LSGrid::get_slack_weights_solver, DocLSGrid::get_slack_weights_solver.c_str(), py::return_value_policy::reference)
        .def("get_ac_slack_weights_solver", &LSGrid::get_ac_slack_weights_solver, DocLSGrid::get_ac_slack_weights_solver.c_str(), py::return_value_policy::reference)
        .def("get_dc_slack_weights_solver", &LSGrid::get_dc_slack_weights_solver, DocLSGrid::get_dc_slack_weights_solver.c_str(), py::return_value_policy::reference)
        .def("get_slack_col_solver", &LSGrid::get_slack_col_solver, DocLSGrid::get_slack_col_solver.c_str())
        .def("get_slack_absorbed_solver", &LSGrid::get_slack_absorbed_solver, DocLSGrid::get_slack_absorbed_solver.c_str())
        .def("get_controller_q_solver", &LSGrid::get_controller_q_solver, py::return_value_policy::reference, DocLSGrid::get_controller_q_solver.c_str())
        .def("get_controller_kind_solver", &LSGrid::get_controller_kind_solver, py::return_value_policy::reference, DocLSGrid::get_controller_kind_solver.c_str())
        .def("get_controller_elem_id_solver", &LSGrid::get_controller_elem_id_solver, py::return_value_policy::reference, DocLSGrid::get_controller_elem_id_solver.c_str())
        .def("get_controller_q_col_solver", &LSGrid::get_controller_q_col_solver, py::return_value_policy::reference, DocLSGrid::get_controller_q_col_solver.c_str())

        .def("get_p_buses_solver", &LSGrid::get_p_buses_solver, py::return_value_policy::reference, DocLSGrid::get_p_buses_solver.c_str())
        .def("get_p_rows_solver", &LSGrid::get_p_rows_solver, py::return_value_policy::reference, DocLSGrid::get_p_rows_solver.c_str())
        .def("get_q_buses_solver", &LSGrid::get_q_buses_solver, py::return_value_policy::reference, DocLSGrid::get_q_buses_solver.c_str())
        .def("get_q_rows_solver", &LSGrid::get_q_rows_solver, py::return_value_policy::reference, DocLSGrid::get_q_rows_solver.c_str())
        .def("get_theta_buses_solver", &LSGrid::get_theta_buses_solver, py::return_value_policy::reference, DocLSGrid::get_theta_buses_solver.c_str())
        .def("get_theta_cols_solver", &LSGrid::get_theta_cols_solver, py::return_value_policy::reference, DocLSGrid::get_theta_cols_solver.c_str())
        .def("get_vm_buses_solver", &LSGrid::get_vm_buses_solver, py::return_value_policy::reference, DocLSGrid::get_vm_buses_solver.c_str())
        .def("get_vm_cols_solver", &LSGrid::get_vm_cols_solver, py::return_value_policy::reference, DocLSGrid::get_vm_cols_solver.c_str())
        .def("get_hvdc_droop_data_solver", &LSGrid::get_hvdc_droop_data_solver, DocLSGrid::get_hvdc_droop_data_solver.c_str())

        .def("get_Ybus", &LSGrid::get_Ybus, DocLSGrid::get_Ybus.c_str())
        .def("get_dcYbus", &LSGrid::get_dcYbus, DocLSGrid::get_dcYbus.c_str())
        .def("get_Sbus", &LSGrid::get_Sbus, DocLSGrid::get_Sbus.c_str())
        .def("get_dcSbus", &LSGrid::get_dcSbus, DocLSGrid::get_dcSbus.c_str())
        .def("get_Ybus_solver", &LSGrid::get_Ybus_solver, DocLSGrid::get_Ybus_solver.c_str(), py::return_value_policy::reference)
        .def("get_dcYbus_solver", &LSGrid::get_dcYbus_solver, DocLSGrid::get_dcYbus_solver.c_str(), py::return_value_policy::reference)
        .def("get_Sbus_solver", &LSGrid::get_Sbus_solver, DocLSGrid::get_Sbus_solver.c_str(), py::return_value_policy::reference)
        .def("get_dcSbus_solver", &LSGrid::get_dcSbus_solver, DocLSGrid::get_dcSbus_solver.c_str(), py::return_value_policy::reference)

        .def("check_solution", &LSGrid::check_solution, DocLSGrid::check_solution.c_str())

        .def("get_loads_res", &LSGrid::get_loads_res, DocLSGrid::get_loads_res.c_str(), py::return_value_policy::reference)
        .def("get_loads_status", &LSGrid::get_loads_status, DocLSGrid::get_loads_status.c_str(), py::return_value_policy::reference)
        .def("get_shunts_res", &LSGrid::get_shunts_res, DocLSGrid::get_shunts_res.c_str(), py::return_value_policy::reference)
        .def("get_shunts_status", &LSGrid::get_shunts_status, DocLSGrid::get_shunts_status.c_str(), py::return_value_policy::reference)
        .def("get_gen_res", &LSGrid::get_gen_res, DocLSGrid::get_gen_res.c_str(), py::return_value_policy::reference)
        .def("get_gen_status", &LSGrid::get_gen_status, DocLSGrid::get_gen_status.c_str(), py::return_value_policy::reference)
        .def("get_line_res1", &LSGrid::get_line_res1, DocLSGrid::get_line_res1.c_str(), py::return_value_policy::reference)
        .def("get_line_res2", &LSGrid::get_line_res2, DocLSGrid::get_line_res2.c_str(), py::return_value_policy::reference)
        .def("get_lines_status", &LSGrid::get_lines_status, DocLSGrid::get_lines_status.c_str(), py::return_value_policy::reference)
        .def("get_lines_status_side1", &LSGrid::get_lines_status_side1, DocLSGrid::get_lines_status_side1.c_str(), py::return_value_policy::reference)
        .def("get_lines_status_side2", &LSGrid::get_lines_status_side2, DocLSGrid::get_lines_status_side2.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res1", &LSGrid::get_trafo_res1, DocLSGrid::get_trafo_res1.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res2", &LSGrid::get_trafo_res2, DocLSGrid::get_trafo_res2.c_str(), py::return_value_policy::reference)
        .def("get_trafo_status", &LSGrid::get_trafo_status, DocLSGrid::get_trafo_status.c_str(), py::return_value_policy::reference)
        .def("get_trafo_status_side1", &LSGrid::get_trafo_status_side1, DocLSGrid::get_trafo_status_side1.c_str(), py::return_value_policy::reference)
        .def("get_trafo_status_side2", &LSGrid::get_trafo_status_side2, DocLSGrid::get_trafo_status_side2.c_str(), py::return_value_policy::reference)
        .def("get_storages_res", &LSGrid::get_storages_res, DocLSGrid::get_storages_res.c_str(), py::return_value_policy::reference)
        .def("get_storages_status", &LSGrid::get_storages_status, DocLSGrid::get_storages_status.c_str(), py::return_value_policy::reference)
        .def("get_sgens_res", &LSGrid::get_sgens_res, DocLSGrid::get_sgens_res.c_str(), py::return_value_policy::reference)
        .def("get_sgens_status", &LSGrid::get_sgens_status, DocLSGrid::get_sgens_status.c_str(), py::return_value_policy::reference)

        .def("get_gen_theta", &LSGrid::get_gen_theta, DocLSGrid::get_gen_theta.c_str(), py::return_value_policy::reference)
        .def("get_load_theta", &LSGrid::get_load_theta, DocLSGrid::get_load_theta.c_str(), py::return_value_policy::reference)
        .def("get_shunt_theta", &LSGrid::get_shunt_theta, DocLSGrid::get_shunt_theta.c_str(), py::return_value_policy::reference)
        .def("get_storage_theta", &LSGrid::get_storage_theta, DocLSGrid::get_storage_theta.c_str(), py::return_value_policy::reference)
        .def("get_line_theta1", &LSGrid::get_line_theta1, DocLSGrid::get_line_theta1.c_str(), py::return_value_policy::reference)
        .def("get_line_theta2", &LSGrid::get_line_theta2, DocLSGrid::get_line_theta2.c_str(), py::return_value_policy::reference)
        .def("get_trafo_theta1", &LSGrid::get_trafo_theta1, DocLSGrid::get_trafo_theta1.c_str(), py::return_value_policy::reference)
        .def("get_trafo_theta2", &LSGrid::get_trafo_theta2, DocLSGrid::get_trafo_theta2.c_str(), py::return_value_policy::reference)

        .def("get_all_shunt_buses", &LSGrid::get_all_shunt_buses_numpy, DocLSGrid::get_all_shunt_buses.c_str(), py::return_value_policy::reference)
        .def("get_loads_res_full", &LSGrid::get_loads_res_full, DocLSGrid::get_loads_res_full.c_str(), py::return_value_policy::reference)
        .def("get_shunts_res_full", &LSGrid::get_shunts_res_full, DocLSGrid::get_shunts_res_full.c_str(), py::return_value_policy::reference)
        .def("get_gen_res_full", &LSGrid::get_gen_res_full, DocLSGrid::get_gen_res_full.c_str(), py::return_value_policy::reference)
        .def("get_line_res1_full", &LSGrid::get_line_res1_full, DocLSGrid::get_line_res1_full.c_str(), py::return_value_policy::reference)
        .def("get_line_res2_full", &LSGrid::get_line_res2_full, DocLSGrid::get_line_res2_full.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res1_full", &LSGrid::get_trafo_res1_full, DocLSGrid::get_trafo_res1_full.c_str(), py::return_value_policy::reference)
        .def("get_trafo_res2_full", &LSGrid::get_trafo_res2_full, DocLSGrid::get_trafo_res2_full.c_str(), py::return_value_policy::reference)
        .def("get_storages_res_full", &LSGrid::get_storages_res_full, DocLSGrid::get_storages_res_full.c_str(), py::return_value_policy::reference)
        .def("get_sgens_res_full", &LSGrid::get_sgens_res_full, DocLSGrid::get_sgens_res_full.c_str(), py::return_value_policy::reference)
        .def("get_dcline_res1_full", &LSGrid::get_dcline_res1_full, DocLSGrid::get_dcline_res1_full.c_str(), py::return_value_policy::reference)
        .def("get_dcline_res2_full", &LSGrid::get_dcline_res2_full, DocLSGrid::get_dcline_res2_full.c_str(), py::return_value_policy::reference)

        .def("get_shunt_target_p", &LSGrid::get_shunt_target_p, DocLSGrid::get_shunt_target_p.c_str(), py::return_value_policy::reference)
        .def("get_load_target_p", &LSGrid::get_load_target_p, DocLSGrid::get_load_target_p.c_str(), py::return_value_policy::reference)
        .def("get_gen_target_p", &LSGrid::get_gen_target_p, DocLSGrid::get_gen_target_p.c_str(), py::return_value_policy::reference)
        .def("get_sgen_target_p", &LSGrid::get_sgen_target_p, DocLSGrid::get_sgen_target_p.c_str(), py::return_value_policy::reference)
        .def("get_storage_target_p", &LSGrid::get_storage_target_p, DocLSGrid::get_storage_target_p.c_str(), py::return_value_policy::reference)

        // do something with the grid
        .def("deactivate_result_computation", &LSGrid::deactivate_result_computation, DocLSGrid::deactivate_result_computation.c_str())
        .def("reactivate_result_computation", &LSGrid::reactivate_result_computation, DocLSGrid::reactivate_result_computation.c_str())
        .def("dc_pf", &LSGrid::dc_pf, DocLSGrid::dc_pf.c_str())
        .def("ac_pf", &LSGrid::ac_pf, DocLSGrid::ac_pf.c_str())
        .def("unset_changes", &LSGrid::unset_changes, DocLSGrid::unset_changes.c_str())
        // cache reuse (on by default, per solver family)
        .def("allow_ac_cache_reuse", &LSGrid::allow_ac_cache_reuse, py::arg("allowed"), DocLSGrid::allow_ac_cache_reuse.c_str())
        .def("allow_dc_cache_reuse", &LSGrid::allow_dc_cache_reuse, py::arg("allowed"), DocLSGrid::allow_dc_cache_reuse.c_str())
        .def("allow_cache_reuse", &LSGrid::allow_cache_reuse, py::arg("allowed"), DocLSGrid::allow_cache_reuse.c_str())
        .def("get_allow_ac_cache_reuse", &LSGrid::get_allow_ac_cache_reuse, DocLSGrid::get_allow_ac_cache_reuse.c_str())
        .def("get_allow_dc_cache_reuse", &LSGrid::get_allow_dc_cache_reuse, DocLSGrid::get_allow_dc_cache_reuse.c_str())
        .def("get_allow_cache_reuse", &LSGrid::get_allow_cache_reuse, DocLSGrid::get_allow_cache_reuse.c_str())
        .def("prevent_ac_cache_reuse", &LSGrid::prevent_ac_cache_reuse, DocLSGrid::prevent_ac_cache_reuse.c_str())
        .def("prevent_dc_cache_reuse", &LSGrid::prevent_dc_cache_reuse, DocLSGrid::prevent_dc_cache_reuse.c_str())
        .def("prevent_cache_reuse", &LSGrid::prevent_cache_reuse, DocLSGrid::prevent_cache_reuse.c_str())
        .def("tell_recompute_ybus", &LSGrid::tell_recompute_ybus, DocLSGrid::_internal_do_not_use.c_str())
        .def("tell_recompute_sbus", &LSGrid::tell_recompute_sbus, DocLSGrid::_internal_do_not_use.c_str())
        .def("tell_solver_need_reset", &LSGrid::tell_solver_need_reset, DocLSGrid::tell_solver_need_reset.c_str())
        .def("tell_ybus_change_sparsity_pattern", &LSGrid::tell_ybus_change_sparsity_pattern, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_ac_algo_controler", &LSGrid::get_ac_algo_controler, DocLSGrid::get_ac_algo_controler.c_str(), py::return_value_policy::reference)
        .def("get_dc_algo_controler", &LSGrid::get_dc_algo_controler, DocLSGrid::get_dc_algo_controler.c_str(), py::return_value_policy::reference)
        // .def("get_solver_control",  &LSGrid::get_algo_controler, "DEPRECATED use 'get_algo_controler'", py::return_value_policy::reference)
        .def("compute_newton", &LSGrid::ac_pf, DocLSGrid::ac_pf.c_str())
        // get_ptdf/get_ptdf_solver/get_lodf/get_Bf/get_Bf_solver all return their
        // matrix BY VALUE (freshly computed, not a reference to persistent state):
        // no return_value_policy::reference here, since that would wrap the numpy
        // array around the returned temporary's memory, which is freed as soon as
        // this call returns (dangling on the Python side). The default policy
        // (copy/move into a Python-owned array) is the only safe choice for a
        // by-value return.
        .def("get_ptdf", &LSGrid::get_ptdf, DocLSGrid::get_ptdf.c_str())
        .def("get_ptdf_solver", &LSGrid::get_ptdf_solver, DocLSGrid::get_ptdf_solver.c_str())
        .def("get_lodf", &LSGrid::get_lodf, DocLSGrid::get_lodf.c_str())
        .def("get_Bf", &LSGrid::get_Bf, DocLSGrid::get_Bf.c_str())
        .def("get_Bf_solver", &LSGrid::get_Bf_solver, DocLSGrid::get_Bf_solver.c_str())

        // apply action faster (optimized for grid2op representation)
        .def("update_gens_p", &LSGrid::update_gens_p, DocLSGrid::update_gens_p.c_str())
        .def("update_sgens_p", &LSGrid::update_sgens_p, DocLSGrid::update_sgens_p.c_str())
        .def("update_gens_v", &LSGrid::update_gens_v, DocLSGrid::update_gens_v.c_str())
        .def("update_loads_p", &LSGrid::update_loads_p, DocLSGrid::update_loads_p.c_str())
        .def("update_loads_q", &LSGrid::update_loads_q, DocLSGrid::update_loads_q.c_str())
        .def("update_topo", &LSGrid::update_topo, DocLSGrid::update_topo.c_str())
        .def("update_storages_p", &LSGrid::update_storages_p, DocLSGrid::update_storages_p.c_str())

        // auxiliary functions
        .def("set_n_sub", &LSGrid::set_n_sub, DocLSGrid::set_n_sub.c_str())
        .def("get_n_sub", &LSGrid::get_n_sub, DocLSGrid::get_n_sub.c_str())
        .def("set_max_nb_bus_per_sub", &LSGrid::set_max_nb_bus_per_sub, DocLSGrid::set_max_nb_bus_per_sub.c_str())
        .def("set_load_pos_topo_vect", &LSGrid::set_load_pos_topo_vect, DocLSGrid::set_load_pos_topo_vect.c_str())
        .def("set_gen_pos_topo_vect", &LSGrid::set_gen_pos_topo_vect, DocLSGrid::set_gen_pos_topo_vect.c_str())
        .def("set_line_pos1_topo_vect", &LSGrid::set_line_pos1_topo_vect, DocLSGrid::set_line_pos1_topo_vect.c_str())
        .def("set_line_pos2_topo_vect", &LSGrid::set_line_pos2_topo_vect, DocLSGrid::set_line_pos2_topo_vect.c_str())
        .def("set_trafo_pos1_topo_vect", &LSGrid::set_trafo_pos1_topo_vect, DocLSGrid::set_trafo_pos1_topo_vect.c_str())
        .def("set_trafo_pos2_topo_vect", &LSGrid::set_trafo_pos2_topo_vect, DocLSGrid::set_trafo_pos2_topo_vect.c_str())
        .def("set_storage_pos_topo_vect", &LSGrid::set_storage_pos_topo_vect, DocLSGrid::set_storage_pos_topo_vect.c_str())
        .def("set_load_to_subid", &LSGrid::set_load_to_subid, DocLSGrid::set_load_to_subid.c_str())
        .def("set_gen_to_subid", &LSGrid::set_gen_to_subid, DocLSGrid::set_gen_to_subid.c_str())
        .def("set_shunt_to_subid", &LSGrid::set_shunt_to_subid, DocLSGrid::set_shunt_to_subid.c_str())
        .def("set_line_to_sub1_id", &LSGrid::set_line_to_sub1_id, DocLSGrid::set_line_to_sub1_id.c_str())
        .def("set_line_to_sub2_id", &LSGrid::set_line_to_sub2_id, DocLSGrid::set_line_to_sub2_id.c_str())
        .def("set_trafo_to_sub1_id", &LSGrid::set_trafo_to_sub1_id, DocLSGrid::set_trafo_to_sub1_id.c_str())
        .def("set_trafo_to_sub2_id", &LSGrid::set_trafo_to_sub2_id, DocLSGrid::set_trafo_to_sub2_id.c_str())
        .def("set_storage_to_subid", &LSGrid::set_storage_to_subid, DocLSGrid::set_storage_to_subid.c_str())

        // debug functions (might disappear without further notice)
        .def("debug_get_Bp_python", &LSGrid::debug_get_Bp_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("debug_get_Bpp_python", &LSGrid::debug_get_Bpp_python, DocLSGrid::_internal_do_not_use.c_str());
}
