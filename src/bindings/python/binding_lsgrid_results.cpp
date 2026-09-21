// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


// The results of LSGrid (per bus, per element, per solver family) -- one of the four translation units the LSGrid bindings are spread
// over, see the note at the top of binding_lsgrid.cpp.

#include "binding_declarations.hpp"
#include "LSGrid.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

void bind_lsgrid_results(py::class_<LSGrid> & cls) {
    cls
        // get back the results
        .def("get_V", &LSGrid::get_V, DocLSGrid::get_V.c_str())
        .def("get_Va", &LSGrid::get_Va, DocLSGrid::get_Va.c_str())
        .def("get_Vm", &LSGrid::get_Vm, DocLSGrid::get_Vm.c_str())
        .def("get_V_solver", &LSGrid::get_V_solver, DocLSGrid::get_V_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_Va_solver", &LSGrid::get_Va_solver, DocLSGrid::get_Va_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_Vm_solver", &LSGrid::get_Vm_solver, DocLSGrid::get_Vm_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_J_solver", &LSGrid::get_J_python_solver, DocLSGrid::get_J_python_solver.c_str(), py::return_value_policy::reference_internal)

        .def("id_me_to_ac_solver", &LSGrid::id_me_to_ac_solver_numpy, DocLSGrid::id_me_to_ac_solver.c_str(), py::return_value_policy::reference_internal)
        .def("id_ac_solver_to_me", &LSGrid::id_ac_solver_to_me_numpy, DocLSGrid::id_ac_solver_to_me.c_str(), py::return_value_policy::reference_internal)
        .def("id_me_to_dc_solver", &LSGrid::id_me_to_dc_solver_numpy, DocLSGrid::id_me_to_dc_solver.c_str(), py::return_value_policy::reference_internal)
        .def("id_dc_solver_to_me", &LSGrid::id_dc_solver_to_me_numpy, DocLSGrid::id_dc_solver_to_me.c_str(), py::return_value_policy::reference_internal)
        .def("total_bus", &LSGrid::total_bus, DocLSGrid::total_bus.c_str())
        .def("nb_connected_bus", &LSGrid::nb_connected_bus, DocLSGrid::nb_connected_bus.c_str())

        .def("get_pv", &LSGrid::get_pv_numpy, DocLSGrid::get_pv.c_str(), py::return_value_policy::reference_internal)
        .def("get_pq", &LSGrid::get_pq_numpy, DocLSGrid::get_pq.c_str(), py::return_value_policy::reference_internal)
        .def("get_slack_ids", &LSGrid::get_slack_ids_numpy, DocLSGrid::get_slack_ids.c_str(), py::return_value_policy::reference_internal)
        .def("get_slack_ids_dc", &LSGrid::get_slack_ids_dc_numpy, DocLSGrid::get_slack_ids_dc.c_str(), py::return_value_policy::reference_internal)
        .def("get_slack_weights", &LSGrid::get_slack_weights, DocLSGrid::get_slack_weights.c_str(), py::return_value_policy::reference_internal)
        .def("get_pv_solver", &LSGrid::get_pv_solver_numpy, DocLSGrid::get_pv_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_pq_solver", &LSGrid::get_pq_solver_numpy, DocLSGrid::get_pq_solver.c_str(), py::return_value_policy::reference_internal)
        // per-family variants: the AC and the DC solver each keep their own split
        .def("get_ac_pv_solver", &LSGrid::get_ac_pv_solver_numpy, DocLSGrid::get_ac_pv_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_dc_pv_solver", &LSGrid::get_dc_pv_solver_numpy, DocLSGrid::get_dc_pv_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_ac_pq_solver", &LSGrid::get_ac_pq_solver_numpy, DocLSGrid::get_ac_pq_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_dc_pq_solver", &LSGrid::get_dc_pq_solver_numpy, DocLSGrid::get_dc_pq_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_slack_ids_solver", &LSGrid::get_slack_ids_solver_numpy, DocLSGrid::get_slack_ids_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_slack_ids_dc_solver", &LSGrid::get_slack_ids_dc_solver_numpy, DocLSGrid::get_slack_ids_dc_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_slack_weights_solver", &LSGrid::get_slack_weights_solver, DocLSGrid::get_slack_weights_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_ac_slack_weights_solver", &LSGrid::get_ac_slack_weights_solver, DocLSGrid::get_ac_slack_weights_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_dc_slack_weights_solver", &LSGrid::get_dc_slack_weights_solver, DocLSGrid::get_dc_slack_weights_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_slack_col_solver", &LSGrid::get_slack_col_solver, DocLSGrid::get_slack_col_solver.c_str())
        .def("get_slack_absorbed_solver", &LSGrid::get_slack_absorbed_solver, DocLSGrid::get_slack_absorbed_solver.c_str())
        .def("get_controller_q_solver", &LSGrid::get_controller_q_solver, py::return_value_policy::reference_internal, DocLSGrid::get_controller_q_solver.c_str())
        .def("get_controller_kind_solver", &LSGrid::get_controller_kind_solver, py::return_value_policy::reference_internal, DocLSGrid::get_controller_kind_solver.c_str())
        .def("get_controller_elem_id_solver", &LSGrid::get_controller_elem_id_solver, py::return_value_policy::reference_internal, DocLSGrid::get_controller_elem_id_solver.c_str())
        .def("get_controller_q_col_solver", &LSGrid::get_controller_q_col_solver, py::return_value_policy::reference_internal, DocLSGrid::get_controller_q_col_solver.c_str())

        .def("get_p_buses_solver", &LSGrid::get_p_buses_solver, py::return_value_policy::reference_internal, DocLSGrid::get_p_buses_solver.c_str())
        .def("get_p_rows_solver", &LSGrid::get_p_rows_solver, py::return_value_policy::reference_internal, DocLSGrid::get_p_rows_solver.c_str())
        .def("get_q_buses_solver", &LSGrid::get_q_buses_solver, py::return_value_policy::reference_internal, DocLSGrid::get_q_buses_solver.c_str())
        .def("get_q_rows_solver", &LSGrid::get_q_rows_solver, py::return_value_policy::reference_internal, DocLSGrid::get_q_rows_solver.c_str())
        .def("get_theta_buses_solver", &LSGrid::get_theta_buses_solver, py::return_value_policy::reference_internal, DocLSGrid::get_theta_buses_solver.c_str())
        .def("get_theta_cols_solver", &LSGrid::get_theta_cols_solver, py::return_value_policy::reference_internal, DocLSGrid::get_theta_cols_solver.c_str())
        .def("get_vm_buses_solver", &LSGrid::get_vm_buses_solver, py::return_value_policy::reference_internal, DocLSGrid::get_vm_buses_solver.c_str())
        .def("get_vm_cols_solver", &LSGrid::get_vm_cols_solver, py::return_value_policy::reference_internal, DocLSGrid::get_vm_cols_solver.c_str())
        .def("get_hvdc_droop_data_solver", &LSGrid::get_hvdc_droop_data_solver, DocLSGrid::get_hvdc_droop_data_solver.c_str())

        .def("get_Ybus", &LSGrid::get_Ybus, DocLSGrid::get_Ybus.c_str())
        .def("get_dcYbus", &LSGrid::get_dcYbus, DocLSGrid::get_dcYbus.c_str())
        .def("get_Sbus", &LSGrid::get_Sbus, DocLSGrid::get_Sbus.c_str())
        .def("get_dcSbus", &LSGrid::get_dcSbus, DocLSGrid::get_dcSbus.c_str())
        .def("get_Ybus_solver", &LSGrid::get_Ybus_solver, DocLSGrid::get_Ybus_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_dcYbus_solver", &LSGrid::get_dcYbus_solver, DocLSGrid::get_dcYbus_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_Sbus_solver", &LSGrid::get_Sbus_solver, DocLSGrid::get_Sbus_solver.c_str(), py::return_value_policy::reference_internal)
        .def("get_dcSbus_solver", &LSGrid::get_dcSbus_solver, DocLSGrid::get_dcSbus_solver.c_str(), py::return_value_policy::reference_internal)

        .def("check_solution", &LSGrid::check_solution, DocLSGrid::check_solution.c_str())

        .def("get_loads_res", &LSGrid::get_loads_res, DocLSGrid::get_loads_res.c_str(), py::return_value_policy::reference_internal)
        .def("get_loads_status", &LSGrid::get_loads_status, DocLSGrid::get_loads_status.c_str(), py::return_value_policy::reference_internal)
        .def("get_shunts_res", &LSGrid::get_shunts_res, DocLSGrid::get_shunts_res.c_str(), py::return_value_policy::reference_internal)
        .def("get_shunts_status", &LSGrid::get_shunts_status, DocLSGrid::get_shunts_status.c_str(), py::return_value_policy::reference_internal)
        .def("get_gen_res", &LSGrid::get_gen_res, DocLSGrid::get_gen_res.c_str(), py::return_value_policy::reference_internal)
        .def("get_gen_status", &LSGrid::get_gen_status, DocLSGrid::get_gen_status.c_str(), py::return_value_policy::reference_internal)
        .def("get_line_res1", &LSGrid::get_line_res1, DocLSGrid::get_line_res1.c_str(), py::return_value_policy::reference_internal)
        .def("get_line_res2", &LSGrid::get_line_res2, DocLSGrid::get_line_res2.c_str(), py::return_value_policy::reference_internal)
        .def("get_lines_status", &LSGrid::get_lines_status, DocLSGrid::get_lines_status.c_str(), py::return_value_policy::reference_internal)
        .def("get_lines_status_side1", &LSGrid::get_lines_status_side1, DocLSGrid::get_lines_status_side1.c_str(), py::return_value_policy::reference_internal)
        .def("get_lines_status_side2", &LSGrid::get_lines_status_side2, DocLSGrid::get_lines_status_side2.c_str(), py::return_value_policy::reference_internal)
        .def("get_trafo_res1", &LSGrid::get_trafo_res1, DocLSGrid::get_trafo_res1.c_str(), py::return_value_policy::reference_internal)
        .def("get_trafo_res2", &LSGrid::get_trafo_res2, DocLSGrid::get_trafo_res2.c_str(), py::return_value_policy::reference_internal)
        .def("get_trafo_status", &LSGrid::get_trafo_status, DocLSGrid::get_trafo_status.c_str(), py::return_value_policy::reference_internal)
        .def("get_trafo_status_side1", &LSGrid::get_trafo_status_side1, DocLSGrid::get_trafo_status_side1.c_str(), py::return_value_policy::reference_internal)
        .def("get_trafo_status_side2", &LSGrid::get_trafo_status_side2, DocLSGrid::get_trafo_status_side2.c_str(), py::return_value_policy::reference_internal)
        .def("get_storages_res", &LSGrid::get_storages_res, DocLSGrid::get_storages_res.c_str(), py::return_value_policy::reference_internal)
        .def("get_storages_status", &LSGrid::get_storages_status, DocLSGrid::get_storages_status.c_str(), py::return_value_policy::reference_internal)
        .def("get_sgens_res", &LSGrid::get_sgens_res, DocLSGrid::get_sgens_res.c_str(), py::return_value_policy::reference_internal)
        .def("get_sgens_status", &LSGrid::get_sgens_status, DocLSGrid::get_sgens_status.c_str(), py::return_value_policy::reference_internal)

        .def("get_gen_theta", &LSGrid::get_gen_theta, DocLSGrid::get_gen_theta.c_str(), py::return_value_policy::reference_internal)
        .def("get_load_theta", &LSGrid::get_load_theta, DocLSGrid::get_load_theta.c_str(), py::return_value_policy::reference_internal)
        .def("get_shunt_theta", &LSGrid::get_shunt_theta, DocLSGrid::get_shunt_theta.c_str(), py::return_value_policy::reference_internal)
        .def("get_storage_theta", &LSGrid::get_storage_theta, DocLSGrid::get_storage_theta.c_str(), py::return_value_policy::reference_internal)
        .def("get_line_theta1", &LSGrid::get_line_theta1, DocLSGrid::get_line_theta1.c_str(), py::return_value_policy::reference_internal)
        .def("get_line_theta2", &LSGrid::get_line_theta2, DocLSGrid::get_line_theta2.c_str(), py::return_value_policy::reference_internal)
        .def("get_trafo_theta1", &LSGrid::get_trafo_theta1, DocLSGrid::get_trafo_theta1.c_str(), py::return_value_policy::reference_internal)
        .def("get_trafo_theta2", &LSGrid::get_trafo_theta2, DocLSGrid::get_trafo_theta2.c_str(), py::return_value_policy::reference_internal)

        .def("get_all_shunt_buses", &LSGrid::get_all_shunt_buses_numpy, DocLSGrid::get_all_shunt_buses.c_str(), py::return_value_policy::reference_internal)
        .def("get_loads_res_full", &LSGrid::get_loads_res_full, DocLSGrid::get_loads_res_full.c_str(), py::return_value_policy::reference_internal)
        .def("get_shunts_res_full", &LSGrid::get_shunts_res_full, DocLSGrid::get_shunts_res_full.c_str(), py::return_value_policy::reference_internal)
        .def("get_gen_res_full", &LSGrid::get_gen_res_full, DocLSGrid::get_gen_res_full.c_str(), py::return_value_policy::reference_internal)
        .def("get_line_res1_full", &LSGrid::get_line_res1_full, DocLSGrid::get_line_res1_full.c_str(), py::return_value_policy::reference_internal)
        .def("get_line_res2_full", &LSGrid::get_line_res2_full, DocLSGrid::get_line_res2_full.c_str(), py::return_value_policy::reference_internal)
        .def("get_trafo_res1_full", &LSGrid::get_trafo_res1_full, DocLSGrid::get_trafo_res1_full.c_str(), py::return_value_policy::reference_internal)
        .def("get_trafo_res2_full", &LSGrid::get_trafo_res2_full, DocLSGrid::get_trafo_res2_full.c_str(), py::return_value_policy::reference_internal)
        .def("get_storages_res_full", &LSGrid::get_storages_res_full, DocLSGrid::get_storages_res_full.c_str(), py::return_value_policy::reference_internal)
        .def("get_sgens_res_full", &LSGrid::get_sgens_res_full, DocLSGrid::get_sgens_res_full.c_str(), py::return_value_policy::reference_internal)
        .def("get_dcline_res1_full", &LSGrid::get_dcline_res1_full, DocLSGrid::get_dcline_res1_full.c_str(), py::return_value_policy::reference_internal)
        .def("get_dcline_res2_full", &LSGrid::get_dcline_res2_full, DocLSGrid::get_dcline_res2_full.c_str(), py::return_value_policy::reference_internal)

        .def("get_shunt_target_p", &LSGrid::get_shunt_target_p, DocLSGrid::get_shunt_target_p.c_str(), py::return_value_policy::reference_internal)
        .def("get_load_target_p", &LSGrid::get_load_target_p, DocLSGrid::get_load_target_p.c_str(), py::return_value_policy::reference_internal)
        .def("get_gen_target_p", &LSGrid::get_gen_target_p, DocLSGrid::get_gen_target_p.c_str(), py::return_value_policy::reference_internal)
        .def("get_sgen_target_p", &LSGrid::get_sgen_target_p, DocLSGrid::get_sgen_target_p.c_str(), py::return_value_policy::reference_internal)
        .def("get_storage_target_p", &LSGrid::get_storage_target_p, DocLSGrid::get_storage_target_p.c_str(), py::return_value_policy::reference_internal);
}
