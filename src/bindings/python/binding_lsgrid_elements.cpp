// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


// The element accessors, mutators and names of LSGrid -- one of the four translation units the LSGrid bindings are spread
// over, see the note at the top of binding_lsgrid.cpp.

#include "binding_declarations.hpp"
#include "LSGrid.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

void bind_lsgrid_elements(py::class_<LSGrid> & cls) {
    cls
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
        .def("set_gen_p_limits", &LSGrid::set_gen_p_limits,
             py::arg("p_min_mw"), py::arg("p_max_mw"),
             "Active power limits (MW) of the generators, OPTIONAL and never enforced -- the "
             "same shape as the current limits above: one entry per generator, NaN where a "
             "machine has none, and two empty vectors to drop them again (`GenInfo.min_p_mw` "
             "/ `max_p_mw` then read NaN).\n\n"
             "What they are for: the distributed slack is solved INSIDE the Newton system, by "
             "fixed participation factors that know nothing about limits, so a participating "
             "machine's converged active power -- its target plus its share of the imbalance "
             "-- can land beyond what it can deliver. That is a physical violation, and these "
             "are what the batch algorithms' `compute_physical_violations` compares against "
             "(LOW_P / HIGH_P on the GENERATOR).")
        .def("set_gen_names", &LSGrid::set_gen_names, DocLSGrid::set_gen_names.c_str())
        .def("set_load_names", &LSGrid::set_load_names, DocLSGrid::set_load_names.c_str())
        .def("set_storage_names", &LSGrid::set_storage_names, DocLSGrid::set_storage_names.c_str())
        .def("set_sgen_names", &LSGrid::set_sgen_names, DocLSGrid::set_sgen_names.c_str())
        .def("set_shunt_names", &LSGrid::set_shunt_names, DocLSGrid::set_shunt_names.c_str())
        .def("set_substation_names", &LSGrid::set_substation_names, DocLSGrid::set_substation_names.c_str())
        .def("get_substation_names", &LSGrid::get_substation_names, DocLSGrid::get_substation_names.c_str())

        // deprecated no-ops since 1.0.0: a bus is in the solved system iff an active element
        // sits on it, so there is no separate switch left for these to flip. Kept so that
        // existing loaders (pandapower / powermodels) and backends keep importing.
        .def("deactivate_bus", &LSGrid::deactivate_bus_python, DocLSGrid::deactivate_bus.c_str())
        .def("reactivate_bus", &LSGrid::reactivate_bus_python, DocLSGrid::reactivate_bus.c_str())

        .def("deactivate_powerline", &LSGrid::deactivate_powerline, DocLSGrid::deactivate_powerline.c_str())
        .def("reactivate_powerline", &LSGrid::reactivate_powerline, DocLSGrid::reactivate_powerline.c_str())
        .def("deactivate_powerline_side1", &LSGrid::deactivate_powerline_side1, DocLSGrid::deactivate_powerline_side1.c_str())
        .def("deactivate_powerline_side2", &LSGrid::deactivate_powerline_side2, DocLSGrid::deactivate_powerline_side2.c_str())
        .def("reactivate_powerline_side1", &LSGrid::reactivate_powerline_side1, DocLSGrid::reactivate_powerline_side1.c_str())
        .def("reactivate_powerline_side2", &LSGrid::reactivate_powerline_side2, DocLSGrid::reactivate_powerline_side2.c_str())
        .def("change_bus1_powerline", &LSGrid::change_bus1_powerline_python, DocLSGrid::change_bus1_powerline.c_str())
        .def("change_bus2_powerline", &LSGrid::change_bus2_powerline_python, DocLSGrid::change_bus2_powerline.c_str())
        .def("get_bus1_powerline", &LSGrid::get_bus1_powerline, DocLSGrid::get_bus1_powerline.c_str(), py::return_value_policy::reference_internal)
        .def("get_bus2_powerline", &LSGrid::get_bus2_powerline, DocLSGrid::get_bus2_powerline.c_str(), py::return_value_policy::reference_internal)

        .def("deactivate_trafo", &LSGrid::deactivate_trafo, DocLSGrid::deactivate_trafo.c_str())
        .def("reactivate_trafo", &LSGrid::reactivate_trafo, DocLSGrid::reactivate_trafo.c_str())
        .def("deactivate_trafo_side1", &LSGrid::deactivate_trafo_side1, DocLSGrid::deactivate_trafo_side1.c_str())
        .def("deactivate_trafo_side2", &LSGrid::deactivate_trafo_side2, DocLSGrid::deactivate_trafo_side2.c_str())
        .def("reactivate_trafo_side1", &LSGrid::reactivate_trafo_side1, DocLSGrid::reactivate_trafo_side1.c_str())
        .def("reactivate_trafo_side2", &LSGrid::reactivate_trafo_side2, DocLSGrid::reactivate_trafo_side2.c_str())
        .def("change_bus1_trafo", &LSGrid::change_bus1_trafo_python, DocLSGrid::change_bus1_trafo.c_str())
        .def("change_bus2_trafo", &LSGrid::change_bus2_trafo_python, DocLSGrid::change_bus2_trafo.c_str())
        .def("get_bus1_trafo", &LSGrid::get_bus1_trafo, DocLSGrid::get_bus1_trafo.c_str(), py::return_value_policy::reference_internal)
        .def("get_bus2_trafo", &LSGrid::get_bus2_trafo, DocLSGrid::get_bus2_trafo.c_str(), py::return_value_policy::reference_internal)
        .def("change_ratio_trafo", &LSGrid::change_ratio_trafo, DocLSGrid::change_ratio_trafo.c_str())
        .def("change_shift_trafo", &LSGrid::change_shift_trafo, DocLSGrid::change_shift_trafo.c_str())
        .def("change_shift_trafo_deg", &LSGrid::change_shift_trafo_deg, DocLSGrid::change_shift_trafo_deg.c_str())
        .def("set_trafo_shift_dependent_rx", &LSGrid::set_trafo_shift_dependent_rx,
            py::arg("enable"), py::arg("alpha_rad"), py::arg("rx_corr_pct"),
            DocLSGrid::set_trafo_shift_dependent_rx.c_str())
        .def("deactivate_load", &LSGrid::deactivate_load, DocLSGrid::deactivate_load.c_str())
        .def("reactivate_load", &LSGrid::reactivate_load, DocLSGrid::reactivate_load.c_str())
        .def("change_bus_load", &LSGrid::change_bus_load_python, DocLSGrid::change_bus_load.c_str())
        .def("get_bus_load", &LSGrid::get_bus_load, DocLSGrid::get_bus_load.c_str(), py::return_value_policy::reference_internal)
        .def("change_p_load", &LSGrid::change_p_load, DocLSGrid::change_p_load.c_str())
        .def("change_q_load", &LSGrid::change_q_load, DocLSGrid::change_q_load.c_str())

        .def("deactivate_gen", &LSGrid::deactivate_gen, DocLSGrid::deactivate_gen.c_str())
        .def("reactivate_gen", &LSGrid::reactivate_gen, DocLSGrid::reactivate_gen.c_str())
        .def("change_bus_gen", &LSGrid::change_bus_gen_python, DocLSGrid::change_bus_gen.c_str())
        .def("get_bus_gen", &LSGrid::get_bus_gen, DocLSGrid::get_bus_gen.c_str(), py::return_value_policy::reference_internal)
        .def("change_p_gen", &LSGrid::change_p_gen, DocLSGrid::change_p_gen.c_str())
        .def("change_v_gen", &LSGrid::change_v_gen, DocLSGrid::change_v_gen.c_str())
        .def("set_gen_regulated_bus", &LSGrid::set_gen_regulated_bus, DocLSGrid::set_gen_regulated_bus.c_str())
        .def("set_gen_reactive_key", &LSGrid::set_gen_reactive_key, DocLSGrid::set_gen_reactive_key.c_str())
        .def("deactivate_svc", &LSGrid::deactivate_svc, DocLSGrid::deactivate_svc.c_str())
        .def("reactivate_svc", &LSGrid::reactivate_svc, DocLSGrid::reactivate_svc.c_str())
        .def("change_bus_svc", &LSGrid::change_bus_svc_python, DocLSGrid::change_bus_svc.c_str())
        .def("get_bus_svc", &LSGrid::get_bus_svc, DocLSGrid::get_bus_svc.c_str())
        .def("set_svc_names", &LSGrid::set_svc_names, DocLSGrid::set_svc_names.c_str())

        .def("deactivate_shunt", &LSGrid::deactivate_shunt, DocLSGrid::deactivate_shunt.c_str())
        .def("reactivate_shunt", &LSGrid::reactivate_shunt, DocLSGrid::reactivate_shunt.c_str())
        .def("change_bus_shunt", &LSGrid::change_bus_shunt_python, DocLSGrid::change_bus_shunt.c_str())
        .def("get_bus_shunt", &LSGrid::get_bus_shunt, DocLSGrid::get_bus_shunt.c_str(), py::return_value_policy::reference_internal)
        .def("change_p_shunt", &LSGrid::change_p_shunt, DocLSGrid::change_p_shunt.c_str())
        .def("change_q_shunt", &LSGrid::change_q_shunt, DocLSGrid::change_q_shunt.c_str())

        .def("deactivate_sgen", &LSGrid::deactivate_sgen, DocLSGrid::deactivate_sgen.c_str())
        .def("reactivate_sgen", &LSGrid::reactivate_sgen, DocLSGrid::reactivate_sgen.c_str())
        .def("change_bus_sgen", &LSGrid::change_bus_sgen_python, DocLSGrid::change_bus_sgen.c_str())
        .def("get_bus_sgen", &LSGrid::get_bus_sgen, DocLSGrid::get_bus_sgen.c_str(), py::return_value_policy::reference_internal)
        .def("change_p_sgen", &LSGrid::change_p_sgen, DocLSGrid::change_p_sgen.c_str())
        .def("change_q_sgen", &LSGrid::change_q_sgen, DocLSGrid::change_q_sgen.c_str())

        .def("deactivate_storage", &LSGrid::deactivate_storage, DocLSGrid::deactivate_storage.c_str())
        .def("reactivate_storage", &LSGrid::reactivate_storage, DocLSGrid::reactivate_storage.c_str())
        .def("change_bus_storage", &LSGrid::change_bus_storage_python, DocLSGrid::change_bus_storage.c_str())
        .def("get_bus_storage", &LSGrid::get_bus_storage, DocLSGrid::get_bus_storage.c_str(), py::return_value_policy::reference_internal)
        .def("change_p_storage", &LSGrid::change_p_storage, DocLSGrid::change_p_storage.c_str())
        .def("change_q_storage", &LSGrid::change_q_storage, DocLSGrid::change_q_storage.c_str())
        .def("change_v_storage", &LSGrid::change_v_storage, DocLSGrid::change_v_storage.c_str())

        .def("deactivate_dcline", &LSGrid::deactivate_dcline, DocLSGrid::deactivate_dcline.c_str())
        .def("deactivate_dcline_side1", &LSGrid::deactivate_dcline_side1, DocLSGrid::deactivate_dcline_side1.c_str())
        .def("deactivate_dcline_side2", &LSGrid::deactivate_dcline_side2, DocLSGrid::deactivate_dcline_side2.c_str())
        .def("reactivate_dcline", &LSGrid::reactivate_dcline, DocLSGrid::reactivate_dcline.c_str())
        .def("change_p_dcline", &LSGrid::change_p_dcline, DocLSGrid::change_p_dcline.c_str())
        .def("change_v1_dcline", &LSGrid::change_v1_dcline, DocLSGrid::change_v1_dcline.c_str())
        .def("change_v2_dcline", &LSGrid::change_v2_dcline, DocLSGrid::change_v2_dcline.c_str())
        .def("change_bus1_dcline", &LSGrid::change_bus1_dcline, DocLSGrid::change_bus1_dcline.c_str())
        .def("change_bus2_dcline", &LSGrid::change_bus2_dcline, DocLSGrid::change_bus2_dcline.c_str())
        .def("get_bus1_dcline", &LSGrid::get_bus1_dcline, DocLSGrid::get_bus1_dcline.c_str(), py::return_value_policy::reference_internal)
        .def("get_bus2_dcline", &LSGrid::get_bus2_dcline, DocLSGrid::get_bus2_dcline.c_str(), py::return_value_policy::reference_internal)
        .def("set_status_droop_hvdc", &LSGrid::set_status_droop_hvdc, DocLSGrid::set_status_droop_hvdc.c_str())
        .def("get_status_droop_hvdc", &LSGrid::get_status_droop_hvdc, DocLSGrid::get_status_droop_hvdc.c_str())
        .def("get_status_droop_hvdc_vect", &LSGrid::get_status_droop_hvdc_vect, DocLSGrid::get_status_droop_hvdc_vect.c_str());
}
