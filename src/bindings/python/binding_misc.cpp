// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "DataConverter.hpp"
#include "Utils.hpp"
#include "AlgoConfig.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

void bind_misc(py::module_& m) {
    py::class_<PandaPowerConverter>(m, "PandaPowerConverter", DocMisc::PandaPowerConverter.c_str())
        .def(py::init<>())
        .def("set_f_hz", &PandaPowerConverter::set_f_hz, DocMisc::set_f_hz.c_str())
        .def("set_sn_mva", &PandaPowerConverter::set_sn_mva, DocMisc::set_sn_mva.c_str())
        .def("get_line_param_legacy", &PandaPowerConverter::get_line_param_legacy, DocMisc::get_line_param_legacy.c_str())
        .def("get_line_param", &PandaPowerConverter::get_line_param, DocMisc::get_line_param.c_str())
        .def("get_trafo_param_pp3", &PandaPowerConverter::get_trafo_param_pp3, DocMisc::get_trafo_param_pp3.c_str())
        .def("get_trafo_param_pp2", &PandaPowerConverter::get_trafo_param_pp2, DocMisc::get_trafo_param_pp2.c_str());

    py::class_<AlgoConfig>(m, "AlgoConfig", DocMisc::AlgoConfig.c_str())
        .def(py::init<>())
        .def_readwrite("int_params",  &AlgoConfig::int_params, DocMisc::int_params.c_str())
        .def_readwrite("real_params", &AlgoConfig::real_params, DocMisc::real_params.c_str());

    py::class_<AlgoControl>(m, "AlgoControl", DocMisc::AlgoControl.c_str())
        .def(py::init<>())
        .def("has_dimension_changed", &AlgoControl::has_dimension_changed, DocMisc::has_dimension_changed.c_str())
        .def("has_pv_changed", &AlgoControl::has_pv_changed, DocMisc::has_pv_changed.c_str())
        .def("has_pq_changed", &AlgoControl::has_pq_changed, DocMisc::has_pq_changed.c_str())
        .def("has_slack_participate_changed", &AlgoControl::has_slack_participate_changed, DocMisc::has_slack_participate_changed.c_str())
        .def("need_reset_solver", &AlgoControl::need_reset_solver, DocMisc::need_reset_solver.c_str())
        .def("need_recompute_sbus", &AlgoControl::need_recompute_sbus, DocMisc::need_recompute_sbus.c_str())
        .def("need_recompute_ybus", &AlgoControl::need_recompute_ybus, DocMisc::need_recompute_ybus.c_str())
        .def("ybus_change_sparsity_pattern", &AlgoControl::ybus_change_sparsity_pattern, DocMisc::ybus_change_sparsity_pattern.c_str())
        .def("has_slack_weight_changed", &AlgoControl::has_slack_weight_changed, DocMisc::has_slack_weight_changed.c_str())
        .def("has_v_changed", &AlgoControl::has_v_changed, DocMisc::has_v_changed.c_str())
        .def("has_ybus_some_coeffs_zero", &AlgoControl::has_ybus_some_coeffs_zero, DocMisc::has_ybus_some_coeffs_zero.c_str())
        .def("has_one_el_changed_bus", &AlgoControl::has_one_el_changed_bus, DocMisc::has_one_el_changed_bus.c_str());
}
