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

using namespace ls2g;

void bind_misc(py::module_& m) {
    py::class_<PandaPowerConverter>(m, "PandaPowerConverter")
        .def(py::init<>())
        .def("set_f_hz", &PandaPowerConverter::set_f_hz)
        .def("set_sn_mva", &PandaPowerConverter::set_sn_mva)
        .def("get_line_param_legacy", &PandaPowerConverter::get_line_param_legacy)
        .def("get_line_param", &PandaPowerConverter::get_line_param)
        .def("get_trafo_param_pp3", &PandaPowerConverter::get_trafo_param_pp3)
        .def("get_trafo_param_pp2", &PandaPowerConverter::get_trafo_param_pp2);

    py::class_<AlgoConfig>(m, "AlgoConfig",
        "Serialisable configuration blob for Newton-Raphson algorithm parameters.\n"
        "Stores scaling/refactor policy type and all associated parameters.\n"
        "Obtain via solver.get_config() or gridmodel.get_ac_algo_config().")
        .def(py::init<>())
        .def_readwrite("int_params",  &AlgoConfig::int_params,
            "Integer parameters: [ScalingPolicyType, RefactorPolicyType, ls_max_iter, refactor_every_n]")
        .def_readwrite("real_params", &AlgoConfig::real_params,
            "Real parameters: [max_dVa, max_dVm, ls_c, ls_rho, iw_mu_min, iw_mu_max]");

    py::class_<AlgoControl>(m, "AlgoControl", "TODO")
        .def(py::init<>())
        .def("has_dimension_changed", &AlgoControl::has_dimension_changed, "TODO")
        .def("has_pv_changed", &AlgoControl::has_pv_changed, "TODO")
        .def("has_pq_changed", &AlgoControl::has_pq_changed, "TODO")
        .def("has_slack_participate_changed", &AlgoControl::has_slack_participate_changed, "TODO")
        .def("need_reset_solver", &AlgoControl::need_reset_solver, "TODO")
        .def("need_recompute_sbus", &AlgoControl::need_recompute_sbus, "TODO")
        .def("need_recompute_ybus", &AlgoControl::need_recompute_ybus, "TODO")
        .def("ybus_change_sparsity_pattern", &AlgoControl::ybus_change_sparsity_pattern, "TODO")
        .def("has_slack_weight_changed", &AlgoControl::has_slack_weight_changed, "TODO")
        .def("has_v_changed", &AlgoControl::has_v_changed, "TODO")
        .def("has_ybus_some_coeffs_zero", &AlgoControl::has_ybus_some_coeffs_zero, "TODO")
        .def("has_one_el_changed_bus", &AlgoControl::has_one_el_changed_bus, "TODO");
}
