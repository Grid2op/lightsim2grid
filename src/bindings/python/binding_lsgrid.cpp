// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


// The LSGrid bindings are spread over five translation units -- this one (the
// class itself, its construction, the solve and the solver-side accessors),
// binding_lsgrid_state.cpp (pickle and the binary format),
// binding_lsgrid_elements.cpp (the element accessors, mutators and names),
// binding_lsgrid_results.cpp (the results) and binding_lsgrid_topology.cpp (the
// grid2op vectors, the substation ids and the switches). As one file, GCC needed
// ~2.7 GB to compile it at -O3, which the 4 GB CI containers could not always
// spare; 2.2 GB of that was the pickle alone (see binding_lsgrid_state.cpp for
// why, and how it is now ~0.9 GB). Each unit takes a py::class_<LSGrid> & and
// appends its group; the split is by topic only.
#include "binding_declarations.hpp"
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
    // pickle and the binary format: binding_lsgrid_state.cpp
    bind_lsgrid_state(lsgrid_cls);
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
        .def("init_storages_full", &LSGrid::init_storages_full, DocLSGrid::init_storages_full.c_str())
        .def("init_sgens", &LSGrid::init_sgens, DocLSGrid::init_sgens.c_str())
        .def("init_dclines", &LSGrid::init_dclines, DocLSGrid::init_dclines.c_str())
        .def("init_hvdc_lines", &LSGrid::init_hvdc_lines, DocLSGrid::init_hvdc_lines.c_str())
        .def("init_svcs", &LSGrid::init_svcs, DocLSGrid::init_svcs.c_str())
        .def("add_gen_slackbus", &LSGrid::add_gen_slackbus, DocLSGrid::add_gen_slackbus.c_str())
        .def("remove_gen_slackbus", &LSGrid::remove_gen_slackbus, DocLSGrid::remove_gen_slackbus.c_str())
        .def("add_storage_slackbus", &LSGrid::add_storage_slackbus, DocLSGrid::add_storage_slackbus.c_str())
        .def("remove_storage_slackbus", &LSGrid::remove_storage_slackbus, DocLSGrid::remove_storage_slackbus.c_str())
        .def("get_bus_vn_kv", &LSGrid::get_bus_vn_kv, DocLSGrid::get_bus_vn_kv.c_str(), py::return_value_policy::reference_internal)
        // NB no return_value_policy::reference: get_bus_status() now BUILDS the vector from the
        // per-bus element counts and returns it by value, so pybind must copy (the default).
        .def("get_bus_status", &LSGrid::get_bus_status, DocLSGrid::get_bus_status.c_str())
        .def("set_bus_voltage_limits", &LSGrid::set_bus_voltage_limits, DocLSGrid::set_bus_voltage_limits.c_str())
        .def("get_bus_vmin_kv", &LSGrid::get_bus_vmin_kv, DocLSGrid::get_bus_vmin_kv.c_str(), py::return_value_policy::reference_internal)
        .def("get_bus_vmax_kv", &LSGrid::get_bus_vmax_kv, DocLSGrid::get_bus_vmax_kv.c_str(), py::return_value_policy::reference_internal)

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
        .def("get_ac_algo_controler", &LSGrid::get_ac_algo_controler, DocLSGrid::get_ac_algo_controler.c_str(), py::return_value_policy::reference_internal)
        .def("get_dc_algo_controler", &LSGrid::get_dc_algo_controler, DocLSGrid::get_dc_algo_controler.c_str(), py::return_value_policy::reference_internal)
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

        // debug functions (might disappear without further notice)
        .def("debug_get_Bp_python", &LSGrid::debug_get_Bp_python, DocLSGrid::_internal_do_not_use.c_str())
        .def("debug_get_Bpp_python", &LSGrid::debug_get_Bpp_python, DocLSGrid::_internal_do_not_use.c_str());

    // the other groups, each in its own translation unit (see the note above)
    bind_lsgrid_elements(lsgrid_cls);
    bind_lsgrid_results(lsgrid_cls);
    bind_lsgrid_topology(lsgrid_cls);
}
