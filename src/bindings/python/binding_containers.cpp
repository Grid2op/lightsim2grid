// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "pickle_helpers.hpp"
#include "element_container/GeneratorContainer.hpp"
#include "element_container/SGenContainer.hpp"
#include "element_container/LoadContainer.hpp"
#include "element_container/ShuntContainer.hpp"
#include "element_container/TrafoContainer.hpp"
#include "element_container/LineContainer.hpp"
#include "element_container/DCLineContainer.hpp"
#include "SubstationContainer.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

void bind_containers(py::module_& m) {
    auto gen_cls = py::class_<GeneratorContainer>(m, "GeneratorContainer", DocIterator::GeneratorContainer.c_str())
        .def("__len__", [](const GeneratorContainer & data) { return data.nb(); })
        .def("__getitem__", [](const GeneratorContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const GeneratorContainer & data)  {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def("get_bus_id", &GeneratorContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>());
    add_pickle(gen_cls, "GeneratorContainer");

    py::class_<GenInfo>(m, "GenInfo", DocIterator::GenInfo.c_str())
        .def_readonly("id", &GenInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &GenInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &GenInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &GenInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &GenInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &GenInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("is_slack", &GenInfo::is_slack, DocIterator::is_slack.c_str())
        .def_readonly("slack_weight", &GenInfo::slack_weight, DocIterator::slack_weight.c_str())
        .def_readonly("voltage_regulator_on", &GenInfo::voltage_regulator_on, "TODO")
        .def_readonly("target_p_mw", &GenInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_vm_pu", &GenInfo::target_vm_pu, DocIterator::target_vm_pu.c_str())
        .def_readonly("target_q_mvar", &GenInfo::target_q_mvar, "TODO")
        .def_readonly("min_q_mvar", &GenInfo::min_q_mvar, DocIterator::min_q_mvar.c_str())
        .def_readonly("max_q_mvar", &GenInfo::max_q_mvar, DocIterator::max_q_mvar.c_str())
        .def_readonly("has_res", &GenInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &GenInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &GenInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &GenInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &GenInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        .def_readonly("voltage_level_id", &GenInfo::sub_id, DocIterator::sub_id.c_str());

    auto sgen_cls = py::class_<SGenContainer>(m, "SGenContainer", DocIterator::SGenContainer.c_str())
        .def("__len__", [](const SGenContainer & data) { return data.nb(); })
        .def("__getitem__", [](const SGenContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const SGenContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def("get_bus_id", &SGenContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>());
    add_pickle(sgen_cls, "SGenContainer");

    py::class_<SGenInfo>(m, "SGenInfo", DocIterator::SGenInfo.c_str())
        .def_readonly("id", &SGenInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &SGenInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &SGenInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &SGenInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &SGenInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &SGenInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("min_q_mvar", &SGenInfo::min_q_mvar, DocIterator::min_q_mvar.c_str())
        .def_readonly("max_q_mvar", &SGenInfo::max_q_mvar, DocIterator::max_q_mvar.c_str())
        .def_readonly("min_p_mw", &SGenInfo::min_p_mw, DocIterator::min_p_mw.c_str())
        .def_readonly("max_p_mw", &SGenInfo::max_p_mw, DocIterator::max_p_mw.c_str())
        .def_readonly("target_p_mw", &SGenInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &SGenInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &SGenInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &SGenInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &SGenInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &SGenInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &SGenInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        .def_readonly("voltage_level_id", &SGenInfo::sub_id, DocIterator::sub_id.c_str());

    auto load_cls = py::class_<LoadContainer>(m, "LoadContainer", DocIterator::LoadContainer.c_str())
        .def("__len__", [](const LoadContainer & data) { return data.nb(); })
        .def("__getitem__", [](const LoadContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const LoadContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def("get_bus_id", &LoadContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>());
    add_pickle(load_cls, "LoadContainer");

    py::class_<LoadInfo>(m, "LoadInfo", DocIterator::LoadInfo.c_str())
        .def_readonly("id", &LoadInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &LoadInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &LoadInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &LoadInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &LoadInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &LoadInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("target_p_mw", &LoadInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &LoadInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &LoadInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &LoadInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &LoadInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &LoadInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &LoadInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        .def_readonly("voltage_level_id", &LoadInfo::sub_id, DocIterator::sub_id.c_str());

    auto shunt_cls = py::class_<ShuntContainer>(m, "ShuntContainer", DocIterator::ShuntContainer.c_str())
        .def("__len__", [](const ShuntContainer & data) { return data.nb(); })
        .def("__getitem__", [](const ShuntContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const ShuntContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def("get_bus_id", &ShuntContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>());
    add_pickle(shunt_cls, "ShuntContainer");

    py::class_<ShuntInfo>(m, "ShuntInfo", DocIterator::ShuntInfo.c_str())
        .def_readonly("id", &ShuntInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &ShuntInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &ShuntInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &ShuntInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &ShuntInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &ShuntInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("target_p_mw", &ShuntInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &ShuntInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &ShuntInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &ShuntInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &ShuntInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &ShuntInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &ShuntInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        .def_readonly("voltage_level_id", &ShuntInfo::sub_id, DocIterator::sub_id.c_str());

    auto trafo_cls = py::class_<TrafoContainer>(m, "TrafoContainer", DocIterator::TrafoContainer.c_str())
        .def("__len__", [](const TrafoContainer & data) { return data.nb(); })
        .def("__getitem__", [](const TrafoContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const TrafoContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def_property_readonly("ignore_tap_side_for_shift", &TrafoContainer::ignore_tap_side_for_shift,
            R"mydelimiter(
            Whether ignore the tap side is ignored when using the
            'shift' attribute (should be True for pandapower,
            where it is ignored and False otherwise).)mydelimiter")
        .def("get_bus_id_side_1", &TrafoContainer::get_bus_id_side_1_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def("get_bus_id_side_2", &TrafoContainer::get_bus_id_side_2_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def("get_yac_eff_11", [](const TrafoContainer & t) -> Eigen::Ref<const CplxVect> { return t.yac_eff_11(); }, "TODO doc", py::keep_alive<0, 1>())
        .def("get_yac_eff_12", [](const TrafoContainer & t) -> Eigen::Ref<const CplxVect> { return t.yac_eff_12(); }, "TODO doc", py::keep_alive<0, 1>())
        .def("get_yac_eff_21", [](const TrafoContainer & t) -> Eigen::Ref<const CplxVect> { return t.yac_eff_21(); }, "TODO doc", py::keep_alive<0, 1>())
        .def("get_yac_eff_22", [](const TrafoContainer & t) -> Eigen::Ref<const CplxVect> { return t.yac_eff_22(); }, "TODO doc", py::keep_alive<0, 1>());
    add_pickle(trafo_cls, "TrafoContainer");

    py::class_<TrafoInfo>(m, "TrafoInfo", DocIterator::TrafoInfo.c_str())
        .def_readonly("id", &TrafoInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &TrafoInfo::name, DocIterator::name.c_str())
        .def_readonly("sub1_id", &TrafoInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("sub2_id", &TrafoInfo::sub_2_id, DocIterator::sub_id.c_str())
        .def_readonly("pos1_topo_vect", &TrafoInfo::pos_1_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("pos2_topo_vect", &TrafoInfo::pos_2_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected_global", &TrafoInfo::connected_global, DocIterator::connected.c_str())
        .def_readonly("connected1", &TrafoInfo::connected_1, DocIterator::connected.c_str())
        .def_readonly("connected2", &TrafoInfo::connected_2, DocIterator::connected.c_str())
        .def_readonly("bus1_id", &TrafoInfo::bus_1_id, DocIterator::bus_hv_id.c_str())
        .def_readonly("bus2_id", &TrafoInfo::bus_2_id, DocIterator::bus_lv_id.c_str())
        .def_readonly("r_pu", &TrafoInfo::r_pu, DocIterator::r_pu.c_str())
        .def_readonly("x_pu", &TrafoInfo::x_pu, DocIterator::x_pu.c_str())
        .def_readonly("h1_pu", &TrafoInfo::h1_pu, DocIterator::h_pu.c_str())
        .def_readonly("h2_pu", &TrafoInfo::h2_pu, DocIterator::h_pu.c_str())
        .def_readonly("is_tap_side_1", &TrafoInfo::is_tap_side1, DocIterator::is_tap_hv_side.c_str())
        .def_readonly("ratio", &TrafoInfo::ratio, DocIterator::ratio.c_str())
        .def_readonly("shift_rad", &TrafoInfo::shift_rad, DocIterator::shift_rad.c_str())
        .def_readonly("has_res", &TrafoInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p1_mw", &TrafoInfo::res_p1_mw, DocIterator::res_p_hv_mw.c_str())
        .def_readonly("res_q1_mvar", &TrafoInfo::res_q1_mvar, DocIterator::res_q_hv_mvar.c_str())
        .def_readonly("res_v1_kv", &TrafoInfo::res_v1_kv, DocIterator::res_v_hv_kv.c_str())
        .def_readonly("res_a1_ka", &TrafoInfo::res_a1_ka, DocIterator::res_a_hv_ka.c_str())
        .def_readonly("res_p2_mw", &TrafoInfo::res_p2_mw, DocIterator::res_p_lv_mw.c_str())
        .def_readonly("res_q2_mvar", &TrafoInfo::res_q2_mvar, DocIterator::res_q_lv_mvar.c_str())
        .def_readonly("res_v2_kv", &TrafoInfo::res_v2_kv, DocIterator::res_v_lv_kv.c_str())
        .def_readonly("res_a2_ka", &TrafoInfo::res_a2_ka, DocIterator::res_a_lv_ka.c_str())
        .def_readonly("res_theta1_deg", &TrafoInfo::res_theta1_deg, DocIterator::res_theta_hv_deg.c_str())
        .def_readonly("res_theta2_deg", &TrafoInfo::res_theta2_deg, DocIterator::res_theta_lv_deg.c_str())
        .def_readonly("yac_11", &TrafoInfo::yac_11, "TODO doc")
        .def_readonly("yac_12", &TrafoInfo::yac_12, "TODO doc")
        .def_readonly("yac_21", &TrafoInfo::yac_21, "TODO doc")
        .def_readonly("yac_22", &TrafoInfo::yac_22, "TODO doc")
        .def_readonly("yac_eff_11", &TrafoInfo::yac_eff_11, "TODO doc")
        .def_readonly("yac_eff_12", &TrafoInfo::yac_eff_12, "TODO doc")
        .def_readonly("yac_eff_21", &TrafoInfo::yac_eff_21, "TODO doc")
        .def_readonly("yac_eff_22", &TrafoInfo::yac_eff_22, "TODO doc")
        .def_readonly("ydc_11", &TrafoInfo::ydc_11, "TODO doc")
        .def_readonly("ydc_12", &TrafoInfo::ydc_12, "TODO doc")
        .def_readonly("ydc_21", &TrafoInfo::ydc_21, "TODO doc")
        .def_readonly("ydc_22", &TrafoInfo::ydc_22, "TODO doc")
        .def_readonly("voltage_level1_id", &TrafoInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("voltage_level2_id", &TrafoInfo::sub_2_id, DocIterator::sub_id.c_str());

    auto line_cls = py::class_<LineContainer>(m, "LineContainer", DocIterator::LineContainer.c_str())
        .def("__len__", [](const LineContainer & data) { return data.nb(); })
        .def("__getitem__", [](const LineContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const LineContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def("get_bus_id_side_1", &LineContainer::get_bus_id_side_1_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def("get_bus_id_side_2", &LineContainer::get_bus_id_side_2_numpy, "TODO doc", py::keep_alive<0, 1>())
        .def("get_yac_eff_11", [](const LineContainer & l) -> Eigen::Ref<const CplxVect> { return l.yac_eff_11(); }, "TODO doc", py::keep_alive<0, 1>())
        .def("get_yac_eff_12", [](const LineContainer & l) -> Eigen::Ref<const CplxVect> { return l.yac_eff_12(); }, "TODO doc", py::keep_alive<0, 1>())
        .def("get_yac_eff_21", [](const LineContainer & l) -> Eigen::Ref<const CplxVect> { return l.yac_eff_21(); }, "TODO doc", py::keep_alive<0, 1>())
        .def("get_yac_eff_22", [](const LineContainer & l) -> Eigen::Ref<const CplxVect> { return l.yac_eff_22(); }, "TODO doc", py::keep_alive<0, 1>());
    add_pickle(line_cls, "LineContainer");

    py::class_<LineInfo>(m, "LineInfo", DocIterator::LineInfo.c_str())
        .def_readonly("id", &LineInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &LineInfo::name, DocIterator::name.c_str())
        .def_readonly("sub1_id", &LineInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("sub2_id", &LineInfo::sub_2_id, DocIterator::sub_id.c_str())
        .def_readonly("pos1_topo_vect", &LineInfo::pos_1_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("pos2_topo_vect", &LineInfo::pos_2_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected_global", &LineInfo::connected_global, DocIterator::connected.c_str())
        .def_readonly("connected1", &LineInfo::connected_1, DocIterator::connected.c_str())
        .def_readonly("connected2", &LineInfo::connected_2, DocIterator::connected.c_str())
        .def_readonly("bus1_id", &LineInfo::bus_1_id, DocIterator::bus_or_id.c_str())
        .def_readonly("bus2_id", &LineInfo::bus_2_id, DocIterator::bus_ex_id.c_str())
        .def_readonly("r_pu", &LineInfo::r_pu, DocIterator::r_pu.c_str())
        .def_readonly("x_pu", &LineInfo::x_pu, DocIterator::x_pu.c_str())
        .def_readonly("h1_pu", &LineInfo::h1_pu, DocIterator::h_pu.c_str())
        .def_readonly("h2_pu", &LineInfo::h2_pu, DocIterator::h_pu.c_str())
        .def_readonly("has_res", &LineInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p1_mw", &LineInfo::res_p1_mw, DocIterator::res_p_or_mw.c_str())
        .def_readonly("res_q1_mvar", &LineInfo::res_q1_mvar, DocIterator::res_q_or_mvar.c_str())
        .def_readonly("res_v1_kv", &LineInfo::res_v1_kv, DocIterator::res_v_or_kv.c_str())
        .def_readonly("res_a1_ka", &LineInfo::res_a1_ka, DocIterator::res_a_or_ka.c_str())
        .def_readonly("res_p2_mw", &LineInfo::res_p2_mw, DocIterator::res_p_ex_mw.c_str())
        .def_readonly("res_q2_mvar", &LineInfo::res_q2_mvar, DocIterator::res_q_ex_mvar.c_str())
        .def_readonly("res_v2_kv", &LineInfo::res_v2_kv, DocIterator::res_v_ex_kv.c_str())
        .def_readonly("res_a2_ka", &LineInfo::res_a2_ka, DocIterator::res_a_ex_ka.c_str())
        .def_readonly("res_theta1_deg", &LineInfo::res_theta1_deg, DocIterator::res_theta_or_deg.c_str())
        .def_readonly("res_theta2_deg", &LineInfo::res_theta2_deg, DocIterator::res_theta_ex_deg.c_str())
        .def_readonly("yac_11", &LineInfo::yac_11, "TODO doc")
        .def_readonly("yac_12", &LineInfo::yac_12, "TODO doc")
        .def_readonly("yac_21", &LineInfo::yac_21, "TODO doc")
        .def_readonly("yac_22", &LineInfo::yac_22, "TODO doc")
        .def_readonly("yac_eff_11", &LineInfo::yac_eff_11, "TODO doc")
        .def_readonly("yac_eff_12", &LineInfo::yac_eff_12, "TODO doc")
        .def_readonly("yac_eff_21", &LineInfo::yac_eff_21, "TODO doc")
        .def_readonly("yac_eff_22", &LineInfo::yac_eff_22, "TODO doc")
        .def_readonly("ydc_11", &LineInfo::ydc_11, "TODO doc")
        .def_readonly("ydc_12", &LineInfo::ydc_12, "TODO doc")
        .def_readonly("ydc_21", &LineInfo::ydc_21, "TODO doc")
        .def_readonly("ydc_22", &LineInfo::ydc_22, "TODO doc")
        .def_readonly("voltage_level1_id", &LineInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("voltage_level2_id", &LineInfo::sub_2_id, DocIterator::sub_id.c_str());

    auto dcline_cls = py::class_<DCLineContainer>(m, "DCLineContainer", DocIterator::DCLineContainer.c_str())
        .def("__len__", [](const DCLineContainer & data) { return data.nb(); })
        .def("__getitem__", [](const DCLineContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const DCLineContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def("get_bus_id_side_1", &DCLineContainer::get_bus_id_side_1_numpy)
        .def("get_bus_id_side_2", &DCLineContainer::get_bus_id_side_2_numpy);
    add_pickle(dcline_cls, "DCLineContainer");

    py::class_<DCLineInfo>(m, "DCLineInfo", DocIterator::DCLineInfo.c_str())
        .def_readonly("id", &DCLineInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &DCLineInfo::name, DocIterator::name.c_str())
        .def_readonly("sub1_id", &DCLineInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("sub2_id", &DCLineInfo::sub_2_id, DocIterator::sub_id.c_str())
        .def_readonly("pos1_topo_vect", &DCLineInfo::pos_1_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("pos2_topo_vect", &DCLineInfo::pos_2_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected_global", &DCLineInfo::connected_global, DocIterator::connected.c_str())
        .def_readonly("connected1", &DCLineInfo::connected_1, DocIterator::connected.c_str())
        .def_readonly("connected2", &DCLineInfo::connected_2, DocIterator::connected.c_str())
        .def_readonly("bus1_id", &DCLineInfo::bus_1_id, DocIterator::bus_or_id.c_str())
        .def_readonly("bus2_id", &DCLineInfo::bus_2_id, DocIterator::bus_ex_id.c_str())
        .def_readonly("target_p1_mw", &DCLineInfo::target_p_1_mw, DocIterator::target_p_or_mw.c_str())
        .def_readonly("p2_mw", &DCLineInfo::p_2_mw, DocIterator::target_p_or_mw.c_str())
        .def_readonly("target_vm1_pu", &DCLineInfo::target_vm_1_pu, DocIterator::target_vm_or_pu.c_str())
        .def_readonly("target_vm2_pu", &DCLineInfo::target_vm_2_pu, DocIterator::target_vm_ex_pu.c_str())
        .def_readonly("loss_pct", &DCLineInfo::loss_pct, DocIterator::loss_pct.c_str())
        .def_readonly("loss_mw", &DCLineInfo::loss_mw, DocIterator::loss_mw.c_str())
        .def_readonly("gen1", &DCLineInfo::gen_side_1, DocIterator::gen_or.c_str())
        .def_readonly("gen2", &DCLineInfo::gen_side_2, DocIterator::gen_ex.c_str())
        .def_readonly("has_res", &DCLineInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p1_mw", &DCLineInfo::res_p1_mw, DocIterator::res_p_or_mw_dcline.c_str())
        .def_readonly("res_p2_mw", &DCLineInfo::res_p2_mw, DocIterator::res_p_ex_mw_dcline.c_str())
        .def_readonly("res_q1_mvar", &DCLineInfo::res_q1_mvar, DocIterator::res_q_or_mvar_dcline.c_str())
        .def_readonly("res_q2_mvar", &DCLineInfo::res_q2_mvar, DocIterator::res_q_ex_mvar_dcline.c_str())
        .def_readonly("res_v1_kv", &DCLineInfo::res_v1_kv, DocIterator::res_v_or_kv_dcline.c_str())
        .def_readonly("res_v2_kv", &DCLineInfo::res_v2_kv, DocIterator::res_v_ex_kv_dcline.c_str())
        .def_readonly("res_theta1_deg", &DCLineInfo::res_theta1_deg, DocIterator::res_theta_or_deg_dcline.c_str())
        .def_readonly("res_theta2_deg", &DCLineInfo::res_theta2_deg, DocIterator::res_theta_ex_deg_dcline.c_str())
        .def_readonly("voltage_level1_id", &DCLineInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("voltage_level2_id", &DCLineInfo::sub_2_id, DocIterator::sub_id.c_str());

    auto sub_cls = py::class_<SubstationContainer>(m, "SubstationContainer", "TODO")
        .def("__len__", [](const SubstationContainer & data) { return data.nb(); })
        .def("__getitem__", [](const SubstationContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const SubstationContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>());
    add_pickle(sub_cls, "SubstationContainer");

    py::class_<SubstationInfo>(m, "SubstationInfo", "TODO")
        .def_readonly("id", &SubstationInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &SubstationInfo::name, DocIterator::name.c_str())
        .def_readonly("nb_max_busbars", &SubstationInfo::nb_max_busbars, DocIterator::name.c_str())
        .def_readonly("vn_kv", &SubstationInfo::vn_kv, DocIterator::name.c_str());
}
