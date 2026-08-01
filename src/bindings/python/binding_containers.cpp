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
#include "element_container/GeneratorContainer.hpp"
#include "element_container/SGenContainer.hpp"
#include "element_container/SvcContainer.hpp"
#include "element_container/LoadContainer.hpp"
#include "element_container/StorageContainer.hpp"
#include "element_container/ShuntContainer.hpp"
#include "element_container/TrafoContainer.hpp"
#include "element_container/LineContainer.hpp"
#include "element_container/HvdcLineContainer.hpp"
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
    add_binary_serialization(gen_cls);

    py::class_<GenInfo>(m, "GenInfo", DocIterator::GenInfo.c_str())
        .def_readonly("id", &GenInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &GenInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &GenInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &GenInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &GenInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &GenInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("is_slack", &GenInfo::is_slack, DocIterator::is_slack.c_str())
        .def_readonly("slack_weight", &GenInfo::slack_weight, DocIterator::slack_weight.c_str())
        .def_readonly("voltage_regulator_on", &GenInfo::voltage_regulator_on, DocIterator::voltage_regulator_on.c_str())
        .def_readonly("target_p_mw", &GenInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_vm_pu", &GenInfo::target_vm_pu, DocIterator::target_vm_pu.c_str())
        .def_readonly("target_q_mvar", &GenInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("min_q_mvar", &GenInfo::min_q_mvar, DocIterator::min_q_mvar.c_str())
        .def_readonly("max_q_mvar", &GenInfo::max_q_mvar, DocIterator::max_q_mvar.c_str())
        .def_readonly("regulated_bus_id", &GenInfo::regulated_bus_id, DocIterator::regulated_bus_id.c_str())
        .def_readonly("has_res", &GenInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &GenInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &GenInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &GenInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &GenInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        .def_readonly("voltage_level_id", &GenInfo::sub_id, DocIterator::sub_id.c_str());

    auto svc_cls = py::class_<SvcContainer>(m, "SvcContainer", "Container of all the Static Var Compensators (SVC) of the grid.")
        .def("__len__", [](const SvcContainer & data) { return data.nb(); })
        .def("__getitem__", [](const SvcContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const SvcContainer & data)  {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>());
    add_pickle(svc_cls, "SvcContainer");
    add_binary_serialization(svc_cls);
    // exposed (also) so TestSerializedEnumValues can pin the integer values,
    // which are serialized verbatim in binary files and pickles
    py::enum_<SvcContainer::RegulationMode>(svc_cls, "RegulationMode", "Regulation mode of an SVC (values follow the IIDM model of powsybl)")
        .value("OFF", SvcContainer::RegulationMode::OFF)
        .value("VOLTAGE", SvcContainer::RegulationMode::VOLTAGE)
        .value("REACTIVE_POWER", SvcContainer::RegulationMode::REACTIVE_POWER);

    py::class_<SvcInfo>(m, "SvcInfo", "Information about one Static Var Compensator (SVC).")
        .def_readonly("id", &SvcInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &SvcInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &SvcInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &SvcInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &SvcInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &SvcInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("regulation_mode", &SvcInfo::regulation_mode, "0 = OFF, 1 = VOLTAGE, 2 = REACTIVE_POWER")
        .def_readonly("target_vm_pu", &SvcInfo::target_vm_pu, "voltage setpoint (pu of the regulated bus), VOLTAGE mode")
        .def_readonly("target_q_mvar", &SvcInfo::target_q_mvar, "reactive setpoint (MVAr), REACTIVE_POWER mode")
        .def_readonly("slope_pu", &SvcInfo::slope_pu, "voltage/reactive slope (pu), 0 = no droop")
        .def_readonly("b_min", &SvcInfo::b_min, "minimum susceptance (stored, NEVER enforced)")
        .def_readonly("b_max", &SvcInfo::b_max, "maximum susceptance (stored, NEVER enforced)")
        .def_readonly("regulated_bus_id", &SvcInfo::regulated_bus_id, "grid bus id whose voltage is regulated (== bus_id for local control)")
        .def_readonly("has_res", &SvcInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &SvcInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &SvcInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &SvcInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &SvcInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        .def_readonly("voltage_level_id", &SvcInfo::sub_id, DocIterator::sub_id.c_str());

    auto sgen_cls = py::class_<SGenContainer>(m, "SGenContainer", DocIterator::SGenContainer.c_str())
        .def("__len__", [](const SGenContainer & data) { return data.nb(); })
        .def("__getitem__", [](const SGenContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const SGenContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def("get_bus_id", &SGenContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>());
    add_pickle(sgen_cls, "SGenContainer");
    add_binary_serialization(sgen_cls);

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
    add_binary_serialization(load_cls);

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

    // storage units share the LoadContainer / LoadInfo documentation (they are PQ
    // injections in the load convention) but are exposed as their own type.
    auto storage_cls = py::class_<StorageContainer>(m, "StorageContainer", DocIterator::LoadContainer.c_str())
        .def("__len__", [](const StorageContainer & data) { return data.nb(); })
        .def("__getitem__", [](const StorageContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const StorageContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def("get_bus_id", &StorageContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>());
    add_pickle(storage_cls, "StorageContainer");
    add_binary_serialization(storage_cls);

    py::class_<StorageInfo>(m, "StorageInfo", DocIterator::LoadInfo.c_str())
        .def_readonly("id", &StorageInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &StorageInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &StorageInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &StorageInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &StorageInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &StorageInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("target_p_mw", &StorageInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_q_mvar", &StorageInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("has_res", &StorageInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &StorageInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &StorageInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &StorageInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &StorageInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        .def_readonly("voltage_level_id", &StorageInfo::sub_id, DocIterator::sub_id.c_str());

    auto shunt_cls = py::class_<ShuntContainer>(m, "ShuntContainer", DocIterator::ShuntContainer.c_str())
        .def("__len__", [](const ShuntContainer & data) { return data.nb(); })
        .def("__getitem__", [](const ShuntContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const ShuntContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def("get_bus_id", &ShuntContainer::get_bus_id_numpy, "TODO doc", py::keep_alive<0, 1>());
    add_pickle(shunt_cls, "ShuntContainer");
    add_binary_serialization(shunt_cls);

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
    add_binary_serialization(trafo_cls);

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
        .def_readonly("limit_a1_ka", &TrafoInfo::limit_a1_ka, "Current limit, hv side, in kA (NaN if not set, see `LSGrid.set_trafo_current_limit_side1`).")
        .def_readonly("limit_a2_ka", &TrafoInfo::limit_a2_ka, "Current limit, lv side, in kA (NaN if not set, see `LSGrid.set_trafo_current_limit_side2`).")
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
    add_binary_serialization(line_cls);

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
        .def_readonly("bus1_id", &LineInfo::bus_1_id, DocIterator::bus_1_id.c_str())
        .def_readonly("bus2_id", &LineInfo::bus_2_id, DocIterator::bus_2_id.c_str())
        .def_readonly("r_pu", &LineInfo::r_pu, DocIterator::r_pu.c_str())
        .def_readonly("x_pu", &LineInfo::x_pu, DocIterator::x_pu.c_str())
        .def_readonly("h1_pu", &LineInfo::h1_pu, DocIterator::h_pu.c_str())
        .def_readonly("h2_pu", &LineInfo::h2_pu, DocIterator::h_pu.c_str())
        .def_readonly("has_res", &LineInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p1_mw", &LineInfo::res_p1_mw, DocIterator::res_p_1_mw.c_str())
        .def_readonly("res_q1_mvar", &LineInfo::res_q1_mvar, DocIterator::res_q_1_mvar.c_str())
        .def_readonly("res_v1_kv", &LineInfo::res_v1_kv, DocIterator::res_v_1_kv.c_str())
        .def_readonly("res_a1_ka", &LineInfo::res_a1_ka, DocIterator::res_a_1_ka.c_str())
        .def_readonly("res_p2_mw", &LineInfo::res_p2_mw, DocIterator::res_p_2_mw.c_str())
        .def_readonly("res_q2_mvar", &LineInfo::res_q2_mvar, DocIterator::res_q_2_mvar.c_str())
        .def_readonly("res_v2_kv", &LineInfo::res_v2_kv, DocIterator::res_v_2_kv.c_str())
        .def_readonly("res_a2_ka", &LineInfo::res_a2_ka, DocIterator::res_a_2_ka.c_str())
        .def_readonly("limit_a1_ka", &LineInfo::limit_a1_ka, "Current limit, origin side, in kA (NaN if not set, see `LSGrid.set_line_current_limit_side1`).")
        .def_readonly("limit_a2_ka", &LineInfo::limit_a2_ka, "Current limit, extremity side, in kA (NaN if not set, see `LSGrid.set_line_current_limit_side2`).")
        .def_readonly("res_theta1_deg", &LineInfo::res_theta1_deg, DocIterator::res_theta_1_deg.c_str())
        .def_readonly("res_theta2_deg", &LineInfo::res_theta2_deg, DocIterator::res_theta_2_deg.c_str())
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

    auto cs_info_cls = py::class_<ConverterStationInfo>(m, "ConverterStationInfo", DocIterator::ConverterStationInfo.c_str());
    // exposed (also) so TestSerializedEnumValues can pin the integer values,
    // which are serialized verbatim in binary files and pickles
    py::enum_<ConverterStationContainer::ConverterType>(cs_info_cls, "ConverterType", "Type of an hvdc converter station")
        .value("VSC", ConverterStationContainer::ConverterType::VSC)
        .value("LCC", ConverterStationContainer::ConverterType::LCC);
    cs_info_cls
        .def_readonly("id", &ConverterStationInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &ConverterStationInfo::name, DocIterator::name.c_str())
        .def_readonly("sub_id", &ConverterStationInfo::sub_id, DocIterator::sub_id.c_str())
        .def_readonly("pos_topo_vect", &ConverterStationInfo::pos_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected", &ConverterStationInfo::connected, DocIterator::connected.c_str())
        .def_readonly("bus_id", &ConverterStationInfo::bus_id, DocIterator::bus_id.c_str())
        .def_readonly("converter_type", &ConverterStationInfo::converter_type, DocIterator::converter_type.c_str())
        .def_readonly("loss_factor", &ConverterStationInfo::loss_factor, DocIterator::loss_factor.c_str())
        .def_readonly("voltage_regulator_on", &ConverterStationInfo::voltage_regulator_on, DocIterator::voltage_regulator_on.c_str())
        .def_readonly("target_p_mw", &ConverterStationInfo::target_p_mw, DocIterator::target_p_mw.c_str())
        .def_readonly("target_vm_pu", &ConverterStationInfo::target_vm_pu, DocIterator::target_vm_pu.c_str())
        .def_readonly("target_q_mvar", &ConverterStationInfo::target_q_mvar, DocIterator::target_q_mvar.c_str())
        .def_readonly("min_q_mvar", &ConverterStationInfo::min_q_mvar, DocIterator::min_q_mvar.c_str())
        .def_readonly("max_q_mvar", &ConverterStationInfo::max_q_mvar, DocIterator::max_q_mvar.c_str())
        .def_readonly("power_factor", &ConverterStationInfo::power_factor, DocIterator::power_factor.c_str())
        .def_readonly("has_res", &ConverterStationInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p_mw", &ConverterStationInfo::res_p_mw, DocIterator::res_p_mw.c_str())
        .def_readonly("res_q_mvar", &ConverterStationInfo::res_q_mvar, DocIterator::res_q_mvar.c_str())
        .def_readonly("res_theta_deg", &ConverterStationInfo::res_theta_deg, DocIterator::res_theta_deg.c_str())
        .def_readonly("res_v_kv", &ConverterStationInfo::res_v_kv, DocIterator::res_v_kv.c_str())
        .def_readonly("voltage_level_id", &ConverterStationInfo::sub_id, DocIterator::sub_id.c_str());

    auto dcline_cls = py::class_<HvdcLineContainer>(m, "HvdcLineContainer", DocIterator::DCLineContainer.c_str())
        .def("__len__", [](const HvdcLineContainer & data) { return data.nb(); })
        .def("__getitem__", [](const HvdcLineContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const HvdcLineContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>())
        .def("get_bus_id_side_1", &HvdcLineContainer::get_bus_id_side_1_numpy)
        .def("get_bus_id_side_2", &HvdcLineContainer::get_bus_id_side_2_numpy);
    add_pickle(dcline_cls, "HvdcLineContainer");
    add_binary_serialization(dcline_cls);
    // exposed (also) so TestSerializedEnumValues can pin the integer values,
    // which are serialized verbatim in binary files and pickles
    py::enum_<HvdcLineContainer::ConvertersMode>(dcline_cls, "ConvertersMode", "Which side of an hvdc line is the rectifier (the other being the inverter)")
        .value("SIDE_1_RECTIFIER", HvdcLineContainer::ConvertersMode::SIDE_1_RECTIFIER)
        .value("SIDE_2_RECTIFIER", HvdcLineContainer::ConvertersMode::SIDE_2_RECTIFIER);

    py::class_<HvdcLineInfo>(m, "HvdcLineInfo", DocIterator::DCLineInfo.c_str())
        .def_readonly("id", &HvdcLineInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &HvdcLineInfo::name, DocIterator::name.c_str())
        .def_readonly("sub1_id", &HvdcLineInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("sub2_id", &HvdcLineInfo::sub_2_id, DocIterator::sub_id.c_str())
        .def_readonly("pos1_topo_vect", &HvdcLineInfo::pos_1_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("pos2_topo_vect", &HvdcLineInfo::pos_2_topo_vect, DocIterator::pos_topo_vect.c_str())
        .def_readonly("connected_global", &HvdcLineInfo::connected_global, DocIterator::connected.c_str())
        .def_readonly("connected1", &HvdcLineInfo::connected_1, DocIterator::connected.c_str())
        .def_readonly("connected2", &HvdcLineInfo::connected_2, DocIterator::connected.c_str())
        .def_readonly("bus1_id", &HvdcLineInfo::bus_1_id, DocIterator::bus_1_id.c_str())
        .def_readonly("bus2_id", &HvdcLineInfo::bus_2_id, DocIterator::bus_2_id.c_str())
        .def_readonly("target_p1_mw", &HvdcLineInfo::target_p_1_mw, DocIterator::target_p_1_mw_dcline.c_str())
        .def_readonly("p2_mw", &HvdcLineInfo::p_2_mw, DocIterator::target_p_2_mw_dcline.c_str())
        .def_readonly("target_vm1_pu", &HvdcLineInfo::target_vm_1_pu, DocIterator::target_vm_1_pu_dcline.c_str())
        .def_readonly("target_vm2_pu", &HvdcLineInfo::target_vm_2_pu, DocIterator::target_vm_2_pu_dcline.c_str())
        .def_readonly("loss_pct", &HvdcLineInfo::loss_pct, DocIterator::loss_pct.c_str())
        .def_readonly("loss_mw", &HvdcLineInfo::loss_mw, DocIterator::loss_mw.c_str())
        .def_readonly("converters_mode", &HvdcLineInfo::converters_mode, DocIterator::converters_mode.c_str())
        .def_readonly("p_setpoint_mw", &HvdcLineInfo::p_setpoint_mw, DocIterator::p_setpoint_mw.c_str())
        .def_readonly("r_ohm", &HvdcLineInfo::r_ohm, DocIterator::r_ohm.c_str())
        .def_readonly("nominal_v_kv", &HvdcLineInfo::nominal_v_kv, DocIterator::nominal_v_kv.c_str())
        .def_readonly("droop_enabled", &HvdcLineInfo::droop_enabled, DocIterator::droop_enabled.c_str())
        .def_readonly("droop_p0_mw", &HvdcLineInfo::droop_p0_mw, DocIterator::droop_p0_mw.c_str())
        .def_readonly("droop_k_mw_per_rad", &HvdcLineInfo::droop_k_mw_per_rad, DocIterator::droop_k_mw_per_rad.c_str())
        .def_readonly("pmax_1to2_mw", &HvdcLineInfo::pmax_1to2_mw, DocIterator::pmax_1to2_mw.c_str())
        .def_readonly("pmax_2to1_mw", &HvdcLineInfo::pmax_2to1_mw, DocIterator::pmax_2to1_mw.c_str())
        .def_readonly("status_droop", &HvdcLineInfo::status_droop, DocIterator::status_droop.c_str())
        .def_readonly("station1", &HvdcLineInfo::station_side_1, DocIterator::station_side_1.c_str())
        .def_readonly("station2", &HvdcLineInfo::station_side_2, DocIterator::station_side_2.c_str())
        .def_readonly("has_res", &HvdcLineInfo::has_res, DocIterator::has_res.c_str())
        .def_readonly("res_p1_mw", &HvdcLineInfo::res_p1_mw, DocIterator::res_p_1_mw_dcline.c_str())
        .def_readonly("res_p2_mw", &HvdcLineInfo::res_p2_mw, DocIterator::res_p_2_mw_dcline.c_str())
        .def_readonly("res_q1_mvar", &HvdcLineInfo::res_q1_mvar, DocIterator::res_q_1_mvar_dcline.c_str())
        .def_readonly("res_q2_mvar", &HvdcLineInfo::res_q2_mvar, DocIterator::res_q_2_mvar_dcline.c_str())
        .def_readonly("res_v1_kv", &HvdcLineInfo::res_v1_kv, DocIterator::res_v_1_kv_dcline.c_str())
        .def_readonly("res_v2_kv", &HvdcLineInfo::res_v2_kv, DocIterator::res_v_2_kv_dcline.c_str())
        .def_readonly("res_theta1_deg", &HvdcLineInfo::res_theta1_deg, DocIterator::res_theta_1_deg_dcline.c_str())
        .def_readonly("res_theta2_deg", &HvdcLineInfo::res_theta2_deg, DocIterator::res_theta_2_deg_dcline.c_str())
        .def_readonly("voltage_level1_id", &HvdcLineInfo::sub_1_id, DocIterator::sub_id.c_str())
        .def_readonly("voltage_level2_id", &HvdcLineInfo::sub_2_id, DocIterator::sub_id.c_str());

    auto sub_cls = py::class_<SubstationContainer>(m, "SubstationContainer", "TODO")
        .def("__len__", [](const SubstationContainer & data) { return data.nb(); })
        .def("__getitem__", [](const SubstationContainer & data, int k){return data[k]; } )
        .def("__iter__", [](const SubstationContainer & data) {
                return py::make_iterator(data.begin(), data.end());
            }, py::keep_alive<0, 1>());
    add_pickle(sub_cls, "SubstationContainer");
    add_binary_serialization(sub_cls);

    py::class_<SubstationInfo>(m, "SubstationInfo", "TODO")
        .def_readonly("id", &SubstationInfo::id, DocIterator::id.c_str())
        .def_readonly("name", &SubstationInfo::name, DocIterator::name.c_str())
        .def_readonly("nb_max_busbars", &SubstationInfo::nb_max_busbars, DocIterator::nb_max_busbars.c_str())
        .def_readonly("vn_kv", &SubstationInfo::vn_kv, DocIterator::vn_kv.c_str());
}
