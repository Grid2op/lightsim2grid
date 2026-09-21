// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


// The topology side of LSGrid: the grid2op vectors, the substation ids, the switches -- one of the four translation units the LSGrid bindings are spread
// over, see the note at the top of binding_lsgrid.cpp.

#include "binding_declarations.hpp"
#include "LSGrid.hpp"
#include "help_fun_msg.hpp"

using namespace ls2g;

void bind_lsgrid_topology(py::class_<LSGrid> & cls) {
    cls
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
        .def("set_sgen_to_subid", &LSGrid::set_sgen_to_subid, DocLSGrid::set_sgen_to_subid.c_str())
        .def("set_svc_to_subid", &LSGrid::set_svc_to_subid, DocLSGrid::set_svc_to_subid.c_str())
        .def("set_dcline_to_sub1_id", &LSGrid::set_dcline_to_sub1_id, DocLSGrid::set_dcline_to_sub1_id.c_str())
        .def("set_dcline_to_sub2_id", &LSGrid::set_dcline_to_sub2_id, DocLSGrid::set_dcline_to_sub2_id.c_str())

        // detailed topology (switches inside each substation)
        .def("init_detailed_topology", &LSGrid::init_detailed_topology, DocLSGrid::init_detailed_topology.c_str())
        .def("has_detailed_topology", &LSGrid::has_detailed_topology, DocLSGrid::has_detailed_topology.c_str())
        .def("set_switch_names", &LSGrid::set_switch_names, DocLSGrid::set_switch_names.c_str())
        .def("set_busbar_section_names", &LSGrid::set_busbar_section_names, DocLSGrid::set_busbar_section_names.c_str())
        .def("set_load_to_node_id", &LSGrid::set_load_to_node_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_gen_to_node_id", &LSGrid::set_gen_to_node_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_sgen_to_node_id", &LSGrid::set_sgen_to_node_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_storage_to_node_id", &LSGrid::set_storage_to_node_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_shunt_to_node_id", &LSGrid::set_shunt_to_node_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_svc_to_node_id", &LSGrid::set_svc_to_node_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_line_to_node1_id", &LSGrid::set_line_to_node1_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_line_to_node2_id", &LSGrid::set_line_to_node2_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_trafo_to_node1_id", &LSGrid::set_trafo_to_node1_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_trafo_to_node2_id", &LSGrid::set_trafo_to_node2_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_dcline_to_node1_id", &LSGrid::set_dcline_to_node1_id, DocLSGrid::set_to_node_id.c_str())
        .def("set_dcline_to_node2_id", &LSGrid::set_dcline_to_node2_id, DocLSGrid::set_to_node_id.c_str())
        .def("project_switches", &LSGrid::project_switches, DocLSGrid::project_switches.c_str())
        .def("set_switch_open", &LSGrid::set_switch_open, DocLSGrid::set_switch_open.c_str())
        .def("update_switches", &LSGrid::update_switches, DocLSGrid::update_switches.c_str())
        // the two containers are views holding a pointer into this grid: keep it alive
        .def("get_switches", &LSGrid::get_switches, DocLSGrid::get_switches.c_str(), py::keep_alive<0, 1>())
        .def("get_busbar_sections", &LSGrid::get_busbar_sections, DocLSGrid::get_busbar_sections.c_str(), py::keep_alive<0, 1>())
        .def("get_node_bus", &LSGrid::get_node_bus, DocLSGrid::get_node_bus.c_str())
        .def("get_substation_topology", &LSGrid::get_substation_topology, DocLSGrid::get_substation_topology.c_str(), py::return_value_policy::reference_internal)
;
}
