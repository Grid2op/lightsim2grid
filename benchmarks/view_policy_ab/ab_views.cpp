// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// A/B of the return value policy of LSGrid's result getters (docs/dev_notes/result_view_lifetime.md):
// the very same C++ getters bound three times, with `reference` (what binding_lsgrid.cpp does),
// `reference_internal` (the proposed fix: each array keeps the grid alive) and `copy` (the other
// way out, for scale). Everything else is identical: same pybind11, same casters, same object.
// Built and driven by run_ab.py.

#include <memory>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "LSGrid.hpp"

namespace py = pybind11;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;

namespace {

// n buses on a ring with chords every 7 buses, a small load on each bus, a generator on every
// 5th one (the first the slack): every result array has a size proportional to n, and the
// solve converges at any n
std::unique_ptr<LSGrid> make_grid(int n_bus)
{
    if (n_bus < 10) throw std::runtime_error("make_grid: at least 10 buses");
    std::unique_ptr<LSGrid> grid(new LSGrid());
    grid->set_sn_mva(100.);
    grid->set_init_vm_pu(1.);
    const int n_line = 2 * n_bus;
    grid->init_bus(static_cast<unsigned int>(n_bus), 1, RealVect::Constant(n_bus, 138.), n_line, 0);
    Eigen::VectorXi from(n_line), to(n_line);
    for (int i = 0; i < n_bus; ++i) {
        from(i) = i;
        to(i) = (i + 1) % n_bus;
        from(n_bus + i) = i;
        to(n_bus + i) = (i + 7) % n_bus;
    }
    grid->init_powerlines(RealVect::Constant(n_line, 0.01), RealVect::Constant(n_line, 0.1),
                          CplxVect::Zero(n_line), from, to);
    Eigen::VectorXi load_bus(n_bus);
    for (int i = 0; i < n_bus; ++i) load_bus(i) = i;
    grid->init_loads(RealVect::Constant(n_bus, 1.), RealVect::Constant(n_bus, .3), load_bus);
    const int n_gen = n_bus / 5;
    Eigen::VectorXi gen_bus(n_gen);
    for (int g = 0; g < n_gen; ++g) gen_bus(g) = 5 * g;
    grid->init_generators(RealVect::Constant(n_gen, 5.), RealVect::Constant(n_gen, 1.02),
                          RealVect::Constant(n_gen, -1000.), RealVect::Constant(n_gen, 1000.), gen_bus);
    grid->add_gen_slackbus(0, 1.);
    const CplxVect V = grid->ac_pf(CplxVect::Constant(n_bus, {1., 0.}), 30, 1e-8);
    if (V.size() != n_bus) throw std::runtime_error("make_grid: the powerflow did not converge");
    return grid;
}

template<class Getter>
void bind_three(py::class_<LSGrid> & cls, const char * name, Getter getter)
{
    const std::string base(name);
    cls.def((base + "_reference").c_str(), getter, py::return_value_policy::reference);
    cls.def((base + "_reference_internal").c_str(), getter, py::return_value_policy::reference_internal);
    cls.def((base + "_copy").c_str(), getter, py::return_value_policy::copy);
}

}  // namespace

PYBIND11_MODULE(ab_views, m)
{
    // module_local: lightsim2grid's own LSGrid binding may live in the same process
    py::class_<LSGrid> cls(m, "Grid", py::module_local());
    bind_three(cls, "get_gen_target_p", &LSGrid::get_gen_target_p);  // one array
    bind_three(cls, "get_gen_res", &LSGrid::get_gen_res);            // tuple of 3
    bind_three(cls, "get_line_res1", &LSGrid::get_line_res1);        // tuple of 4
    bind_three(cls, "get_line_res1_full", &LSGrid::get_line_res1_full);  // tuple of 5
    cls.def("nop", [](const LSGrid &) {});  // the dispatch alone, to subtract
    m.def("make_grid", &make_grid, py::arg("n_bus"));
}
