// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// The state of LSGrid, for python: pickle and the binary format -- one of the
// translation units the LSGrid bindings are spread over, see the note at the top
// of binding_lsgrid.cpp.
//
// Why the pickle is not the generic add_pickle() every container uses. That
// helper casts the whole StateRes tuple through pybind11, and LSGrid::StateRes
// nests every container's StateRes (ten tuples of a dozen vectors each, plus a
// vector of tuples for the detailed topology): instantiating those casters in
// both directions cost GCC ~2.2 GB at compile time, in one translation unit,
// which is what made the CI containers kill the build. The containers' own
// casters already exist, compiled once in binding_containers.cpp, and are
// reachable at run time as their __getstate__ / __setstate__: the LSGrid pickle
// below builds the same nested tuple by calling them, and rebuilds the C++ state
// the way pickle.loads itself does (`cls.__new__(cls)` then `__setstate__`),
// which is also what test_state_poisoning.py does by hand. This unit only casts
// registered class objects, strings and small vectors: ~0.9 GB.
//
// The pickled layout is EXACTLY the one add_pickle() produced --
// `(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, state)` with `state` the
// LSGrid::StateRes tuple, every container slot its own nested StateRes tuple --
// and tests poke inside it (test_state_poisoning.py, test_check_grid.py), so keep
// it that way.

#include "binding_declarations.hpp"
#include "binary_helpers.hpp"
#include "LSGrid.hpp"

using namespace ls2g;

namespace {

// the nested StateRes tuple of a container, as its own pickle emits it
// (element [3] of `(major, medium, minor, state)`)
template<class Container>
py::object container_state(const Container & container)
{
    py::object obj = py::cast(container);  // a copy, pickled right away
    py::tuple full = obj.attr("__getstate__")();
    return full[3];
}

// the C++ StateRes of a container, from the nested tuple its pickle emits:
// through the container's own __setstate__, the way pickle.loads restores it
template<class Container>
typename Container::StateRes container_state_from_py(py::handle inner)
{
    py::object cls = py::type::of<Container>();
    py::object obj = cls.attr("__new__")(cls);
    obj.attr("__setstate__")(py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, inner));
    return obj.cast<Container &>().get_state();
}

void check_version(const py::tuple & py_state, const char * class_name)
{
    if (py_state.size() != 4) {
        throw std::runtime_error(std::string("Invalid state size when loading ") + class_name);
    }
    if (py_state[0].cast<std::string>() != VERSION_MAJOR)
        throw std::runtime_error(std::string("Invalid state size when loading ") + class_name +
            ": wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
    if (py_state[1].cast<std::string>() != VERSION_MEDIUM)
        throw std::runtime_error(std::string("Invalid state size when loading ") + class_name +
            ": wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
    if (py_state[2].cast<std::string>() != VERSION_MINOR)
        throw std::runtime_error(std::string("Invalid state size when loading ") + class_name +
            ": wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
}

using AlgoConfigState = std::tuple<std::vector<int>, std::vector<double> >;

}  // anonymous namespace

void bind_lsgrid_state(py::class_<LSGrid> & cls)
{
    cls.def(py::pickle(
        [](const LSGrid & grid) {
            const LSGrid::StateRes st = grid.get_state();
            py::tuple state = py::make_tuple(
                std::get<LSGrid::VERSION_MAJOR_ID>(st),
                std::get<LSGrid::VERSION_MEDIUM_ID>(st),
                std::get<LSGrid::VERSION_MINOR_ID>(st),
                std::get<LSGrid::LS_TO_ORIG_ID>(st),
                std::get<LSGrid::INIT_VM_PU_ID>(st),
                std::get<LSGrid::SN_MVA_ID>(st),
                container_state(grid.get_substations()),
                container_state(grid.get_lines()),
                container_state(grid.get_shunts()),
                container_state(grid.get_trafos()),
                container_state(grid.get_generators()),
                container_state(grid.get_loads()),
                container_state(grid.get_static_generators()),
                container_state(grid.get_storages()),
                container_state(grid.get_dclines()),
                container_state(grid.get_svcs()),
                std::get<LSGrid::AC_ALGO_NAME_ID>(st),
                std::get<LSGrid::DC_ALGO_NAME_ID>(st),
                std::get<LSGrid::AC_ALGO_CONFIG_ID>(st),
                std::get<LSGrid::DC_ALGO_CONFIG_ID>(st),
                std::get<LSGrid::INIT_KWARGS_KEYS_ID>(st),
                std::get<LSGrid::INIT_KWARGS_VALUES_ID>(st),
                std::get<LSGrid::BUS_FUSION_REP_ID>(st));
            return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, state);
        },
        [](py::tuple py_state) {
            check_version(py_state, "LSGrid");
            py::tuple state = py_state[3].cast<py::tuple>();
            if (state.size() != std::tuple_size<LSGrid::StateRes>::value) {
                throw std::runtime_error("Invalid state size when loading LSGrid: the inner state has " +
                                         std::to_string(state.size()) + " entries, " +
                                         std::to_string(std::tuple_size<LSGrid::StateRes>::value) + " expected.");
            }
            LSGrid::StateRes st;
            std::get<LSGrid::VERSION_MAJOR_ID>(st) = state[LSGrid::VERSION_MAJOR_ID].cast<std::string>();
            std::get<LSGrid::VERSION_MEDIUM_ID>(st) = state[LSGrid::VERSION_MEDIUM_ID].cast<std::string>();
            std::get<LSGrid::VERSION_MINOR_ID>(st) = state[LSGrid::VERSION_MINOR_ID].cast<std::string>();
            std::get<LSGrid::LS_TO_ORIG_ID>(st) = state[LSGrid::LS_TO_ORIG_ID].cast<std::vector<int> >();
            std::get<LSGrid::INIT_VM_PU_ID>(st) = state[LSGrid::INIT_VM_PU_ID].cast<real_type>();
            std::get<LSGrid::SN_MVA_ID>(st) = state[LSGrid::SN_MVA_ID].cast<real_type>();
            std::get<LSGrid::SUBSTATION_ID>(st) = container_state_from_py<SubstationContainer>(state[LSGrid::SUBSTATION_ID]);
            std::get<LSGrid::LINE_ID>(st) = container_state_from_py<LineContainer>(state[LSGrid::LINE_ID]);
            std::get<LSGrid::SHUNT_ID>(st) = container_state_from_py<ShuntContainer>(state[LSGrid::SHUNT_ID]);
            std::get<LSGrid::TRAFO_ID>(st) = container_state_from_py<TrafoContainer>(state[LSGrid::TRAFO_ID]);
            std::get<LSGrid::GEN_ID>(st) = container_state_from_py<GeneratorContainer>(state[LSGrid::GEN_ID]);
            std::get<LSGrid::LOAD_ID>(st) = container_state_from_py<LoadContainer>(state[LSGrid::LOAD_ID]);
            std::get<LSGrid::SGEN_ID>(st) = container_state_from_py<SGenContainer>(state[LSGrid::SGEN_ID]);
            std::get<LSGrid::STORAGE_ID>(st) = container_state_from_py<StorageContainer>(state[LSGrid::STORAGE_ID]);
            std::get<LSGrid::HVDC_ID>(st) = container_state_from_py<HvdcLineContainer>(state[LSGrid::HVDC_ID]);
            std::get<LSGrid::SVC_ID>(st) = container_state_from_py<SvcContainer>(state[LSGrid::SVC_ID]);
            std::get<LSGrid::AC_ALGO_NAME_ID>(st) = state[LSGrid::AC_ALGO_NAME_ID].cast<std::string>();
            std::get<LSGrid::DC_ALGO_NAME_ID>(st) = state[LSGrid::DC_ALGO_NAME_ID].cast<std::string>();
            std::get<LSGrid::AC_ALGO_CONFIG_ID>(st) = state[LSGrid::AC_ALGO_CONFIG_ID].cast<AlgoConfigState>();
            std::get<LSGrid::DC_ALGO_CONFIG_ID>(st) = state[LSGrid::DC_ALGO_CONFIG_ID].cast<AlgoConfigState>();
            std::get<LSGrid::INIT_KWARGS_KEYS_ID>(st) = state[LSGrid::INIT_KWARGS_KEYS_ID].cast<std::vector<std::string> >();
            std::get<LSGrid::INIT_KWARGS_VALUES_ID>(st) = state[LSGrid::INIT_KWARGS_VALUES_ID].cast<std::vector<std::string> >();
            std::get<LSGrid::BUS_FUSION_REP_ID>(st) = state[LSGrid::BUS_FUSION_REP_ID].cast<std::vector<int> >();
            LSGrid res;
            res.set_state(st);
            return res;
        }
    ));

    add_binary_serialization(cls);
    // Companion to load_binary(): loads the grid *data* without re-selecting the
    // solver it was saved with. Lets a grid saved with a solver that is not
    // available here (a plugin that has not been loaded, or a built-in needing an
    // optional backend this build lacks) still be loaded -- the grid keeps the
    // default solvers, and you pick one yourself with change_algorithm().
    cls.def_static("load_binary_without_algorithm", [](const std::string& path) {
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
}
