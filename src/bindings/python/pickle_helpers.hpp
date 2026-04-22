// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef PICKLE_HELPERS_HPP
#define PICKLE_HELPERS_HPP

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// Helper: attach __getstate__/__setstate__ pickle support to any container
// that exposes get_state()/set_state() and a nested StateRes type.
template<typename T>
void add_pickle(py::class_<T>& cls, const char* class_name) {
    cls.def(py::pickle(
        [](const T& obj) {
            return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, obj.get_state());
        },
        [class_name](py::tuple py_state) {
            if (py_state.size() != 4) {
                std::cout << class_name << ".__setstate__ : state size " << py_state.size() << std::endl;
                throw std::runtime_error(std::string("Invalid state size when loading ") + class_name);
            }
            T res{};
            std::string major = py_state[0].cast<std::string>();
            if (major != VERSION_MAJOR)
                throw std::runtime_error(std::string("Invalid state size when loading ") + class_name +
                    ": wrong lightsim2grid MAJOR.minor.patch version (you can only load pickle from same lightsim2grid version)");
            std::string minor = py_state[1].cast<std::string>();
            if (minor != VERSION_MEDIUM)
                throw std::runtime_error(std::string("Invalid state size when loading ") + class_name +
                    ": wrong lightsim2grid major.MINOR.patch version (you can only load pickle from same lightsim2grid version)");
            std::string patch = py_state[2].cast<std::string>();
            if (patch != VERSION_MINOR)
                throw std::runtime_error(std::string("Invalid state size when loading ") + class_name +
                    ": wrong lightsim2grid major.minor.PATCH version (you can only load pickle from same lightsim2grid version)");
            auto state = py_state[3].cast<typename T::StateRes>();
            res.set_state(state);
            return res;
        }
    ));
}

#endif // PICKLE_HELPERS_HPP
