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
#include <tuple>
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// A state is a tuple of the states of every container, itself nested tuples. Handing the whole
// std::tuple to pybind11 in one go instantiates a single enormous tuple_caster: it costs several
// GB of compiler memory (OOM-kills gcc 8) and overflows MSVC's 64K type-record limit (C1067).
// Converting it one top-level element at a time gives the very same python tuple from many
// small casters. (The dummy array stands in for a C++17 fold expression.)
namespace pickle_detail {
template<typename T> struct is_tuple : std::false_type {};
template<typename... Ts> struct is_tuple<std::tuple<Ts...>> : std::true_type {};

// leaf: whatever pybind11 knows how to convert (vectors, Eigen arrays, scalars, strings, ...)
template<typename T>
typename std::enable_if<!is_tuple<T>::value, py::object>::type to_py(const T & v) { return py::cast(v); }
template<typename T>
typename std::enable_if<!is_tuple<T>::value, T>::type from_py(const py::handle & h) { return h.cast<T>(); }

// tuple: recurse, one element at a time
template<typename Tuple, std::size_t... I>
py::object tuple_to_py(const Tuple & t, std::index_sequence<I...>);
template<typename Tuple, std::size_t... I>
Tuple tuple_from_py(const py::tuple & h, std::index_sequence<I...>);

template<typename T>
typename std::enable_if<is_tuple<T>::value, py::object>::type to_py(const T & t) {
    return tuple_to_py(t, std::make_index_sequence<std::tuple_size<T>::value>{});
}
template<typename T>
typename std::enable_if<is_tuple<T>::value, T>::type from_py(const py::handle & h) {
    constexpr std::size_t n = std::tuple_size<T>::value;
    py::tuple py_t = py::reinterpret_borrow<py::tuple>(h);
    if (!PyTuple_Check(h.ptr()) || py_t.size() != n) throw std::runtime_error("Invalid state size when loading");
    return tuple_from_py<T>(py_t, std::make_index_sequence<n>{});
}

template<typename Tuple, std::size_t... I>
py::object tuple_to_py(const Tuple & t, std::index_sequence<I...>) {
    py::tuple res(sizeof...(I));
    int dummy[] = {0, (PyTuple_SetItem(res.ptr(), I, to_py(std::get<I>(t)).release().ptr()), 0)...};
    (void)dummy;
    return std::move(res);
}
template<typename Tuple, std::size_t... I>
Tuple tuple_from_py(const py::tuple & h, std::index_sequence<I...>) {
    // braced init guarantees left-to-right evaluation
    return Tuple{from_py<typename std::tuple_element<I, Tuple>::type>(h[I])...};
}
}  // namespace pickle_detail

// Helper: attach __getstate__/__setstate__ pickle support to any container
// that exposes get_state()/set_state() and a nested StateRes type.
template<typename T>
void add_pickle(py::class_<T>& cls, const char* class_name) {
    cls.def(py::pickle(
        [](const T& obj) {
            return py::make_tuple(VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, pickle_detail::to_py(obj.get_state()));
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
            auto state = pickle_detail::from_py<typename T::StateRes>(py_state[3]);
            res.set_state(state);
            return res;
        }
    ));
}

#endif // PICKLE_HELPERS_HPP
