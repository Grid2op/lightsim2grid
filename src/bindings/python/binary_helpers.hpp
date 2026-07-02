// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BINARY_HELPERS_HPP
#define BINARY_HELPERS_HPP

#include <pybind11/pybind11.h>
#include "BinaryArchive.hpp"

namespace py = pybind11;

// Helper: attach save_binary()/load_binary() python methods to any container
// that exposes get_state()/set_state() and a nested StateRes type (the same
// contract used by add_pickle in pickle_helpers.hpp). This is an additive,
// faster alternative to pickle -- NOT a replacement, and NOT cross-version
// compatible (a version mismatch is a hard failure, see BinaryArchive.hpp).
//
// These lambdas call ls2g::save_binary_generic/load_binary_generic directly
// (rather than T::save_binary/T::load_binary): VERSION_MAJOR/MEDIUM/MINOR are
// only defined via target_compile_definitions on the python bindings target,
// not on lightsim2grid_core, so this is the single, consistent translation
// unit where those macros carry the real version (mirrors pickle_helpers.hpp).
template<typename T>
void add_binary_serialization(py::class_<T>& cls) {
    cls.def("save_binary", [](const T& obj, const std::string& path) {
        ls2g::save_binary_generic(obj, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
    }, py::arg("path"),
       "Save this object's state to a fast custom binary file (additive alternative "
       "to pickle: faster, but NOT portable across lightsim2grid versions).");
    cls.def_static("load_binary", [](const std::string& path) {
        return ls2g::load_binary_generic<T>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
    }, py::arg("path"),
       "Load an object previously saved with save_binary(). Raises RuntimeError on "
       "a version mismatch or a corrupted / truncated file.");
}

#endif // BINARY_HELPERS_HPP
