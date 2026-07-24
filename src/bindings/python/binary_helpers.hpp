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
// that exposes get_state()/set_state(), a nested StateRes type (the same
// contract used by add_pickle in pickle_helpers.hpp) and binary_type_tag().
// This is an additive, faster alternative to pickle -- NOT a replacement.
// Files are readable by any lightsim2grid version sharing the same
// BINARY_FORMAT_VERSION (see BinaryArchive.hpp), and ill-formed input
// (garbage, truncation, corrupted sizes, wrong object type, trailing bytes)
// raises RuntimeError instead of crashing or over-allocating. Loading a whole
// grid additionally runs LSGrid::check_grid() (via set_state), so a byte-wise
// well-formed but inconsistent grid -- an out-of-range bus / substation /
// topology-vector index, or a non-finite value -- is rejected with an
// IndexError (out-of-range index) or a RuntimeError (structural inconsistency).
//
// These lambdas call ls2g::save_binary_generic/load_binary_generic directly
// (rather than T::save_binary/T::load_binary): VERSION_MAJOR/MEDIUM/MINOR are
// only defined via target_compile_definitions on the python bindings target,
// not on lightsim2grid_core, so this is the single, consistent translation
// unit where those macros carry the real version (mirrors pickle_helpers.hpp).
template<typename T>
void add_binary_serialization(py::class_<T>& cls) {
    cls.def("save_binary", [](const T& obj, const std::string& path, bool atomic) {
        ls2g::save_binary_generic(obj, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
    }, py::arg("path"), py::arg("atomic") = true,
       "Save this object's state to a fast custom binary file (additive alternative "
       "to pickle). By default (atomic=True) the write is atomic: an existing file "
       "at that path is only replaced once the new content has been written "
       "completely (an interrupted save never destroys a previous file). Pass "
       "atomic=False to write the destination directly instead -- marginally faster "
       "(skips one temporary file + rename), without that protection. The file stays "
       "readable by any lightsim2grid version sharing the same binary format number.");
    cls.def_static("load_binary", [](const std::string& path) {
        return ls2g::load_binary_generic<T>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
    }, py::arg("path"),
       "Load an object previously saved with save_binary(). Raises RuntimeError on "
       "an incompatible binary format, a wrong object type, or a corrupted / "
       "truncated file (including corrupted internal sizes: no attempt is made to "
       "allocate more data than the file actually contains). Loading a whole grid "
       "additionally validates its consistency (see check_grid): a byte-wise "
       "well-formed but inconsistent grid raises IndexError (out-of-range index) "
       "or RuntimeError (structural inconsistency).");
}

#endif // BINARY_HELPERS_HPP
