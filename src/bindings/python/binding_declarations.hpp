// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BINDING_DECLARATIONS_HPP
#define BINDING_DECLARATIONS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void bind_enums(py::module_& m);
void bind_solvers(py::module_& m);
void bind_containers(py::module_& m);
void bind_misc(py::module_& m);
void bind_gridmodel(py::module_& m);
void bind_batch(py::module_& m);

#endif // BINDING_DECLARATIONS_HPP
