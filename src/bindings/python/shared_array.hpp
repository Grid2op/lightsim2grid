// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LS2G_BINDINGS_SHARED_ARRAY_H
#define LS2G_BINDINGS_SHARED_ARRAY_H

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace ls2g_bindings {

/**
 * A read-only numpy array over a matrix a SharedMat holds (see
 * core/batch_algorithm/SharedMat.hpp), without a copy. The array's base is a capsule
 * owning a copy of the shared_ptr, so the buffer lives as long as the array does --
 * whatever the batch object does next (compute, clear, or being deleted).
 *
 * Read-only as the views of a const reference were: the algorithm may still reuse
 * this buffer in place once no array holds it, so writing through it is not a way
 * to change a result.
 */
template<class M>
pybind11::array shared_mat_to_array(std::shared_ptr<const M> mat)
{
    namespace py = pybind11;
    using Scalar = typename M::Scalar;
    const M & m = *mat;
    std::unique_ptr<std::shared_ptr<const M> > owner(new std::shared_ptr<const M>(std::move(mat)));
    py::capsule base(owner.get(), [](void * p){ delete static_cast<std::shared_ptr<const M> *>(p); });
    owner.release();  // the capsule deletes it from now on

    const py::ssize_t item = static_cast<py::ssize_t>(sizeof(Scalar));
    const py::ssize_t outer = static_cast<py::ssize_t>(m.outerStride()) * item;
    const py::ssize_t inner = static_cast<py::ssize_t>(m.innerStride()) * item;
    py::array res(py::dtype::of<Scalar>(),
                  {static_cast<py::ssize_t>(m.rows()), static_cast<py::ssize_t>(m.cols())},
                  {M::IsRowMajor ? outer : inner, M::IsRowMajor ? inner : outer},
                  m.data(), base);
    res.attr("setflags")(py::arg("write") = false);
    return res;
}

}  // namespace ls2g_bindings

#endif  // LS2G_BINDINGS_SHARED_ARRAY_H
