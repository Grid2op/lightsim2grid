// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SPARSE_HELPERS_HPP
#define SPARSE_HELPERS_HPP

#include <memory>
#include <utility>
#include <Eigen/Sparse>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// Eigen sparse matrices to scipy, without copying their arrays.
//
// pybind11's own caster for an Eigen::SparseMatrix copies its three compressed arrays
// (values, inner indices, outer starts) into new numpy arrays, whatever the return value
// policy says. The helpers below build the scipy matrix on the C++ buffers themselves:
// `base` is set as the numpy base of every array, and is what keeps the buffers alive.

// The scipy.sparse csc_matrix (csr_matrix for a row major `mat`) whose arrays are the ones
// of `mat`, which must be compressed. Read only unless `writeable`.
template<typename SpMat>
py::object sparse_view_to_scipy(const SpMat & mat, py::handle base, bool writeable)
{
    using Scalar = typename SpMat::Scalar;
    using StorageIndex = typename SpMat::StorageIndex;
    const py::ssize_t nnz = static_cast<py::ssize_t>(mat.nonZeros());
    // an empty matrix may have no value buffer at all: numpy then allocates a (0,) array
    // of its own, which is just as good
    py::array_t<Scalar> data(nnz, mat.valuePtr(), base);
    py::array_t<StorageIndex> indices(nnz, mat.innerIndexPtr(), base);
    py::array_t<StorageIndex> indptr(static_cast<py::ssize_t>(mat.outerSize()) + 1, mat.outerIndexPtr(), base);
    if(!writeable){
        for(py::array * arr : {static_cast<py::array *>(&data), static_cast<py::array *>(&indices), static_cast<py::array *>(&indptr)}){
            arr->attr("setflags")(py::arg("write") = false);
        }
    }
    py::object matrix_type = py::module_::import("scipy.sparse").attr(SpMat::IsRowMajor ? "csr_matrix" : "csc_matrix");
    return matrix_type(py::make_tuple(std::move(data), std::move(indices), std::move(indptr)),
                       py::arg("shape") = py::make_tuple(mat.rows(), mat.cols()),
                       py::arg("copy") = false);
}

// Hand a matrix the caller no longer needs (typically one a getter just built) to python
// without copying it: it is moved to the heap, and freed with the last array using it.
template<typename SpMat>
py::object sparse_to_scipy(SpMat mat)
{
    mat.makeCompressed();
    std::unique_ptr<SpMat> owned(new SpMat(std::move(mat)));
    py::capsule owner(owned.get(), [](void * ptr){ delete static_cast<SpMat *>(ptr); });
    const SpMat & ref = *owned.release();  // the capsule owns it now
    return sparse_view_to_scipy(ref, owner, true);
}

#endif // SPARSE_HELPERS_HPP
