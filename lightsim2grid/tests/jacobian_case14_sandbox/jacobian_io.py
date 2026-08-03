# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
"""Pickle-free (de)serialization for the jacobian test fixtures.

The fixtures are a nested ``dict`` mixing python scalars, numpy arrays and
scipy CSC sparse matrices. They used to be stored with ``pickle``, but a pickle
produced with numpy >= 2 references ``numpy._core`` and cannot be read back with
numpy < 2 (and vice-versa).

Here everything is stored in a single ``.npz`` archive instead: the ``.npy``
binary format is stable across numpy versions and contains no module
references, so a file saved with any numpy can be read with any other one.
A small JSON manifest describes the structure so it can be rebuilt on load.
"""
import json

import numpy as np
import scipy.sparse as sp


def save_res(path, res):
    """Serialize the nested ``res`` mapping to a single ``.npz`` archive."""
    flat = {}
    counter = [0]

    def _new_key():
        key = f"arr_{counter[0]}"
        counter[0] += 1
        return key

    def _store(arr):
        key = _new_key()
        flat[key] = np.asarray(arr)
        return key

    def _walk(obj):
        if isinstance(obj, dict):
            return {"kind": "dict", "items": {k: _walk(v) for k, v in obj.items()}}
        if sp.issparse(obj):
            mat = obj.tocsc()
            return {
                "kind": "csc",
                "data": _store(mat.data),
                "indices": _store(mat.indices),
                "indptr": _store(mat.indptr),
                "shape": _store(np.asarray(mat.shape, dtype=np.int64)),
            }
        if isinstance(obj, np.ndarray):
            return {"kind": "array", "key": _store(obj)}
        # python / numpy scalar (e.g. the convergence tolerance)
        return {"kind": "scalar", "key": _store(obj)}

    manifest = _walk(res)
    flat["manifest"] = np.frombuffer(json.dumps(manifest).encode("utf-8"), dtype=np.uint8)
    np.savez_compressed(path, **flat)


def load_res(path):
    """Inverse of :func:`save_res`. Never uses pickle (``allow_pickle=False``)."""
    with np.load(path, allow_pickle=False) as npz:
        manifest = json.loads(bytes(npz["manifest"]).decode("utf-8"))

        def _build(node):
            kind = node["kind"]
            if kind == "dict":
                return {k: _build(v) for k, v in node["items"].items()}
            if kind == "csc":
                shape = tuple(int(s) for s in npz[node["shape"]])
                return sp.csc_matrix(
                    (npz[node["data"]], npz[node["indices"]], npz[node["indptr"]]),
                    shape=shape,
                )
            if kind == "array":
                return npz[node["key"]]
            return npz[node["key"]].item()  # scalar

        return _build(manifest)
