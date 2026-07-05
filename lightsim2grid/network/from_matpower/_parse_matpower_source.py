# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import numpy as np


def _as_array(obj):
    # accepts a numpy array, a pandas DataFrame (via ".values") or anything array-like
    if hasattr(obj, "values"):
        arr = np.asarray(obj.values)
    else:
        arr = np.asarray(obj)
    if arr.ndim == 1:
        # a single-row matpower matrix (e.g. exactly one branch) can come back as a
        # flat 1d array, in particular through `scipy.io.loadmat(..., squeeze_me=True)`
        # which squeezes away the (1, ncols) leading dimension. Restore it so that
        # column indexing (`arr[:, SOME_COL]`) always works.
        arr = arr.reshape(1, -1)
    return arr


def _load_from_m(path):
    """Parse a MATPOWER ``.m`` case file using the optional ``matpowercaseframes`` package."""
    try:
        from matpowercaseframes import CaseFrames
    except ImportError as exc_:
        raise ImportError("Reading a matpower \".m\" file requires the optional "
                           "\"matpowercaseframes\" package. Install it with "
                           "`pip install matpowercaseframes` (or `pip install lightsim2grid[matpower]`).") from exc_
    cf = CaseFrames(path)
    return _as_array(cf.bus), _as_array(cf.gen), _as_array(cf.branch), float(cf.baseMVA)


def _load_from_mat(path):
    """Parse a MATPOWER ``.mat`` case file using the optional ``scipy`` package."""
    try:
        import scipy.io
    except ImportError as exc_:
        raise ImportError("Reading a matpower \".mat\" file requires the optional "
                           "\"scipy\" package. Install it with "
                           "`pip install scipy` (or `pip install lightsim2grid[matpower]`).") from exc_
    mat = scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)
    if "mpc" in mat:
        mpc = mat["mpc"]
    else:
        candidates = [k for k in mat.keys() if not k.startswith("__")]
        if len(candidates) != 1:
            raise RuntimeError(f"Could not find a unique matpower case struct in \"{path}\". "
                                f"Found top level variables: {candidates}. Expected a single "
                                f"variable named \"mpc\" (or a single top level struct).")
        mpc = mat[candidates[0]]
    if hasattr(mpc, "bus"):
        # scipy loaded it as a matlab struct (mat_struct instance)
        bus, gen, branch, baseMVA = mpc.bus, mpc.gen, mpc.branch, mpc.baseMVA
    else:
        # scipy loaded it as a structured / record array
        bus, gen, branch, baseMVA = mpc["bus"], mpc["gen"], mpc["branch"], mpc["baseMVA"]
    return _as_array(bus), _as_array(gen), _as_array(branch), float(np.asarray(baseMVA).item())


def _load_from_dict_like(source):
    """Accept a plain dict (e.g. as returned by pypower's ``caseN()`` functions) or any
    object exposing ``.bus`` / ``.gen`` / ``.branch`` / ``.baseMVA`` (e.g. a
    ``matpowercaseframes.CaseFrames`` instance built by the user beforehand).
    """
    if isinstance(source, dict):
        bus, gen, branch, baseMVA = source["bus"], source["gen"], source["branch"], source["baseMVA"]
    else:
        bus, gen, branch, baseMVA = source.bus, source.gen, source.branch, source.baseMVA
    return _as_array(bus), _as_array(gen), _as_array(branch), float(baseMVA)


def load_matpower_data(source):
    """
    Normalize any accepted matpower "source" into plain numpy arrays.

    Parameters
    ----------
    source:
        Either a path (str or os.PathLike) to a ".m" or ".mat" matpower case file,
        or an already parsed matpower case (a dict with "bus" / "gen" / "branch" /
        "baseMVA" keys, e.g. from `pypower`, or any object exposing these as attributes,
        e.g. a `matpowercaseframes.CaseFrames` instance).

    Returns
    -------
    bus, gen, branch: numpy arrays
        The raw matpower matrices (rows = elements, columns = matpower's
        positional column layout, see `lightsim2grid.network.from_matpower._my_const`)
    baseMVA: float
        The system base MVA
    """
    if isinstance(source, (str, os.PathLike)):
        path = os.fspath(source)
        if path.endswith(".m"):
            return _load_from_m(path)
        elif path.endswith(".mat"):
            return _load_from_mat(path)
        else:
            raise RuntimeError(f"Unsupported matpower file extension for \"{path}\", "
                               f"expected a path ending in \".m\" or \".mat\".")
    return _load_from_dict_like(source)
