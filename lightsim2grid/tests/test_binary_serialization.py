# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import struct
import tempfile
import os
import unittest
import warnings

import numpy as np

import grid2op
from lightsim2grid.lightSimBackend import LightSimBackend
from lightsim2grid.network.compare_lsgrid import (
    compare_network_input,
    _compare_loads,
    _compare_lines,
    _compare_generators,
    _compare_shunts,
    _compare_storages,
    _compare_substations,
    _compare_trafos,
    _compare_static_generators,
    _compare_dclines,
    _compare_svcs)


class TestBinarySerialization(unittest.TestCase):
    """Round-trip tests for LSGrid.save_binary / LSGrid.load_binary, the fast
    additive alternative to pickle. Mirrors the structure of test_pickleable.py
    (same fixture, same comparison helpers), so pickle and binary stay in sync
    behaviourally even though they are separate, independent code paths.
    """

    def _aux_test_2sides(self, grid1, grid2, method_name, test_results=False):
        assert len(getattr(grid1, method_name)()) == len(getattr(grid2, method_name)())
        for line1, line2 in zip(getattr(grid1, method_name)(), getattr(grid2, method_name)()):
            for attr_nm in ["id", "name", "sub1_id", "sub2_id",
                            "pos1_topo_vect", "pos2_topo_vect",
                            "connected_global", "connected1", "connected2",
                            "bus1_id", "bus2_id", "has_res",
                            ]:
                assert getattr(line1, attr_nm) == getattr(line2, attr_nm), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"

            for attr_nm in ["r_pu", "x_pu", "h1_pu", "h2_pu",
                            "yac_11", "yac_12", "yac_21", "yac_22",
                            "ydc_11", "ydc_12", "ydc_21", "ydc_22",
                            ]:
                assert np.allclose(getattr(line1, attr_nm), getattr(line2, attr_nm)), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"

            if test_results:
                for attr_nm in ["res_p1_mw", "res_q1_mvar", "res_theta1_deg", "res_v1_kv", "res_a1_ka",
                                "res_p2_mw", "res_q2_mvar", "res_theta2_deg", "res_v2_kv", "res_a2_ka"]:
                    assert np.allclose(getattr(line1, attr_nm), getattr(line2, attr_nm)), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"

    def _aux_test_1side(self, grid1, grid2, method_name, test_results=False,
                        add_attr_int=None,
                        add_attr_float=None):
        li_attr_to_test_int = [
            "id", "name", "sub_id",
            "pos_topo_vect",
            "connected",
            "bus_id",
        ]
        if add_attr_int is not None:
            li_attr_to_test_int += add_attr_int

        li_attr_to_test_float = []
        if add_attr_float is not None:
            li_attr_to_test_float += add_attr_float
        assert len(getattr(grid1, method_name)()) == len(getattr(grid2, method_name)())
        for line1, line2 in zip(getattr(grid1, method_name)(), getattr(grid2, method_name)()):
            for attr_nm in li_attr_to_test_int:
                assert getattr(line1, attr_nm) == getattr(line2, attr_nm), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"

            for attr_nm in li_attr_to_test_float:
                assert np.allclose(getattr(line1, attr_nm), getattr(line2, attr_nm)), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"

            if not test_results:
                continue
            li_attr_to_test = ["res_p_mw", "res_q_mvar", "res_theta_deg", "res_v_kv"]
            for attr_nm in li_attr_to_test:
                if not np.allclose(getattr(line1, attr_nm), getattr(line2, attr_nm)):
                    diff_ = np.abs(getattr(line1, attr_nm) - getattr(line2, attr_nm))
                    raise AssertionError(f"{method_name} error for {attr_nm} for gen id {line1.id} ({line1.name}): {getattr(line1, attr_nm)} "
                                            f"vs {getattr(line2, attr_nm)} -> {diff_}")

    def aux_test_2sides(self, grid1, grid2, test_results=False):
        self._aux_test_2sides(grid1, grid2, "get_lines", test_results)
        self._aux_test_2sides(grid1, grid2, "get_trafos", test_results)

    def aux_test_1side(self, grid1, grid2, test_results=False):
        self._aux_test_1side(grid1, grid2, "get_loads", test_results)
        self._aux_test_1side(grid1, grid2, "get_generators", test_results,
                             add_attr_int=["is_slack", "voltage_regulator_on"],
                             add_attr_float=["slack_weight", "target_vm_pu", "min_q_mvar", "max_q_mvar"])
        self._aux_test_1side(grid1, grid2, "get_storages", test_results)
        self._aux_test_1side(grid1, grid2, "get_shunts", test_results)

    def test_save_load(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_idf_2023", test=True, backend=LightSimBackend())

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary.lsb")
            grid = self.env.backend._grid
            grid.save_binary(path)
            grid_1 = type(grid).load_binary(path)

            assert grid_1.get_algo_type() == grid.get_algo_type()
            assert grid_1.get_dc_solver_type() == grid.get_dc_solver_type()

            self.aux_test_2sides(grid, grid_1)
            self.aux_test_1side(grid, grid_1)
            tmp = compare_network_input(grid, grid_1)
            assert len(tmp) == 0

            nb_bus_total = self.env.n_sub * 2
            max_it = 10
            tol = 1e-8

            # test dc_pf
            V_0 = np.ones(nb_bus_total, dtype=complex)
            V_1 = V_0.copy()
            V_0 = grid.dc_pf(V_0, max_it, tol)
            V_1 = grid_1.dc_pf(V_1, max_it, tol)

            assert np.all(np.abs(V_0 - V_1) <= 1e-7), "dc pf does not lead to same results"
            self.aux_test_2sides(grid, grid_1, True)
            self.aux_test_1side(grid, grid_1, True)

            # test ac_pf
            V_0 = grid.ac_pf(V_0, max_it, tol)
            V_1 = grid_1.ac_pf(V_1, max_it, tol)
            assert np.all(np.abs(V_0 - V_1) <= 1e-7), "ac pf does not lead to same results"
            self.aux_test_2sides(grid, grid_1, True)
            self.aux_test_1side(grid, grid_1, True)

    def _aux_test_binary(self, fun_name, fun_comp):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_idf_2023", test=True, backend=LightSimBackend())
        els = getattr(self.env.backend._grid, fun_name)()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, f"test_binary_{fun_name}.lsb")
            els.save_binary(path)
            els_reloaded = type(els).load_binary(path)

        class Struct:
            pass
        setattr(Struct, fun_name, lambda self: els_reloaded)
        diff_ = fun_comp(Struct(), self.env.backend._grid)
        assert len(diff_) == 0

    def test_binary_loads(self):
        self._aux_test_binary("get_loads", _compare_loads)

    def test_binary_lines(self):
        self._aux_test_binary("get_lines", _compare_lines)

    def test_binary_trafos(self):
        self._aux_test_binary("get_trafos", _compare_trafos)

    def test_binary_storages(self):
        self._aux_test_binary("get_storages", _compare_storages)

    def test_binary_generators(self):
        self._aux_test_binary("get_generators", _compare_generators)

    def test_binary_shunts(self):
        self._aux_test_binary("get_shunts", _compare_shunts)

    def test_binary_substations(self):
        self._aux_test_binary("get_substations", _compare_substations)

    def test_binary_sgens(self):
        self._aux_test_binary("get_static_generators", _compare_static_generators)

    def test_binary_hvdc(self):
        self._aux_test_binary("get_dclines", _compare_dclines)

    def test_binary_svcs(self):
        self._aux_test_binary("get_svcs", _compare_svcs)

    def test_cannot_load_wrong_version(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_idf_2023", test=True, backend=LightSimBackend())
        grid = self.env.backend._grid

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary.lsb")
            grid.save_binary(path)
            with open(path, "rb") as f:
                data = bytearray(f.read())

            # file layout: 4-byte magic "LSB1", then 3 length-prefixed strings
            # (uint32 little-endian length + raw bytes): major, medium, minor
            # version. Doctor the major-version bytes so they cannot match the
            # currently installed lightsim2grid version.
            maj_len = struct.unpack("<I", bytes(data[4:8]))[0]
            for i in range(8, 8 + maj_len):
                data[i] = ord("9")

            bad_path = os.path.join(tmpdir, "test_binary_bad_version.lsb")
            with open(bad_path, "wb") as f:
                f.write(bytes(data))

            with self.assertRaises(RuntimeError):
                type(grid).load_binary(bad_path)

    def test_corrupted_file(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_idf_2023", test=True, backend=LightSimBackend())
        grid = self.env.backend._grid

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary.lsb")
            grid.save_binary(path)
            with open(path, "rb") as f:
                data = f.read()

            # truncated file: valid header, but the tuple body cuts off mid-read
            trunc_path = os.path.join(tmpdir, "test_binary_truncated.lsb")
            with open(trunc_path, "wb") as f:
                f.write(data[:20])
            with self.assertRaises(RuntimeError):
                type(grid).load_binary(trunc_path)

            # completely unrelated / garbage file (bad magic number)
            garbage_path = os.path.join(tmpdir, "test_binary_garbage.lsb")
            with open(garbage_path, "wb") as f:
                f.write(b"this is not a lightsim2grid binary file at all, just some garbage bytes")
            with self.assertRaises(RuntimeError):
                type(grid).load_binary(garbage_path)


if __name__ == "__main__":
    unittest.main()
