# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import struct
import sys
import tempfile
import os
import unittest
import warnings

import numpy as np

import grid2op
from lightsim2grid.lightSimBackend import LightSimBackend

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _exotic_elements_fixture import build_exotic_elements_grid  # noqa: E402

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

    def _aux_test_binary(self, fun_name, fun_comp, grid=None):
        """`grid` defaults to a full l2rpn_idf_2023 grid2op env's LSGrid; pass one
        explicitly (eg `build_exotic_elements_grid()`) to test an element type
        l2rpn_idf_2023 does not carry any of (SVC, HVDC -- see
        `_exotic_elements_fixture.py`)."""
        if grid is None:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self.env = grid2op.make("l2rpn_idf_2023", test=True, backend=LightSimBackend())
            grid = self.env.backend._grid
        els = getattr(grid, fun_name)()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, f"test_binary_{fun_name}.lsb")
            els.save_binary(path)
            els_reloaded = type(els).load_binary(path)

        class Struct:
            pass
        setattr(Struct, fun_name, lambda self: els_reloaded)
        diff_ = fun_comp(Struct(), grid)
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
        # l2rpn_idf_2023 carries no HVDC line at all -- use the exotic-elements
        # fixture instead, which has 3 (VSC droop, VSC no-droop, LCC).
        self._aux_test_binary("get_dclines", _compare_dclines, grid=build_exotic_elements_grid())

    def test_binary_svcs(self):
        # l2rpn_idf_2023 carries no SVC at all -- use the exotic-elements fixture.
        self._aux_test_binary("get_svcs", _compare_svcs, grid=build_exotic_elements_grid())

    def test_binary_algo_config(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_idf_2023", test=True, backend=LightSimBackend())
        grid = self.env.backend._grid

        ac_cfg = grid.get_ac_algo_config()
        ac_cfg.int_params = [1, 2, 3, 4]
        ac_cfg.real_params = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        grid.set_ac_algo_config(ac_cfg)
        # DC solver (DC_SparseLU) does not populate AlgoConfig at all (only the
        # NRAlgo family does): get_config()/set_config() are BaseAlgo's no-op.
        # Read back whatever is actually applied rather than assuming the
        # values above took effect, so this only tests that the round trip
        # preserves whatever the config actually is.
        dc_cfg = grid.get_dc_algo_config()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary.lsb")
            grid.save_binary(path)
            grid_1 = type(grid).load_binary(path)

        ac_cfg_1 = grid_1.get_ac_algo_config()
        assert list(ac_cfg_1.int_params) == list(ac_cfg.int_params)
        assert np.allclose(ac_cfg_1.real_params, ac_cfg.real_params)

        dc_cfg_1 = grid_1.get_dc_algo_config()
        assert list(dc_cfg_1.int_params) == list(dc_cfg.int_params)
        assert np.allclose(dc_cfg_1.real_params, dc_cfg.real_params)

    @staticmethod
    def _header_end(data):
        """Return the offset of the first state byte.

        File header layout (stable across binary format versions): 4-byte
        magic "LSB2", uint32 little-endian binary format version, then 4
        length-prefixed strings (uint32 little-endian length + raw bytes):
        type tag, major / medium / minor writer version.
        """
        off = 4 + 4  # magic + format version
        for _ in range(4):  # type tag + 3 version strings
            (str_len,) = struct.unpack_from("<I", data, off)
            off += 4 + str_len
        return off

    def _make_grid(self, env_name="l2rpn_idf_2023"):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make(env_name, test=True, backend=LightSimBackend())
        return self.env.backend._grid

    def test_cross_version_load(self):
        """Compatibility is decided by the binary *format* version, not the
        lightsim2grid version: a file whose writer-version strings differ
        (but whose format number matches) must load fine."""
        grid = self._make_grid()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary.lsb")
            grid.save_binary(path)
            with open(path, "rb") as f:
                data = bytearray(f.read())

            # doctor the 3 writer-version strings (they follow the type tag)
            off = 4 + 4
            (tag_len,) = struct.unpack_from("<I", bytes(data), off)
            off += 4 + tag_len
            for _ in range(3):
                (str_len,) = struct.unpack_from("<I", bytes(data), off)
                off += 4
                for i in range(off, off + str_len):
                    data[i] = ord("9")
                off += str_len

            other_version_path = os.path.join(tmpdir, "test_binary_other_version.lsb")
            with open(other_version_path, "wb") as f:
                f.write(bytes(data))

            grid_1 = type(grid).load_binary(other_version_path)
            assert len(compare_network_input(grid, grid_1)) == 0

    def test_cannot_load_wrong_format(self):
        """A file with an unsupported binary format number must be rejected
        with a clear error naming both versions and the format numbers."""
        grid = self._make_grid()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary.lsb")
            grid.save_binary(path)
            with open(path, "rb") as f:
                data = bytearray(f.read())

            # doctor the uint32 binary format version right after the magic
            struct.pack_into("<I", data, 4, 999999)

            bad_path = os.path.join(tmpdir, "test_binary_bad_format.lsb")
            with open(bad_path, "wb") as f:
                f.write(bytes(data))

            with self.assertRaises(RuntimeError) as cm:
                type(grid).load_binary(bad_path)
            assert "binary format" in str(cm.exception)
            assert "999999" in str(cm.exception)

    def test_corrupted_count(self):
        """A corrupted element-count field must raise a clean RuntimeError
        (checked against the real file size before any allocation), not a
        MemoryError / out-of-memory crash."""
        grid = self._make_grid()
        loads = grid.get_loads()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary_loads.lsb")
            loads.save_binary(path)
            with open(path, "rb") as f:
                data = bytearray(f.read())

            # the first state field of a LoadContainer is a string-vector
            # (names) whose uint64 count immediately follows the header:
            # replace it with an absurdly huge value
            struct.pack_into("<Q", data, self._header_end(bytes(data)), 2**62)

            bad_path = os.path.join(tmpdir, "test_binary_loads_bad_count.lsb")
            with open(bad_path, "wb") as f:
                f.write(bytes(data))

            with self.assertRaises(RuntimeError) as cm:
                type(loads).load_binary(bad_path)
            assert "truncated or corrupted" in str(cm.exception)

    def test_cannot_load_wrong_type(self):
        """LoadContainer and StorageContainer have identical StateRes layouts:
        without the type tag in the header, loading one as the other would
        silently succeed with an object of the wrong type."""
        grid = self._make_grid()
        loads = grid.get_loads()
        storages = grid.get_storages()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary_loads.lsb")
            loads.save_binary(path)

            with self.assertRaises(RuntimeError) as cm:
                type(storages).load_binary(path)
            assert "LoadContainer" in str(cm.exception)
            assert "StorageContainer" in str(cm.exception)

    def test_trailing_bytes_rejected(self):
        """Bytes remaining after the state has been fully read are treated as
        corruption (they would previously be silently ignored)."""
        grid = self._make_grid()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary.lsb")
            grid.save_binary(path)
            with open(path, "ab") as f:
                f.write(b"some trailing garbage")

            with self.assertRaises(RuntimeError) as cm:
                type(grid).load_binary(path)
            assert "remain after the end" in str(cm.exception)

    def test_atomic_save(self):
        """save_binary writes to a temporary file renamed into place on
        success: no temp file survives a successful save, and a failing save
        leaves a pre-existing good file untouched."""
        grid = self._make_grid()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary.lsb")
            grid.save_binary(path)
            assert not os.path.exists(path + ".lsb_tmp"), "temporary file left behind after a successful save"

            # make the next save fail (the temp file cannot be created
            # because a directory occupies its path), then check the
            # previously saved file is still intact and loadable
            os.mkdir(path + ".lsb_tmp")
            with self.assertRaises(RuntimeError):
                grid.save_binary(path)
            os.rmdir(path + ".lsb_tmp")

            grid_1 = type(grid).load_binary(path)
            assert len(compare_network_input(grid, grid_1)) == 0

    def test_non_atomic_save(self):
        """atomic=False writes the destination directly (the marginally
        faster path, no temporary file): the file must still round-trip."""
        grid = self._make_grid()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_binary_fast.lsb")
            grid.save_binary(path, atomic=False)
            assert not os.path.exists(path + ".lsb_tmp")
            grid_1 = type(grid).load_binary(path)
            assert len(compare_network_input(grid, grid_1)) == 0

            # overwriting an existing file the fast way works too
            grid.save_binary(path, atomic=False)
            grid_1 = type(grid).load_binary(path)
            assert len(compare_network_input(grid, grid_1)) == 0

    def test_corrupted_file(self):
        grid = self._make_grid()

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


class TestCorruptionSweep(unittest.TestCase):
    """Deterministic mini-fuzzer for load_binary: corrupt a valid file at
    EVERY byte offset and check the loader never does anything worse than
    raising a clean RuntimeError.

    Loading a corrupted file has exactly two acceptable outcomes: success
    (the corruption hit payload data and produced a different but structurally
    valid grid) or RuntimeError (the corruption was detected). A segfault, a
    MemoryError (a corrupted count being trusted before any bounds check) or
    any other exception type is a bug in the reader. This test is most potent
    in the ASan+UBSan CI job (.github/workflows/sanitizers.yml), where any
    out-of-bounds read also becomes a hard failure even if it would not have
    crashed a regular build.
    """

    @classmethod
    def setUpClass(cls):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        cls.grid_cls = type(env.backend._grid)
        cls.tmpdir = tempfile.TemporaryDirectory()
        path = os.path.join(cls.tmpdir.name, "sweep_reference.lsb")
        env.backend._grid.save_binary(path)
        with open(path, "rb") as f:
            cls.data = f.read()
        env.close()

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def _check_load(self, corrupted, description):
        bad_path = os.path.join(self.tmpdir.name, "sweep_corrupted.lsb")
        with open(bad_path, "wb") as f:
            f.write(corrupted)
        try:
            self.grid_cls.load_binary(bad_path)
        except RuntimeError:
            pass  # corruption detected and rejected cleanly: fine
        except BaseException as exc:
            self.fail(f"{description}: load_binary raised {type(exc).__name__} "
                      f"instead of RuntimeError: {exc}")

    def test_single_byte_flips(self):
        """Invert one byte at every offset of the file."""
        for offset in range(len(self.data)):
            corrupted = bytearray(self.data)
            corrupted[offset] ^= 0xFF
            self._check_load(bytes(corrupted), f"byte flip at offset {offset}")

    def test_uint64_stamps(self):
        """Stamp 8 bytes of 0xFF at every offset: turns any count/length
        field it lands on into a huge value, the worst case for the
        bounds-checked allocations."""
        for offset in range(len(self.data) - 7):
            corrupted = bytearray(self.data)
            corrupted[offset:offset + 8] = b"\xff" * 8
            self._check_load(bytes(corrupted), f"0xFF uint64 stamp at offset {offset}")

    def test_truncations(self):
        """Cut the file at every length: must always raise RuntimeError
        (a strictly shorter file can never be a complete state)."""
        for length in range(len(self.data)):
            bad_path = os.path.join(self.tmpdir.name, "sweep_truncated.lsb")
            with open(bad_path, "wb") as f:
                f.write(self.data[:length])
            with self.assertRaises(RuntimeError,
                                   msg=f"no error for a file truncated at {length} bytes"):
                self.grid_cls.load_binary(bad_path)


class TestSerializedEnumValues(unittest.TestCase):
    """The integer values of these C++ enums are written verbatim into the
    binary files (and into pickles), so they are part of the serialized
    format even though no StateRes tuple changes when they are renumbered.

    If this test fails: either restore the previous numbering (if the change
    was accidental -- note that for AlgorithmType, inserting a member
    anywhere but at the very end shifts every following value), or, if the
    renumbering is deliberate, bump BINARY_FORMAT_VERSION in
    src/core/BinaryArchive.hpp (and regenerate the reference fixture of
    TestBinaryLayoutUnchanged), then update the expected values here.
    """

    def _check_enum(self, enum_cls, expected):
        actual = {name: int(value) for name, value in enum_cls.__members__.items()}
        assert actual == expected, (
            f"The integer values of {enum_cls.__name__} changed, and they are "
            f"serialized verbatim in binary files / pickles: bump "
            f"BINARY_FORMAT_VERSION in src/core/BinaryArchive.hpp if this is "
            f"deliberate (see this test's docstring). Expected {expected}, "
            f"got {actual}")

    def test_algorithm_type(self):
        from lightsim2grid.lightsim2grid_cpp import AlgorithmType
        self._check_enum(AlgorithmType, {
            "NR_SparseLU": 0, "NR_KLU": 1, "GaussSeidel": 2, "DC_SparseLU": 3,
            "GaussSeidelSynch": 4, "NR_NICSLU": 5,
            "NRSing_SparseLU": 6, "NRSing_KLU": 7, "NRSing_NICSLU": 8,
            "DC_KLU": 9, "DC_NICSLU": 10,
            "NR_CKTSO": 11, "NRSing_CKTSO": 12, "DC_CKTSO": 13,
            "FDPF_XB_SparseLU": 14, "FDPF_BX_SparseLU": 15,
            "FDPF_XB_KLU": 16, "FDPF_BX_KLU": 17,
            "FDPF_XB_NICSLU": 18, "FDPF_BX_NICSLU": 19,
            "FDPF_XB_CKTSO": 20, "FDPF_BX_CKTSO": 21,
            "Custom": 22,
        })

    def test_svc_regulation_mode(self):
        from lightsim2grid.lightsim2grid_cpp import SvcContainer
        self._check_enum(SvcContainer.RegulationMode,
                         {"OFF": 0, "VOLTAGE": 1, "REACTIVE_POWER": 2})

    def test_hvdc_converters_mode(self):
        from lightsim2grid.lightsim2grid_cpp import HvdcLineContainer
        self._check_enum(HvdcLineContainer.ConvertersMode,
                         {"SIDE_1_RECTIFIER": 0, "SIDE_2_RECTIFIER": 1})

    def test_converter_station_type(self):
        from lightsim2grid.lightsim2grid_cpp import ConverterStationInfo
        self._check_enum(ConverterStationInfo.ConverterType,
                         {"VSC": 0, "LCC": 1})


# reference file saved with binary format 3 (see BINARY_FORMAT_VERSION in
# src/core/BinaryArchive.hpp) + a few values it is known to contain, used by
# TestBinaryLayoutUnchanged below. Regenerate (only after a deliberate format
# bump) with: python -m lightsim2grid.tests.test_binary_serialization regen
FIXTURE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "binary_format_fixture",
                            "case14_sandbox_format3.lsb")
FIXTURE_ENV_NAME = "l2rpn_case14_sandbox"
FIXTURE_N_SUB = 14
FIXTURE_N_LOAD = 11
FIXTURE_N_GEN = 6
FIXTURE_N_LINE_PLUS_TRAFO = 20


class TestBinaryLayoutUnchanged(unittest.TestCase):
    """Guards the stability of the binary layout: loads a reference file
    committed to the repository and checks its content.

    If this test starts failing after a code change, that change modified the
    serialized layout (some StateRes tuple, or the encoding in
    BinaryArchive.hpp/.cpp) without bumping BINARY_FORMAT_VERSION: files
    saved by released lightsim2grid versions would be silently mis-read.
    Bump BINARY_FORMAT_VERSION in src/core/BinaryArchive.hpp, then
    regenerate this fixture (see FIXTURE_PATH above).
    """

    def test_reference_file_still_loads(self):
        from lightsim2grid.network import LSGrid
        try:
            grid = LSGrid.load_binary(FIXTURE_PATH)
        except RuntimeError as exc:
            if "binary format" in str(exc):
                self.fail(
                    "The committed reference file uses an older BINARY_FORMAT_VERSION. "
                    "If the format was bumped on purpose, regenerate the fixture "
                    f"({FIXTURE_PATH}) and update the expected values next to it. "
                    f"Original error: {exc}")
            self.fail(
                "The committed reference file no longer loads: the binary layout "
                "changed without bumping BINARY_FORMAT_VERSION (src/core/BinaryArchive.hpp). "
                f"Original error: {exc}")

        # a mis-aligned read can also produce garbage instead of raising:
        # check a few known values of the reference grid
        assert len(grid.get_loads()) == FIXTURE_N_LOAD
        assert len(grid.get_generators()) == FIXTURE_N_GEN
        assert len(grid.get_lines()) + len(grid.get_trafos()) == FIXTURE_N_LINE_PLUS_TRAFO
        for load in grid.get_loads():
            assert load.connected
        # the reference grid must still be able to run a powerflow
        nb_bus_total = FIXTURE_N_SUB * 2
        V = np.ones(nb_bus_total, dtype=complex)
        V = grid.ac_pf(V, 10, 1e-8)
        assert V.shape[0] > 0, "ac powerflow diverged on the reference grid"

    def test_roundtrip_matches_reference_env(self):
        """The fixture must stay equivalent to a freshly-built grid of the
        same environment (guards against regenerating it with wrong data)."""
        from lightsim2grid.network import LSGrid
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make(FIXTURE_ENV_NAME, test=True, backend=LightSimBackend())
        grid_ref = env.backend._grid
        grid_fixture = LSGrid.load_binary(FIXTURE_PATH)
        assert len(compare_network_input(grid_ref, grid_fixture)) == 0


def _regenerate_fixture():
    """Manual helper: regenerate the committed reference file after a
    deliberate BINARY_FORMAT_VERSION bump (never for any other reason)."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env = grid2op.make(FIXTURE_ENV_NAME, test=True, backend=LightSimBackend())
    os.makedirs(os.path.dirname(FIXTURE_PATH), exist_ok=True)
    env.backend._grid.save_binary(FIXTURE_PATH)
    print(f"fixture regenerated: {FIXTURE_PATH}")
    print("update the FIXTURE_* expected values and the file name if the format number changed")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "regen":
        _regenerate_fixture()
    else:
        unittest.main()
