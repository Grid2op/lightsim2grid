# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import unittest
import grid2op
import warnings
import pdb
import numpy as np
from lightsim2grid import LightSimBackend

SparseLUSolver_AVAILBLE = False
try:
    from lightsim2grid_cpp import SparseLUSolver
    SparseLUSolver_AVAILBLE = True
except ImportError:
    # KLU solver is not available, these tests cannot be carried out
    pass


class MakeTests(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make(backend=LightSimBackend(), test=True)

    def tearDown(self) -> None:
        self.env.close()

    def aux_gen_ok(self, el, gen_id, tol, gen_p, gen_q, gen_v):
        assert el.connected is True, f"gen {gen_id} is not connected"
        assert el.bus_id == self.env.backend.gen_to_subid[gen_id], f"gen {gen_id} is connected to wrong bus"
        # assert np.abs(el.target_p_mw - gen_p[gen_id]) <= tol  # do not work on the slack bus
        assert np.abs(el.target_vm_pu - gen_v[gen_id] / self.env.backend.prod_pu_to_kv[gen_id]) <= tol, \
            f"gen {gen_id} has wrong voltage setpoint"
        q_min_ref = self.env.backend.init_pp_backend._grid.gen["min_q_mvar"].values[gen_id]
        assert np.abs(el.min_q_mvar - q_min_ref) <= tol, f"gen {gen_id} has wrong qmin"
        q_max_ref = self.env.backend.init_pp_backend._grid.gen["max_q_mvar"].values[gen_id]
        assert np.abs(el.max_q_mvar - q_max_ref) <= tol, f"gen {gen_id} has wrong qmax"
        assert el.has_res is True, f"gen {gen_id} has no results available"
        assert np.abs(el.res_p_mw - gen_p[gen_id]) <= tol, f"gen {gen_id} has wrong res_p"
        assert np.abs(el.res_q_mvar - gen_q[gen_id]) <= tol, f"gen {gen_id} has wrong res_q"
        assert np.abs(el.res_v_kv - gen_v[gen_id]) <= tol, f"gen {gen_id} has wrong res_v"

    def test_getters_gen(self):
        """test that the generators getter return the right values"""
        data_gen = self.env.backend._grid.get_generators()
        gen_p, gen_q, gen_v = self.env.backend.generators_info()
        assert len(data_gen) == self.env.n_gen
        tol = 1e-5
        for el in data_gen:
            gen_id = el.id
            self.aux_gen_ok(el, gen_id, tol, gen_p, gen_q, gen_v)

            with self.assertRaises(AttributeError):
                # this should not be possible
                el.bus_id = 2

        gen_info = data_gen[0]
        self.aux_gen_ok(gen_info, 0, tol, gen_p, gen_q, gen_v)
        with self.assertRaises(ValueError):
            gen_info = data_gen[-1]

        with self.assertRaises(ValueError):
            gen_info = data_gen[self.env.n_gen]

    def aux_trafo_ok(self, el, trafo_id, tol, data_ref):
        """trafo_id: id for the trafo in grid2Op"""
        p_or, q_or, v_or, a_or = data_ref
        assert el.connected is True, f"gen {trafo_id} is not connected"
        assert el.bus_hv_id == self.env.backend.line_or_to_subid[trafo_id], \
            f"trafo {trafo_id} is connected to wrong bus (hv side)"
        assert el.bus_lv_id == self.env.backend.line_ex_to_subid[trafo_id], \
            f"trafo {trafo_id} is connected to wrong bus (lv side)"
        assert el.has_res, f"trafo {trafo_id} don't have any results"
        assert np.abs(el.res_p_hv_mw - p_or[trafo_id]) <= tol, f"trafo {trafo_id} has wrong p_hv"
        assert np.abs(el.res_q_hv_mvar - q_or[trafo_id]) <= tol, f"trafo {trafo_id} has wrong q_hv"
        assert np.abs(el.res_a_hv_ka - 0.001 * a_or[trafo_id]) <= tol, f"trafo {trafo_id} has wrong a_hv"

    def test_getters_trafo(self):
        """test that the trafo getter return the right values"""
        data_trafos = self.env.backend._grid.get_trafos()
        p_or, q_or, v_or, a_or = self.env.backend.lines_or_info()
        nb_trafo = 5
        nb_line = 15
        data_ref = (p_or, q_or, v_or, a_or)
        assert len(data_trafos) == nb_trafo
        tol = 1e-5
        for el in data_trafos:
            trafo_id = el.id
            self.aux_trafo_ok(el, trafo_id + nb_line, tol, data_ref)

            with self.assertRaises(AttributeError):
                # this should not be possible
                el.bus_id = 2

        gen_info = data_trafos[0]
        self.aux_trafo_ok(gen_info, 0 + nb_line, tol, data_ref)
        with self.assertRaises(ValueError):
            gen_info = data_trafos[-1]

        with self.assertRaises(ValueError):
            gen_info = data_trafos[nb_trafo]

    def aux_line_ok(self, el, line_id, tol, data_ref):
        """trafo_id: id for the trafo in grid2Op"""
        p_or, q_or, v_or, a_or = data_ref
        assert el.connected is True, f"line {line_id} is not connected"
        assert el.bus_or_id == self.env.backend.line_or_to_subid[line_id], \
            f"line {line_id} is connected to wrong bus (or side)"
        assert el.bus_ex_id == self.env.backend.line_ex_to_subid[line_id], \
            f"line {line_id} is connected to wrong bus (ex side)"
        assert el.has_res, f"line {line_id} don't have any results"
        assert np.abs(el.res_p_or_mw - p_or[line_id]) <= tol, f"line {line_id} has wrong p_or"
        assert np.abs(el.res_q_or_mvar - q_or[line_id]) <= tol, f"line {line_id} has wrong q_or"
        assert np.abs(el.res_a_or_ka - 0.001 * a_or[line_id]) <= tol, f"line {line_id} has wrong a_or"

    def test_getters_line(self):
        """test that the line getter return the right values"""
        data_lines = self.env.backend._grid.get_lines()
        p_or, q_or, v_or, a_or = self.env.backend.lines_or_info()
        nb_trafo = 5
        nb_line = 15
        data_ref = (p_or, q_or, v_or, a_or)
        assert len(data_lines) == nb_line
        tol = 1e-5
        for el in data_lines:
            line_id = el.id
            self.aux_line_ok(el, line_id, tol, data_ref)

            with self.assertRaises(AttributeError):
                # this should not be possible
                el.bus_id = 2

        gen_info = data_lines[0]
        self.aux_line_ok(gen_info, 0, tol, data_ref)
        with self.assertRaises(ValueError):
            line_info = data_lines[-1]

        with self.assertRaises(ValueError):
            line_info = data_lines[nb_line]


if __name__ == "__main__":
    unittest.main()
