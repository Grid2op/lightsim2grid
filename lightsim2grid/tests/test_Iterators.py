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
        assert np.abs(el.res_p - gen_p[gen_id]) <= tol, f"gen {gen_id} has wrong res_p"
        assert np.abs(el.res_q - gen_q[gen_id]) <= tol, f"gen {gen_id} has wrong res_q"
        assert np.abs(el.res_v - gen_v[gen_id]) <= tol, f"gen {gen_id} has wrong res_v"

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


if __name__ == "__main__":
    unittest.main()
