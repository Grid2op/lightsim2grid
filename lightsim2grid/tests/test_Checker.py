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

    def aux_test_check(self, qlim=False):
        mismatch = self.env.backend._grid.check_solution(self.env.backend.V, qlim)
        assert mismatch.shape == (28, )
        assert np.all(np.abs(mismatch) <= self.env.backend.tol)

        mismatch = self.env.backend._grid.check_solution(2*self.env.backend.V, qlim)
        assert mismatch.shape == (28, )
        assert np.any(np.abs(mismatch) >= self.env.backend.tol)

    def test_check_withoutqlims(self):
        self.aux_test_check(False)

    def test_check_withqlims(self):
        with self.assertRaises(AssertionError):
            # this should raise because Q lim are not met
            self.aux_test_check(True)

        # and now i check that there are mismatch only due to qlim
        tol = self.env.backend.tol
        tol_gen_q = 1e-5  #because of the algo to compute the gen_q i am using in ls
        qlim = True
        mismatch = self.env.backend._grid.check_solution(self.env.backend.V, qlim)
        _, q_gen, _ = self.env.backend.generators_info()
        gen_qmin = self.env.backend.init_pp_backend._grid.gen["min_q_mvar"].values
        gen_qmax = self.env.backend.init_pp_backend._grid.gen["max_q_mvar"].values

        # check the the stuff that "breaks" the limits are indeed due to issue with qmin / qmax
        mismatch_q = np.imag(mismatch)
        gen_bus = self.env.backend.init_pp_backend._grid.gen["bus"].values
        # I lack q when i should have generated too much
        not_enough = q_gen > gen_qmax
        assert np.all(mismatch_q[gen_bus][not_enough] < 0.)
        assert np.all(np.abs(mismatch_q[gen_bus][not_enough] + (q_gen[not_enough] - gen_qmax[not_enough])) <= tol_gen_q)
        # I have too much q when i could not absorb enough
        too_much = q_gen < gen_qmin
        assert np.all(mismatch_q[gen_bus][too_much] > 0.)
        assert np.all(np.abs(mismatch_q[gen_bus][too_much] + (q_gen[too_much] - gen_qmin[too_much])) <= tol_gen_q)
        # I have no mismatch when generator
        gen_ok = (q_gen >= gen_qmin) & (q_gen <= gen_qmax)
        assert np.all(mismatch_q[gen_bus][gen_ok] <= tol)
        # no problem when there are no generators
        assert np.all(mismatch_q[~np.isin(np.arange(mismatch.shape[0]), gen_bus)] <= tol)

if __name__ == "__main__":
    unittest.main()
