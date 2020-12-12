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


if __name__ == "__main__":
    unittest.main()
