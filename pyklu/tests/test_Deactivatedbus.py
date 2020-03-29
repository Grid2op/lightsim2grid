import os
import unittest
import numpy as np
import pdb
from scipy import sparse
from pyklu.initGridModel import init
import pandapower.networks as pn
import pandapower as pp


class MakeACTests(unittest.TestCase):
    def setUp(self):
        self.net = pn.case14()
        pp.create_bus(self.net, vn_kv=self.net.bus["vn_kv"][0])
        self.model = init(self.net)
        self.model.deactivate_bus(14)

    def test_ac_pf(self):
        V = self.model.ac_pf(np.ones(15, dtype=np.complex_), 10, 1e-8)
        assert V.shape[0] > 0

    def test_dc_pf(self):
        V = self.model.dc_pf(np.ones(15, dtype=np.complex_), 10, 1e-8)
        assert V.shape[0] > 0


if __name__ == "__main__":
    unittest.main()