# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The HVDC machinery must be a strict no-op when no HVDC line is present.

The ``Hvdc`` NRSystem extension is always compiled into the solver aliases; with zero
declared lines it must loop over nothing and leave the solution, the Jacobian sparsity
and the iteration count bit-identical to a model that never touched the HVDC API.
"""

import os
import sys
import unittest
import warnings
import numpy as np
import pandapower.networks as pn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower

_EMPTY_I = np.array([], dtype=np.int32)
_EMPTY_F = np.array([])


def _declare_no_hvdc(model):
    model.init_hvdc_lines(
        _EMPTY_I, _EMPTY_I, [], [], _EMPTY_F, _EMPTY_F, [], [], _EMPTY_F, _EMPTY_F,
        _EMPTY_F, _EMPTY_F, _EMPTY_F, _EMPTY_F, _EMPTY_F, _EMPTY_F, _EMPTY_F, _EMPTY_F,
        [], _EMPTY_F, _EMPTY_F, _EMPTY_F, [], _EMPTY_F, _EMPTY_F, _EMPTY_F, _EMPTY_F)
    model.tell_solver_need_reset()


class TestNoHvdcBitIdentical(unittest.TestCase):
    def _cases(self):
        return [("case14", pn.case14), ("case118", pn.case118)]

    def _models(self, factory):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            net = factory()
            plain = init_from_pandapower(net)
            with_empty = init_from_pandapower(net)
        _declare_no_hvdc(with_empty)
        return net, plain, with_empty

    def test_ac_bit_identical(self):
        for name, factory in self._cases():
            with self.subTest(case=name):
                net, plain, with_empty = self._models(factory)
                nb = net.bus.shape[0]
                Vp = plain.ac_pf(np.ones(nb, dtype=np.complex128), 30, 1e-10)
                Ve = with_empty.ac_pf(np.ones(nb, dtype=np.complex128), 30, 1e-10)
                self.assertGreater(Vp.shape[0], 0)
                self.assertEqual(np.abs(Vp - Ve).max(), 0.0)
                self.assertEqual(plain.get_solver().get_nb_iter(),
                                 with_empty.get_solver().get_nb_iter())
                self.assertEqual(plain.get_solver().get_J().nnz,
                                 with_empty.get_solver().get_J().nnz)

    def test_dc_bit_identical(self):
        for name, factory in self._cases():
            with self.subTest(case=name):
                net, plain, with_empty = self._models(factory)
                nb = net.bus.shape[0]
                Vp = plain.dc_pf(np.ones(nb, dtype=np.complex128), 1, 1.0)
                Ve = with_empty.dc_pf(np.ones(nb, dtype=np.complex128), 1, 1.0)
                self.assertGreater(Vp.shape[0], 0)
                self.assertEqual(np.abs(Vp - Ve).max(), 0.0)


if __name__ == "__main__":
    unittest.main()
