# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import numpy as np
import pdb
import pandapower.networks as pn
from lightsim2grid_cpp import PandaPowerConverter


class MakeTests(unittest.TestCase):
    def setUp(self):
        self.converter = PandaPowerConverter()
        self.tol = 1e-8

    def assert_equal(self, tmp, ref):
        assert np.max(np.abs(tmp - ref)) <= self.tol, f"maximum error {np.max(np.abs(tmp - ref))}"
        assert np.sum(np.abs(tmp - ref)) <= tmp.shape[0] * self.tol, f"average error {np.sum(np.abs(tmp - ref))}"

    def test_case6_data(self):
        net = pn.case6ww()
        self.converter.set_sn_mva(net.sn_mva)  # TODO raise an error if not set !
        self.converter.set_f_hz(net.f_hz)
        line_r, line_x, line_h = self.converter.get_line_param(
            net.line["r_ohm_per_km"].values * net.line["length_km"].values,
            net.line["x_ohm_per_km"].values * net.line["length_km"].values,
            net.line["c_nf_per_km"].values * net.line["length_km"].values,
            net.line["g_us_per_km"].values * net.line["length_km"].values,
            net.bus.loc[net.line["from_bus"]]["vn_kv"],
            net.bus.loc[net.line["to_bus"]]["vn_kv"]
            )
        res_r = np.array([0.001, 0.0005, 0.001, 0.0008, 0.0005, 0.0005, 0.001, 0.0007, 0.0012, 0.0002, 0.002])
        res_x = np.array([0.002, 0.002, 0.003, 0.003, 0.0025, 0.001, 0.003, 0.002, 0.0026, 0.001, 0.004])
        res_h = np.array([4.+0.j, 4.+0.j, 6.+0.j, 6.+0.j, 6.+0.j, 2.+0.j, 4.+0.j, 5.+0.j, 5.+0.j, 2.+0.j, 8.+0.j])
        self.assert_equal(line_r, res_r)
        self.assert_equal(line_x, res_x)
        self.assert_equal(line_h, res_h)

    def test_case30_data(self):
        net = pn.case30()
        self.converter.set_sn_mva(net.sn_mva)  # TODO raise an error if not set !
        self.converter.set_f_hz(net.f_hz)
        line_r, line_x, line_h = self.converter.get_line_param(
            net.line["r_ohm_per_km"].values * net.line["length_km"].values,
            net.line["x_ohm_per_km"].values * net.line["length_km"].values,
            net.line["c_nf_per_km"].values * net.line["length_km"].values,
            net.line["g_us_per_km"].values * net.line["length_km"].values,
            net.bus.loc[net.line["from_bus"]]["vn_kv"],
            net.bus.loc[net.line["to_bus"]]["vn_kv"]
            )
        res_r = np.array([0.0002, 0.0005, 0.    , 0.    , 0.    , 0.    , 0.    , 0.    ,
                           0.0012, 0.0007, 0.0009, 0.0022, 0.0006, 0.0008, 0.0011, 0.0006,
                           0.0003, 0.0009, 0.0003, 0.0003, 0.0007, 0.0001, 0.001 , 0.0001,
                           0.0012, 0.0013, 0.0019, 0.0025, 0.0011, 0.    , 0.0022, 0.0032,
                           0.0024, 0.0006, 0.0005, 0.0002, 0.0006, 0.0001, 0.0005, 0.0003,
                           0.0001])
        res_x = np.array([0.0006, 0.0019, 0.0021, 0.0056, 0.0021, 0.0011, 0.0026, 0.0014,
                           0.0026, 0.0013, 0.002 , 0.002 , 0.0017, 0.0019, 0.0022, 0.0013,
                           0.0007, 0.0021, 0.0008, 0.0007, 0.0015, 0.0002, 0.002 , 0.0004,
                           0.0018, 0.0027, 0.0033, 0.0038, 0.0021, 0.004 , 0.0042, 0.006 ,
                           0.0045, 0.002 , 0.002 , 0.0006, 0.0018, 0.0004, 0.0012, 0.0008,
                           0.0004])
        res_h = np.array([3.+0.j, 2.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,
                           0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 2.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,
                           0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,
                           0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,
                           0.+0.j, 2.+0.j, 2.+0.j, 1.+0.j, 2.+0.j, 0.+0.j, 1.+0.j, 1.+0.j,
                           0.+0.j])
        self.assert_equal(line_r, res_r)
        self.assert_equal(line_x, res_x)
        self.assert_equal(line_h, res_h)

    def test_case118_data(self):
        net = pn.case118()
        self.converter.set_sn_mva(net.sn_mva)  # TODO raise an error if not set !
        self.converter.set_f_hz(net.f_hz)
        line_r, line_x, line_h = self.converter.get_line_param(
            net.line["r_ohm_per_km"].values * net.line["length_km"].values,
            net.line["x_ohm_per_km"].values * net.line["length_km"].values,
            net.line["c_nf_per_km"].values * net.line["length_km"].values,
            net.line["g_us_per_km"].values * net.line["length_km"].values,
            net.bus.loc[net.line["from_bus"]]["vn_kv"],
            net.bus.loc[net.line["to_bus"]]["vn_kv"]
            )
        res_r = np.array([3.030e-04, 1.290e-04, 5.950e-05, 8.820e-05, 4.880e-04, 4.460e-04,
                           8.660e-05, 4.010e-04, 4.280e-04, 4.050e-04, 1.230e-04, 4.440e-04,
                           3.090e-04, 1.870e-04, 6.010e-04, 3.760e-05, 5.460e-05, 1.700e-04,
                           2.940e-04, 1.560e-04, 2.980e-04, 1.120e-04, 6.250e-04, 4.300e-04,
                           4.840e-04, 3.020e-04, 3.500e-04, 2.000e-04, 2.390e-04, 1.390e-04,
                           5.180e-04, 2.380e-04, 2.540e-04, 9.900e-05, 3.930e-04, 8.620e-05,
                           3.870e-04, 2.580e-04, 4.810e-04, 2.230e-04, 1.320e-04, 3.560e-04,
                           1.620e-04, 2.690e-04, 1.830e-04, 2.380e-04, 2.225e-04, 4.540e-04,
                           6.480e-04, 1.780e-04, 1.710e-04, 1.730e-04, 3.970e-04, 1.800e-04,
                           2.770e-04, 1.230e-04, 2.460e-04, 2.150e-04, 1.600e-04, 4.510e-04,
                           4.660e-04, 5.350e-04, 6.050e-04, 9.940e-05, 1.400e-04, 5.300e-04,
                           2.610e-04, 5.300e-04, 7.440e-04, 1.050e-04, 3.906e-04, 2.780e-04,
                           2.200e-04, 2.470e-04, 9.130e-05, 6.150e-04, 1.350e-04, 1.640e-04,
                           2.300e-05, 5.950e-04, 3.290e-04, 1.450e-04, 1.640e-04, 2.120e-04,
                           1.320e-04, 1.760e-05, 4.540e-04, 1.230e-04, 1.119e-04, 2.520e-04,
                           1.200e-04, 1.830e-04, 2.090e-04, 3.420e-04, 1.350e-04, 1.560e-04,
                           2.410e-04, 3.180e-04, 1.913e-04, 2.370e-04, 4.310e-05, 7.990e-05,
                           4.740e-04, 1.080e-04, 3.170e-04, 2.980e-04, 2.290e-04, 1.190e-04,
                           3.800e-04, 7.520e-04, 2.240e-05, 1.100e-04, 4.150e-04, 8.710e-05,
                           2.560e-05, 3.210e-04, 5.930e-04, 4.640e-05, 4.590e-05, 1.840e-04,
                           1.450e-04, 5.550e-04, 4.100e-04, 6.080e-04, 4.130e-04, 2.240e-04,
                           4.000e-04, 3.800e-04, 6.010e-04, 2.440e-05, 1.910e-04, 7.150e-04,
                           7.150e-04, 6.840e-04, 1.790e-04, 2.670e-04, 4.860e-04, 2.030e-04,
                           4.050e-04, 2.630e-04, 2.580e-05, 7.300e-04, 8.690e-04, 1.690e-04,
                           2.750e-05, 4.880e-05, 3.430e-04, 4.740e-04, 3.430e-04, 2.550e-04,
                           5.030e-04, 2.090e-04, 8.250e-04, 8.030e-04, 4.739e-04, 3.170e-04,
                           3.280e-04, 2.640e-05, 1.230e-04, 8.240e-05, 1.720e-05, 9.010e-05,
                           2.030e-04, 2.690e-05, 1.800e-04, 1.800e-04, 4.820e-04, 2.580e-04,
                           2.240e-04, 8.440e-04, 9.850e-04, 3.000e-04, 2.210e-05])
        res_x = np.array([9.990e-04, 4.240e-04, 1.960e-04, 3.550e-04, 1.960e-03, 1.800e-03,
                           4.540e-04, 1.323e-03, 1.410e-03, 1.220e-03, 4.060e-04, 1.480e-03,
                           1.010e-03, 6.160e-04, 1.999e-03, 1.240e-04, 2.440e-04, 4.850e-04,
                           1.050e-03, 7.040e-04, 8.530e-04, 3.665e-04, 1.320e-03, 1.480e-03,
                           1.600e-03, 6.410e-04, 1.230e-03, 1.020e-03, 1.730e-03, 7.120e-04,
                           1.880e-03, 9.970e-04, 8.360e-04, 5.050e-04, 1.581e-03, 3.400e-04,
                           1.272e-03, 8.480e-04, 1.580e-03, 7.320e-04, 4.340e-04, 1.820e-03,
                           5.300e-04, 8.690e-04, 9.340e-04, 1.080e-03, 7.310e-04, 2.060e-03,
                           2.950e-03, 5.800e-04, 5.470e-04, 8.850e-04, 1.790e-03, 8.130e-04,
                           1.262e-03, 5.590e-04, 1.120e-03, 7.070e-04, 5.250e-04, 2.040e-03,
                           1.584e-03, 1.625e-03, 2.290e-03, 3.780e-04, 5.470e-04, 1.830e-03,
                           7.030e-04, 1.830e-03, 2.444e-03, 2.880e-04, 1.813e-03, 7.620e-04,
                           7.550e-04, 6.400e-04, 3.010e-04, 2.030e-03, 6.120e-04, 7.410e-04,
                           1.040e-04, 1.950e-03, 1.400e-03, 4.810e-04, 5.440e-04, 8.340e-04,
                           4.370e-04, 7.980e-05, 1.801e-03, 5.050e-04, 4.930e-04, 1.170e-03,
                           3.940e-04, 8.490e-04, 9.700e-04, 1.590e-03, 4.920e-04, 8.000e-04,
                           1.080e-03, 1.630e-03, 8.550e-04, 9.430e-04, 5.040e-04, 8.600e-04,
                           1.563e-03, 3.310e-04, 1.153e-03, 9.850e-04, 7.550e-04, 5.400e-04,
                           1.244e-03, 2.470e-03, 1.020e-04, 4.970e-04, 1.420e-03, 2.680e-04,
                           9.400e-05, 1.060e-03, 1.680e-03, 5.400e-04, 2.080e-04, 6.050e-04,
                           4.870e-04, 1.830e-03, 1.350e-03, 2.454e-03, 1.681e-03, 9.010e-04,
                           1.356e-03, 1.270e-03, 1.890e-03, 3.050e-04, 6.250e-04, 3.230e-03,
                           3.230e-03, 1.860e-03, 5.050e-04, 7.520e-04, 1.370e-03, 5.880e-04,
                           1.635e-03, 1.220e-03, 3.220e-04, 2.890e-03, 2.910e-03, 7.070e-04,
                           9.550e-05, 1.510e-04, 9.660e-04, 1.340e-03, 9.660e-04, 7.190e-04,
                           2.293e-03, 6.880e-04, 2.510e-03, 2.390e-03, 2.158e-03, 1.450e-03,
                           1.500e-03, 1.350e-04, 5.610e-04, 3.760e-04, 2.000e-04, 9.860e-04,
                           6.820e-04, 3.020e-04, 9.190e-04, 9.190e-04, 2.180e-03, 1.170e-03,
                           1.015e-03, 2.778e-03, 3.240e-03, 1.270e-03, 4.115e-03])
        res_h = np.array([  2.54 +0.j,   1.082+0.j,   0.502+0.j,   0.878+0.j,   4.88 +0.j,
                             4.444+0.j,   1.178+0.j,   3.368+0.j,   3.6  +0.j,  12.4  +0.j,
                             1.034+0.j,   3.68 +0.j,  10.38 +0.j,   1.572+0.j,   4.978+0.j,
                             1.264+0.j,   0.648+0.j,   4.72 +0.j,   2.28 +0.j,   1.87 +0.j,
                             8.174+0.j,   3.796+0.j,   2.58 +0.j,   3.48 +0.j,   4.06 +0.j,
                             1.234+0.j,   2.76 +0.j,   2.76 +0.j,   4.7  +0.j,   1.934+0.j,
                             5.28 +0.j,  10.6  +0.j,   2.14 +0.j,   5.48 +0.j,   4.14 +0.j,
                             0.874+0.j,   3.268+0.j,   2.18 +0.j,   4.06 +0.j,   1.876+0.j,
                             1.11 +0.j,   4.94 +0.j,   5.44 +0.j,   2.3  +0.j,   2.54 +0.j,
                             2.86 +0.j,   1.876+0.j,   5.46 +0.j,   4.72 +0.j,   6.04 +0.j,
                             1.474+0.j,   2.4  +0.j,   4.76 +0.j,   2.16 +0.j,   3.28 +0.j,
                             1.464+0.j,   2.94 +0.j,   1.816+0.j,   5.36 +0.j,   5.41 +0.j,
                             4.07 +0.j,   4.08 +0.j,   6.2  +0.j,   0.986+0.j,   1.434+0.j,
                             4.72 +0.j,   1.844+0.j,   4.72 +0.j,   6.268+0.j,   0.76 +0.j,
                             4.61 +0.j,   2.02 +0.j,   2.   +0.j,   6.2  +0.j,   0.768+0.j,
                             5.18 +0.j,   1.628+0.j,   1.972+0.j,   0.276+0.j,   5.02 +0.j,
                             3.58 +0.j,   1.198+0.j,   1.356+0.j,   2.14 +0.j,   4.44 +0.j,
                             0.21 +0.j,   4.66 +0.j,   1.298+0.j,   1.142+0.j,   2.98 +0.j,
                             1.01 +0.j,   2.16 +0.j,   2.46 +0.j,   4.04 +0.j,   4.98 +0.j,
                             8.64 +0.j,   2.84 +0.j,  17.64 +0.j,   2.16 +0.j,   2.38 +0.j,
                            51.4  +0.j,  90.8  +0.j,   3.99 +0.j,   0.83 +0.j,  11.73 +0.j,
                             2.51 +0.j,   1.926+0.j,   1.426+0.j,   3.194+0.j,   6.32 +0.j,
                             0.268+0.j,   1.318+0.j,   3.66 +0.j,   0.568+0.j,   0.984+0.j,
                             2.7  +0.j,   4.2  +0.j,  42.2  +0.j,   0.55 +0.j,   1.552+0.j,
                             1.222+0.j,   4.66 +0.j,   3.44 +0.j,   6.068+0.j,   4.226+0.j,
                             2.24 +0.j,   3.32 +0.j,   3.16 +0.j,   4.72 +0.j, 116.2  +0.j,
                             1.604+0.j,   8.6  +0.j,   8.6  +0.j,   4.44 +0.j,   1.258+0.j,
                             1.874+0.j,   3.42 +0.j,   1.396+0.j,   4.058+0.j,   3.1  +0.j,
                           123.   +0.j,   7.38 +0.j,   7.3  +0.j,   2.02 +0.j,   0.732+0.j,
                             0.374+0.j,   2.42 +0.j,   3.32 +0.j,   2.42 +0.j,   1.788+0.j,
                             5.98 +0.j,   1.748+0.j,   5.69 +0.j,   5.36 +0.j,   5.646+0.j,
                             3.76 +0.j,   3.88 +0.j,   1.456+0.j,   1.468+0.j,   0.98 +0.j,
                            21.6  +0.j, 104.6  +0.j,   1.738+0.j,  38.   +0.j,   2.48 +0.j,
                             2.48 +0.j,   5.78 +0.j,   3.1  +0.j,   2.682+0.j,   7.092+0.j,
                             8.28 +0.j,  12.2  +0.j,  10.198+0.j])
        self.assert_equal(line_r, res_r)
        self.assert_equal(line_x, res_x)
        self.assert_equal(line_h, res_h)

        pp_net = net
        # fix the missing values
        tap_step_pct = 1.0 * pp_net.trafo["tap_step_percent"].values
        tap_step_pct[~np.isfinite(tap_step_pct)] = 0.

        tap_pos = 1.0 * pp_net.trafo["tap_pos"].values
        tap_pos[~np.isfinite(tap_pos)] = 0.

        is_tap_hv_side = pp_net.trafo["tap_side"].values == "hv"
        is_tap_hv_side[~np.isfinite(is_tap_hv_side)] = True

        tap_angles_ = 1.0 * pp_net.trafo["tap_step_degree"].values
        tap_angles_[~np.isfinite(tap_angles_)] = 0.
        tap_angles_ = np.deg2rad(tap_angles_)

        trafo_r, trafo_x, trafo_b = self.converter.get_trafo_param(tap_step_pct,
                                                                   tap_pos,
                                                                   tap_angles_,  # in radian !
                                                                   is_tap_hv_side,
                                                                   pp_net.bus.loc[pp_net.trafo["hv_bus"]]["vn_kv"],
                                                                   pp_net.bus.loc[pp_net.trafo["lv_bus"]]["vn_kv"],
                                                                   pp_net.trafo["vk_percent"].values,
                                                                   pp_net.trafo["vkr_percent"].values,
                                                                   pp_net.trafo["sn_mva"].values,
                                                                   pp_net.trafo["pfe_kw"].values,
                                                                   pp_net.trafo["i0_percent"].values,
                                                                   )
        trafo_r_res = np.array([0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.81494977e-04,
                               3.39887086e-06, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                               0.00000000e+00, 0.00000000e+00, 1.37295648e-05, 0.00000000e+00,
                               1.73571860e-05])
        trafo_x_res = np.array([2.67000000e-04, 3.82000000e-04, 3.70000000e-04, 2.06930358e-03,
                               4.04933224e-05, 3.88000000e-04, 3.75000000e-04, 3.86000000e-04,
                               2.68000000e-04, 3.70000000e-04, 1.59594718e-04, 3.70000000e-04,
                               2.01181945e-04])
        trafo_h_res = np.array([ 0.        -0.j        ,  0.        -0.j        ,
                                0.        -0.j        ,  4.4602909 -0.00140652j,
                               16.40272367-0.00022869j,  0.        -0.j        ,
                                0.        -0.j        ,  0.        -0.j        ,
                                0.        -0.j        ,  0.        -0.j        ,
                               63.96323106-0.01411497j,  0.        -0.j        ,
                               81.1310369 -0.02879733j])
        self.assert_equal(trafo_r, trafo_r_res)
        self.assert_equal(trafo_x, trafo_x_res)
        self.assert_equal(trafo_b, trafo_h_res)


if __name__ == "__main__":
    unittest.main()