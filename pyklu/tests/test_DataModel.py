import os
import unittest
import numpy as np
import pdb
from scipy import sparse
from pyklu_cpp import DataModel, PandaPowerConverter
import pandapower.networks as pn
import pandapower as pp


class MakeTests(unittest.TestCase):
    def setUp(self):
        self.net_ref = pn.case118()
        self.net_datamodel = pn.case118()

        # compute a powerflow on a net without anything
        pp.runpp(self.net_ref)

        # initialize constant stuff
        self.max_it = 10
        self.tol = 1e-8  # tolerance for the solver
        self.tol_test = 1e-4  # tolerance for the test (2 matrices are equal if the l_1 of their difference is less than this)

        # initialize and use converters
        self.converter = PandaPowerConverter()
        self.converter.set_sn_mva(self.net_datamodel.sn_mva)  # TODO raise an error if not set !
        self.converter.set_f_hz(self.net_datamodel.f_hz)
        self.line_r, self.line_x, self.line_h = \
            self.converter.get_line_param(
                self.net_datamodel.line["r_ohm_per_km"].values * self.net_datamodel.line["length_km"].values,
                self.net_datamodel.line["x_ohm_per_km"].values * self.net_datamodel.line["length_km"].values,
                self.net_datamodel.line["c_nf_per_km"].values * self.net_datamodel.line["length_km"].values,
                self.net_datamodel.line["g_us_per_km"].values * self.net_datamodel.line["length_km"].values,
                self.net_datamodel.bus.loc[self.net_datamodel.line["from_bus"]]["vn_kv"],
                self.net_datamodel.bus.loc[self.net_datamodel.line["to_bus"]]["vn_kv"]
            )
        self.trafo_r, self.trafo_x, self.trafo_b = \
            self.converter.get_trafo_param(self.net_datamodel.trafo["vn_hv_kv"].values,
                                           self.net_datamodel.trafo["vn_lv_kv"].values,
                                           self.net_datamodel.trafo["vk_percent"].values,
                                           self.net_datamodel.trafo["vkr_percent"].values,
                                           self.net_datamodel.trafo["sn_mva"].values,
                                           self.net_datamodel.trafo["pfe_kw"].values,
                                           self.net_datamodel.trafo["i0_percent"].values,
                                           self.net_datamodel.bus.loc[self.net_datamodel.trafo["lv_bus"]]["vn_kv"]
                                           )

        # set up the data model accordingly
        self.model = DataModel()
        tmp_bus_ind = np.argsort(self.net_datamodel.bus.index)
        self.model.init_bus(self.net_datamodel.bus.iloc[tmp_bus_ind]["vn_kv"].values,
                            self.net_datamodel.line.shape[0],
                            self.net_datamodel.trafo.shape[0])

        self.model.init_powerlines(self.line_r, self.line_x, self.line_h,
                                   self.net_datamodel.line["from_bus"].values,
                                   self.net_datamodel.line["to_bus"].values
                                   )

        # init the shunts
        self.model.init_shunt(self.net_datamodel.shunt["p_mw"].values,
                              self.net_datamodel.shunt["q_mvar"].values,
                              self.net_datamodel.shunt["bus"].values
                              )

        tap_step_pct = self.net_datamodel.trafo["tap_step_percent"].values
        tap_step_pct[~np.isfinite(tap_step_pct)] = 0.

        tap_pos = self.net_datamodel.trafo["tap_pos"].values
        tap_pos[~np.isfinite(tap_pos)] = 0.

        is_tap_hv_side = self.net_datamodel.trafo["tap_side"].values == "hv"
        is_tap_hv_side[~np.isfinite(tap_pos)] = True
        self.model.init_trafo(self.trafo_r,
                              self.trafo_x,
                              self.trafo_b,
                              tap_step_pct,
                              tap_pos,
                              is_tap_hv_side,
                              self.net_datamodel.trafo["hv_bus"].values,
                              self.net_datamodel.trafo["lv_bus"].values)

        self.model.init_loads(self.net_datamodel.load["p_mw"].values,
                              self.net_datamodel.load["q_mvar"].values,
                              self.net_datamodel.load["bus"].values
                              )
        self.model.init_generators(self.net_datamodel.gen["p_mw"].values,
                                   self.net_datamodel.gen["vm_pu"].values,
                                   self.net_datamodel.gen["bus"].values
                                   )

        # TODO handle that better maybe
        self.model.add_slackbus(self.net_datamodel.ext_grid["bus"].values[0])

        # this should be the vector use as a starting point of the ac powerflow
        self.V0 = np.array([0.92371316+0.24244381j, 0.95434983+0.26251308j,
                           0.95305717+0.26716798j, 0.94099945+0.33245156j,
                           0.93050132+0.3374374j , 0.94554411+0.29333656j,
                           0.94747004+0.2863516j , 0.91976137+0.42925986j,
                           0.83121514+0.53738077j, 0.79021868+0.69141481j,
                           0.94703065+0.28780142j, 0.94922593+0.28119411j,
                           0.95435892+0.26248003j, 0.95289026+0.26776271j,
                           0.93678483+0.25166284j, 0.95116199+0.2738386j ,
                           0.94299095+0.30077361j, 0.93833685+0.2573965j ,
                           0.93034369+0.2447542j , 0.95398128+0.26384925j,
                           0.94725669+0.28705657j, 0.9345412 +0.32608199j,
                           0.90393663+0.40323105j, 0.91010681+0.39467656j,
                           0.89230065+0.55344336j, 0.84467665+0.56280223j,
                           0.91638352+0.31187377j, 0.94562819+0.29237653j,
                           0.94991823+0.27812238j, 0.91392466+0.38005054j,
                           0.92731294+0.27418919j, 0.91487454+0.30062198j,
                           0.95935742+0.24357759j, 0.9517411 +0.24988972j,
                           0.9593314 +0.24368006j, 0.94985607+0.24119171j,
                           0.95510399+0.2597558j , 0.92958567+0.33995173j,
                           0.96901644+0.20175195j, 0.95305256+0.18052927j,
                           0.97415369+0.17527494j, 0.96487073+0.19811479j,
                           0.9591693 +0.24431734j, 0.9504645 +0.27624978j,
                           0.9427526 +0.30151989j, 0.94220104+0.34968872j,
                           0.91500296+0.37744708j, 0.91811695+0.36980801j,
                           0.94394221+0.39949732j, 0.92562893+0.35058208j,
                           0.94138116+0.3057748j , 0.94629095+0.29022429j,
                           0.9509182 +0.27468398j, 0.91317677+0.27952314j,
                           0.91164381+0.27424363j, 0.91268579+0.27770605j,
                           0.94061208+0.30813248j, 0.94525377+0.29358478j,
                           0.91896615+0.35457892j, 0.89909147+0.41392179j,
                           0.89739679+0.42977202j, 0.90439526+0.42198723j,
                           0.90192434+0.40771215j, 0.88911864+0.43493075j,
                           0.87693705+0.49092403j, 0.91352982+0.51765169j,
                           0.88530602+0.44263976j, 0.86830043+0.4751327j ,
                           0.89633629+0.5175j    , 0.90385573+0.38897405j,
                           0.91063484+0.38786711j, 0.90471675+0.37667971j,
                           0.91279673+0.38585378j, 0.88787865+0.35977145j,
                           0.9093251 +0.39092784j, 0.87330985+0.35577929j,
                           0.89103286+0.46700797j, 0.87908727+0.45486512j,
                           0.87589698+0.46097851j, 0.89625694+0.52756373j,
                           0.86276979+0.48510307j, 0.86964118+0.47267423j,
                           0.85945866+0.49094554j, 0.83840965+0.52608551j,
                           0.81982782+0.54599208j, 0.83511119+0.53130594j,
                           0.85183707+0.55190453j, 0.79180803+0.59391645j,
                           0.75766148+0.66028333j, 0.80798425+0.56337061j,
                           0.80514315+0.55869895j, 0.80807875+0.5719342j ,
                           0.8365277 +0.52907288j, 0.85506824+0.49855292j,
                           0.86493788+0.48122674j, 0.86609865+0.47913447j,
                           0.86312091+0.48447807j, 0.86529541+0.48058356j,
                           0.88467666+0.48728554j, 0.87935138+0.5109111j ,
                           0.84499803+0.51543675j, 0.82171572+0.55179704j,
                           0.90204736+0.45432429j, 0.88815266+0.39246127j,
                           0.88958415+0.37398535j, 0.91423777+0.37929673j,
                           0.89432564+0.3263214j , 0.9188473 +0.3679896j ,
                           0.92118576+0.36209598j, 0.90936353+0.34610254j,
                           0.90609384+0.37335499j, 0.92545851+0.30684125j,
                           0.9461778 +0.3013247j , 0.94194179+0.30404337j,
                           0.94196023+0.30398625j, 0.88520788+0.47584873j,
                           0.95652941+0.25445667j, 0.91599499+0.37503317j])

    def assert_equal(self, tmp, ref):
        assert np.all(tmp.shape == ref.shape), "vector does not have the same shape"
        assert np.max(np.abs(tmp - ref)) <= self.tol_test
        assert np.mean(np.abs(tmp - ref)) <= self.tol_test

    def check_res(self, net):
        # check lines
        por, qor, vor, aor = self.model.get_lineor_res()
        self.assert_equal(por, net.res_line["p_from_mw"].values)
        self.assert_equal(qor, net.res_line["q_from_mvar"].values)
        self.assert_equal(aor, net.res_line["i_from_ka"].values)
        vor_pp = net.bus.loc[net.line["from_bus"].values]["vn_kv"].values * net.res_line["vm_from_pu"].values
        self.assert_equal(vor, vor_pp)

        # check trafo
        plv, qlv, vlv, alv = self.model.get_trafolv_res()
        self.assert_equal(plv, net.res_trafo["p_lv_mw"].values)
        self.assert_equal(qlv, net.res_trafo["q_lv_mvar"].values)
        self.assert_equal(alv, net.res_trafo["i_lv_ka"].values)
        vlv_pp = net.bus.loc[net.trafo["lv_bus"].values]["vn_kv"].values * net.res_trafo["vm_lv_pu"].values
        self.assert_equal(vlv, vlv_pp)

        # check loads
        load_p, load_q, load_v = self.model.get_loads_res()
        self.assert_equal(load_p, net.res_load["p_mw"].values)
        self.assert_equal(load_q, net.res_load["q_mvar"].values)

        # check shunts
        shunt_p, shunt_q, shunt_v = self.model.get_shunts_res()
        self.assert_equal(shunt_p, net.res_shunt["p_mw"].values)
        self.assert_equal(shunt_q, net.res_shunt["q_mvar"].values)

        # check generators
        prod_p, prod_q, prod_v = self.model.get_gen_res()
        self.assert_equal(prod_p, net.res_gen["p_mw"].values)
        # self.assert_equal(prod_q, net.res_gen["q_mvar"].values)
        v_gen_pp = net.bus.loc[net.gen["bus"].values]["vn_kv"].values * net.res_gen["vm_pu"].values
        self.assert_equal(prod_v, v_gen_pp)

    def test_acpf(self):
        """
        Reference without modifying anything
        """
        has_conv = self.model.compute_newton(self.V0, self.max_it, self.tol)
        assert has_conv, "powerflow diverged !"

        self.check_res(self.net_ref)


if __name__ == "__main__":
    unittest.main()