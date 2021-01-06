import unittest
import numpy as np
import pandapower.networks as pn
import pandapower as pp

from lightsim2grid.initGridModel import init
import warnings
import pdb


class BaseTests:
    def setUp(self):
        self.net_ref = pn.case118()
        self.net_datamodel = pn.case118()
        pp.runpp(self.net_datamodel)

        # initialize constant stuff
        self.max_it = 10
        self.tol = 1e-8  # tolerance for the solver
        self.tol_test = 1e-5  # tolerance for the test (2 matrices are equal if the l_1 of their difference is less than this)

        # initialize and use converters
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.model = init(self.net_datamodel)

    def assert_equal(self, tmp, ref, error=""):
        assert np.all(tmp.shape == ref.shape), "vector does not have the same shape"
        assert np.max(np.abs(tmp - ref)) <= self.tol_test, error
        assert np.mean(np.abs(tmp - ref)) <= self.tol_test, error

    def check_res(self, Vfinal, net):
        assert Vfinal.shape[0] > 0, "powerflow diverged !"

        # check lines
        l_is = self.net_ref.line["in_service"]
        por, qor, vor, aor = self.model.get_lineor_res()
        self.assert_equal(por, net.res_line["p_from_mw"].values, "error for p_from")
        self.assert_equal(qor, net.res_line["q_from_mvar"].values, "error for q_from")
        self.assert_equal(aor, net.res_line["i_from_ka"].values, "error for i_from_ka")
        vor_pp = net.bus.loc[net.line["from_bus"].values]["vn_kv"].values * net.res_line["vm_from_pu"].values
        self.assert_equal(vor[l_is], vor_pp[l_is], "error for vor_pp")

        # check trafo
        f_is = self.net_ref.trafo["in_service"]
        plv, qlv, vlv, alv = self.model.get_trafolv_res()
        self.assert_equal(plv, net.res_trafo["p_lv_mw"].values, "error for p_lv_mw")
        self.assert_equal(qlv, net.res_trafo["q_lv_mvar"].values, "error for q_lv_mvar")
        self.assert_equal(alv, net.res_trafo["i_lv_ka"].values, "error for i_lv_ka")
        vlv_pp = net.bus.loc[net.trafo["lv_bus"].values]["vn_kv"].values * net.res_trafo["vm_lv_pu"].values
        self.assert_equal(vlv[f_is], vlv_pp[f_is], "error for vlv_pp")

        # check loads
        l_is = self.net_ref.load["in_service"]
        load_p, load_q, load_v = self.model.get_loads_res()
        self.assert_equal(load_p[l_is], net.res_load["p_mw"].values[l_is], "error for load p_mw")
        self.assert_equal(load_q[l_is], net.res_load["q_mvar"].values[l_is], "error for load q_mvar")

        # check shunts
        s_is = self.net_ref.shunt["in_service"]
        shunt_p, shunt_q, shunt_v = self.model.get_shunts_res()
        self.assert_equal(shunt_p[s_is], net.res_shunt["p_mw"].values[s_is], "error for shunt p_mw")
        self.assert_equal(shunt_q[s_is], net.res_shunt["q_mvar"].values[s_is], "error for shunt q_mvar")

        # check generators
        g_is = self.net_ref.gen["in_service"]
        prod_p, prod_q, prod_v = self.model.get_gen_res()
        # test the slack bus is properly modeled as a generator
        assert np.abs(np.sum(net._ppc["gen"][:, 1]) - np.sum(prod_p)) <= self.tol
        if len(prod_p) != g_is.shape[0]:
            # it means a generator has been added for the slack bus
            prod_p = prod_p[:-1]
            prod_q = prod_q[:-1]
            prod_v = prod_v[:-1]
        self.assert_equal(prod_p[g_is], net.res_gen["p_mw"].values[g_is], "error for gen p_mw")
        self.assert_equal(prod_q, net.res_gen["q_mvar"].values, "error for gen q_mvar")
        v_gen_pp = net.bus.loc[net.gen["bus"].values]["vn_kv"].values * net.res_gen["vm_pu"].values
        self.assert_equal(prod_v[g_is], v_gen_pp[g_is], "error for prod_v")

    def make_v0(self, net):
        V0 = np.full(net.bus.shape[0],
                     fill_value=1.0,
                     dtype=np.complex_)
        all_gen_conn = net.gen["bus"].values
        g_is = net.gen["in_service"].values
        V0[all_gen_conn[g_is]] *= self.net_datamodel.gen.iloc[g_is]["vm_pu"].values
        V0[net.ext_grid["bus"].values] = net.ext_grid["vm_pu"].values * np.exp(
            1j * net.ext_grid["va_degree"].values / 360. * 2 * np.pi)
        return V0

    def run_me_pf(self, V0):
        return self.model.compute_newton(V0, self.max_it, self.tol)

    def run_ref_pf(self, net):
        pp.runpp(net, init="flat")

    def do_i_skip(self, func_name):
        # self.skipTest("dev")
        pass

    def _run_both_pf(self, net):
        V0 = self.make_v0(net)
        self.run_ref_pf(net)
        Vfinal = self.run_me_pf(V0)
        return Vfinal

    def test_function_works(self):
        self.do_i_skip("test_function_works")
        self.model.get_loads_status()
        self.model.get_shunts_status()
        self.model.get_gen_status()
        self.model.get_lines_status()
        self.model.get_trafo_status()

    def test_deactivate_index_out_of_bound(self):
        self.do_i_skip("test_deactivate_index_out_of_bound")
        with self.assertRaises(IndexError):
            self.model.deactivate_load(self.net_datamodel.load.shape[0])
        with self.assertRaises(IndexError):
            # +1 is added for gen because in some cases, gen is added in gridmodel for the slack bus
            self.model.deactivate_gen(self.net_datamodel.gen.shape[0]+1)
        with self.assertRaises(IndexError):
            self.model.deactivate_trafo(self.net_datamodel.trafo.shape[0])
        with self.assertRaises(IndexError):
            self.model.deactivate_powerline(self.net_datamodel.line.shape[0])
        with self.assertRaises(IndexError):
            self.model.deactivate_shunt(self.net_datamodel.shunt.shape[0])

    def test_changebus_index_out_of_bound(self):
        self.do_i_skip("test_changebus_index_out_of_bound")
        with self.assertRaises(IndexError):
            self.model.change_bus_load(self.net_datamodel.load.shape[0], 1)
        with self.assertRaises(IndexError):
            # +1 is added for gen because in some cases, gen is added in gridmodel for the slack bus
            self.model.change_bus_gen(self.net_datamodel.gen.shape[0]+1, 1)
        with self.assertRaises(IndexError):
            self.model.change_bus_shunt(self.net_datamodel.shunt.shape[0], 1)
        with self.assertRaises(IndexError):
            self.model.change_bus_powerline_or(self.net_datamodel.line.shape[0], 1)
        with self.assertRaises(IndexError):
            self.model.change_bus_powerline_ex(self.net_datamodel.line.shape[0], 1)
        with self.assertRaises(IndexError):
            self.model.change_bus_trafo_hv(self.net_datamodel.trafo.shape[0], 1)
        with self.assertRaises(IndexError):
            self.model.change_bus_trafo_lv(self.net_datamodel.trafo.shape[0], 1)

    def test_changebus_newbus_out_of_bound(self):
        self.do_i_skip("test_changebus_newbus_out_of_bound")
        newbusid = self.net_datamodel.bus.shape[0]
        with self.assertRaises(IndexError):
            self.model.change_bus_load(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus_gen(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus_shunt(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus_powerline_or(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus_powerline_ex(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus_trafo_hv(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus_trafo_lv(0, newbusid)

    def test_changesetpoint_out_of_bound(self):
        self.do_i_skip("test_changesetpoint_out_of_bound")
        with self.assertRaises(IndexError):
            self.model.change_p_load(self.net_datamodel.load.shape[0], 1)
        with self.assertRaises(IndexError):
            self.model.change_q_load(self.net_datamodel.load.shape[0], 1)
        with self.assertRaises(IndexError):
            # +1 is added for gen because in some cases, gen is added in gridmodel for the slack bus
            self.model.change_p_gen(self.net_datamodel.gen.shape[0] + 1, 1)
        with self.assertRaises(IndexError):
            # +1 is added for gen because in some cases, gen is added in gridmodel for the slack bus
            self.model.change_v_gen(self.net_datamodel.gen.shape[0] + 1, 1)
        with self.assertRaises(IndexError):
            self.model.change_p_shunt(self.net_datamodel.shunt.shape[0], 1)
        with self.assertRaises(IndexError):
            self.model.change_p_shunt(self.net_datamodel.shunt.shape[0], 1)

    def test_pf(self):
        """
        Reference without modifying anything
        """
        self.do_i_skip("test_pf")
        # compute a powerflow on a net without anything
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_disco_gen(self):
        self.do_i_skip("test_pf_disco_gen")
        self.net_ref.gen["in_service"][0] = False
        self.model.deactivate_gen(0)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_disco_load(self):
        self.do_i_skip("test_pf_disco_load")
        self.net_ref.load["in_service"][0] = False
        self.model.deactivate_load(0)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_disco_line(self):
        self.do_i_skip("test_pf_disco_line")
        self.net_ref.line["in_service"][0] = False
        self.model.deactivate_powerline(0)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_disco_shunt(self):
        self.do_i_skip("test_pf_disco_shunt")
        self.net_ref.shunt["in_service"][0] = False
        self.model.deactivate_shunt(0)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_disco_trafo(self):
        self.do_i_skip("test_pf_disco_trafo")
        self.net_ref.trafo["in_service"][0] = False
        self.model.deactivate_trafo(0)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_reactivate(self):
        # i deactivate everything, run a powerflow, and check that reactivating everything and supposes that the results
        # is the same
        self.do_i_skip("test_reactivate")
        self.run_ref_pf(self.net_ref)
        V0 = self.make_v0(self.net_ref)

        # i disconnect a load, the reconnect it
        self.model.deactivate_load(0)
        Vfinal = self.run_me_pf(V0)
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        self.model.reactivate_load(0)
        Vfinal = self.run_me_pf(V0)
        self.check_res(Vfinal, self.net_ref)

        self.model.deactivate_gen(0)
        Vfinal = self.run_me_pf(V0)
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        self.model.reactivate_gen(0)
        Vfinal = self.run_me_pf(V0)
        self.check_res(Vfinal, self.net_ref)

        self.model.deactivate_powerline(0)
        Vfinal = self.run_me_pf(V0)
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        self.model.reactivate_powerline(0)
        Vfinal = self.run_me_pf(V0)
        self.check_res(Vfinal, self.net_ref)

        self.model.deactivate_trafo(0)
        Vfinal = self.run_me_pf(V0)
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        self.model.reactivate_trafo(0)
        Vfinal = self.run_me_pf(V0)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_gen(self):
        self.do_i_skip("test_pf_changebus_gen")
        self.net_ref.gen["bus"][0] = 2
        self.model.change_bus_gen(0, 2)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_load(self):
        self.do_i_skip("test_pf_changebus_load")
        self.net_ref.load["bus"][0] = 2
        self.model.change_bus_load(0, 2)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_shunt(self):
        self.do_i_skip("test_pf_changebus_shunt")
        self.net_ref.shunt["bus"][0] = 2
        self.model.change_bus_shunt(0, 2)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_lineor(self):
        self.do_i_skip("test_pf_changebus_lineor")
        self.net_ref.line["from_bus"][0] = 2
        self.model.change_bus_powerline_or(0, 2)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_lineex(self):
        self.do_i_skip("test_pf_changebus_lineex")
        self.net_ref.line["to_bus"][0] = 2
        self.model.change_bus_powerline_ex(0, 2)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_trafolv(self):
        self.do_i_skip("test_pf_changebus_trafolv")
        self.net_ref.trafo["lv_bus"][0] = 5  # was 4 initially, and 4 is connected to 5
        self.model.change_bus_trafo_lv(0, 5)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_trafohv(self):
        self.do_i_skip("test_pf_changebus_trafohv")
        self.net_ref.trafo["hv_bus"][0] = 29  # was 7 initially, and 7 is connected to 29
        self.model.change_bus_trafo_hv(0, 29)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeloadp(self):
        self.do_i_skip("test_pf_changeloadp")
        self.net_ref.load["p_mw"][0] = 50
        self.model.change_p_load(0, 50)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeloadq(self):
        self.do_i_skip("test_pf_changeloadq")
        self.net_ref.load["q_mvar"][0] = 50
        self.model.change_q_load(0, 50)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeprodp(self):
        self.do_i_skip("test_pf_changeprodp")
        self.net_ref.gen["p_mw"][0] = 50
        self.model.change_p_gen(0, 50)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeprodv(self):
        self.do_i_skip("test_pf_changeprodv")
        self.net_ref.gen["vm_pu"][0] = 1.06
        self.model.change_v_gen(0, 1.06)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeshuntp(self):
        self.skipTest("not usefull but not working at the moment")
        self.do_i_skip("test_pf_changeshuntp")
        self.net_ref.shunt["p_mw"][0] = 10
        self.model.change_p_shunt(0, 10)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeshuntq(self):
        self.do_i_skip("test_pf_changeshuntq")
        self.net_ref.shunt["q_mvar"][0] = 10
        self.model.change_q_shunt(0, 10)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)


class MakeDCTests(BaseTests, unittest.TestCase):
    def run_me_pf(self, V0):
        return self.model.dc_pf(V0, self.max_it, self.tol)

    def run_ref_pf(self, net):
        pp.rundcpp(net, init="flat")

    def do_i_skip(self, test_nm):
        #self.skipTest("dev")
        pass
        # if test_nm == "test_pf":
        #    pass
        #else:
        #    self.skipTest("dev")

    def check_res(self, Vfinal, net):
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        tmp_bus_ind = np.argsort(net.bus.index)
        va_deg = net.res_bus["va_degree"].values
        # vm_pu = net.res_bus["vm_pu"].values
        # pdb.set_trace()
        self.assert_equal(np.angle(Vfinal), va_deg[tmp_bus_ind] / 180. * np.pi)
        # self.assert_equal(np.abs(Vfinal), vm_pu[tmp_bus_ind])


class MakeACTests(BaseTests, unittest.TestCase):
    def run_me_pf(self, V0):
        return self.model.ac_pf(V0, self.max_it, self.tol)

    def run_ref_pf(self, net):
        pp.runpp(net, init="flat")

    def do_i_skip(self, test_nm):
        pass
        # if test_nm == "test_pf":
        #    pass
        # else:
        #    self.skipTest("dev")


if __name__ == "__main__":
    unittest.main()