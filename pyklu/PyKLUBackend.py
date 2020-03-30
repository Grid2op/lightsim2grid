import warnings
import os
import sys
import copy
import numpy as np

import pdb
try:
    # TODO will be deprecated in future version
    from grid2op.Backend import Backend
    from grid2op.BackendPandaPower import PandaPowerBackend
    grid2op_installed = True
except (ImportError, ModuleNotFoundError) as e:
    print("Error {}".format(e))
    grid2op_installed = False


    class Backend():
        def __init__(self, detailed_infos_for_cascading_failures=False):
            pass

from pyklu.initGridModel import init


class PyKLUBackend(Backend):
    def __init__(self, detailed_infos_for_cascading_failures=False):
        if not grid2op_installed:
            raise NotImplementedError("Impossible to use a Backend if grid2op is not installed.")
        Backend.__init__(self, detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)

        self.prod_pu_to_kv = None
        self.load_pu_to_kv = None
        self.lines_or_pu_to_kv = None
        self.lines_ex_pu_to_kv = None

        self.p_or = None
        self.q_or = None
        self.v_or = None
        self.a_or = None
        self.p_ex = None
        self.q_ex = None
        self.v_ex = None
        self.a_ex = None

        self.load_p = None
        self.load_q = None
        self.load_v = None

        self.prod_p = None
        self.prod_q = None
        self.prod_v = None

        self._pf_init = "flat"
        self._pf_init = "results"
        self._nb_bus_before = 0

        self.thermal_limit_a = None

        self._iref_slack = None
        self._id_bus_added = None
        self._fact_mult_gen = -1
        self._what_object_where = None
        self._number_true_line = -1
        self._corresp_name_fun = {}
        self._get_vector_inj = {}
        self.dim_topo = -1

        self.init_pp_backend = PandaPowerBackend()

        self.V = None
        self.max_it = 10
        self.tol = 1e-8  # tolerance for the solver

        self.nb_bus_total = None
        self.initdc = False

    def load_grid(self, path=None, filename=None):
        self.init_pp_backend.load_grid(path, filename)
        self._grid = init(self.init_pp_backend._grid)

        self.n_line = self.init_pp_backend.n_line
        self.n_gen = self.init_pp_backend.n_gen
        self.n_load = self.init_pp_backend.n_load
        self.n_sub = self.init_pp_backend.n_sub
        self.sub_info = self.init_pp_backend.sub_info
        self.dim_topo = self.init_pp_backend.dim_topo
        self.load_to_subid = self.init_pp_backend.load_to_subid
        self.gen_to_subid = self.init_pp_backend.gen_to_subid
        self.line_or_to_subid = self.init_pp_backend.line_or_to_subid
        self.line_ex_to_subid = self.init_pp_backend.line_ex_to_subid
        self.load_to_sub_pos = self.init_pp_backend.load_to_sub_pos
        self.gen_to_sub_pos = self.init_pp_backend.gen_to_sub_pos
        self.line_or_to_sub_pos = self.init_pp_backend.line_or_to_sub_pos
        self.line_ex_to_sub_pos = self.init_pp_backend.line_ex_to_sub_pos

        self.prod_pu_to_kv = self.init_pp_backend.prod_pu_to_kv
        self.load_pu_to_kv = self.init_pp_backend.load_pu_to_kv
        self.lines_or_pu_to_kv = self.init_pp_backend.lines_or_pu_to_kv
        self.lines_ex_pu_to_kv = self.init_pp_backend.lines_ex_pu_to_kv

        self.name_gen = self.init_pp_backend.name_gen
        self.name_load = self.init_pp_backend.name_load
        self.name_line = self.init_pp_backend.name_line
        self.name_sub = self.init_pp_backend.name_sub
        self._compute_pos_big_topo()
        self.nb_bus_total = self.init_pp_backend._grid.bus.shape[0]

        self.thermal_limit_a = self.init_pp_backend.thermal_limit_a

        # deactive the buses that have been added
        nb_bus_init = self.init_pp_backend._grid.bus.shape[0] // 2
        for i in range(nb_bus_init):
            self._grid.deactivate_bus(i + nb_bus_init)

    def close(self):
        self.init_pp_backend.close()
        self._grid = None

    def apply_action(self, action):
        # change the _injection if needed
        dict_injection, set_status, switch_status, set_topo_vect, switcth_topo_vect, redispatching = action()

        if "load_p" in dict_injection:
            tmp = dict_injection["load_p"]
            for i, val in enumerate(tmp):
                if np.isfinite(val):
                    pass
                    self._grid.change_p_load(i, val)
        if "load_q" in dict_injection:
            tmp = dict_injection["load_q"]
            for i, val in enumerate(tmp):
                if np.isfinite(val):
                    pass
                    self._grid.change_q_load(i, val)
        if "prod_p" in dict_injection:
            tmp = dict_injection["prod_p"]
            for i, val in enumerate(tmp):
                if np.isfinite(val):
                    pass
                    self._grid.change_p_gen(i, val)

        if "prod_v" in dict_injection:
            tmp = dict_injection["prod_v"]
            for i, val in enumerate(tmp):
                if np.isfinite(val):
                    pass
                    self._grid.change_v_gen(i, val / self.prod_pu_to_kv[i])

        # TODO for the rest !!!

        pass

    def runpf(self, is_dc=False):
        try:
            if is_dc:
                if self.V is None:
                    self.V = np.ones(self.nb_bus_total, dtype=np.complex_)
                self.V = self._grid.dc_pf(self.V, self.max_it, self.tol)
                raise NotImplementedError("DC is not implemented at the moment")
            else:
                if self.V is None:
                    # init from dc approx in this case
                    self.V = np.ones(self.nb_bus_total, dtype=np.complex_) * 1.04
                if self.initdc:
                    self.V = self._grid.dc_pf(self.V, self.max_it, self.tol)
                V = self._grid.ac_pf(self.V, self.max_it, self.tol)
                if V.shape[0] == 0:
                    raise RuntimeError("divergence of powerflow")
                # self.V = V
                lpor, lqor, lvor, laor = self._grid.get_lineor_res()
                lpex, lqex, lvex, laex = self._grid.get_lineex_res()
                tpor, tqor, tvor, taor = self._grid.get_trafohv_res()
                tpex, tqex, tvex, taex = self._grid.get_trafolv_res()

                self.p_or = np.concatenate((lpor, tpor))
                self.q_or = np.concatenate((lqor, tqor))
                self.v_or = np.concatenate((lvor, tvor))
                self.a_or = np.concatenate((laor, taor))
                self.a_or *= 1000
                self.p_ex = np.concatenate((lpex, tpex))
                self.q_ex = np.concatenate((lqex, tqex))
                self.v_ex = np.concatenate((lvex, tvex))
                self.a_ex = np.concatenate((laex, taex))
                self.a_ex *= 1000

                self.load_p, self.load_q, self.load_v = self._grid.get_loads_res()
                self.prod_p, self.prod_q, self.prod_v = self._grid.get_gen_res()
                # TODO ! below !!!
                self.prod_q = 1.0 * self.prod_p

        except Exception as e:
            # of the powerflow has not converged, results are Nan
            self.p_or = np.full(self.n_line, dtype=np.float, fill_value=np.NaN)
            self.q_or = self.p_or
            self.v_or = self.p_or
            self.a_or = self.p_or
            self.p_ex = self.p_or
            self.q_ex = self.p_or
            self.v_ex = self.p_or
            self.a_ex = self.p_or
            self._nb_bus_before = None
            self.load_p = np.full(self.n_load, dtype=np.float, fill_value=np.NaN)
            self.load_q = self.load_p
            self.load_v = self.load_p
            self.prod_p = np.full(self.n_load, dtype=np.float, fill_value=np.NaN)
            self.prod_q = self.prod_p
            self.prod_v = self.prod_p
            return False
        return True

    def copy(self):
        mygrid = self._grid
        self._grid = None
        res = copy.deepcopy(self)
        res._grid = init(self.init_pp_backend._grid)
        self._grid = mygrid
        return res

    def get_line_status(self):
        l_s = self._grid.get_lines_status()
        t_s = self._grid.get_trafo_status()
        return np.concatenate((l_s, t_s))

    def get_line_flow(self):
        return self.a_or

    def get_topo_vect(self):
        # TODO
        return np.ones(self.dim_topo, dtype=np.int)

    def generators_info(self):
        return self.prod_p, self.prod_q, self.prod_v

    def loads_info(self):
        return self.load_p, self.load_q, self.load_v

    def lines_or_info(self):
        return self.p_or, self.q_or, self.v_or, self.a_or

    def lines_ex_info(self):
        return self.p_ex, self.q_ex, self.v_ex, self.a_ex

    def shunt_info(self):
        return self._grid.get_shunts_res()

    def _disconnect_line(self, id):
        self._grid.deactivate_powerline(id)