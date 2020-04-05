import warnings
import os
import sys
import copy
import numpy as np

import pdb
try:
    # TODO will be deprecated in future version
    from grid2op.Backend import Backend
    # from grid2op.BackendPandaPower import PandaPowerBackend
    from grid2op.Backend import PandaPowerBackend
    from grid2op.Exceptions import InvalidLineStatus, BackendError, DivergingPowerFlow
    grid2op_installed = True
except (ImportError, ModuleNotFoundError) as e:
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

        self.nb_bus_total = None
        self.initdc = True  # does not really hurt computation time
        self.__nb_powerline = None
        self.__nb_bus_before = None
        self._init_bus_load = None
        self._init_bus_gen = None
        self._init_bus_lor = None
        self._init_bus_lex = None
        self._big_topo_to_obj = None
        self.nb_obj_per_bus = None

        self.next_prod_p = None  # this vector is updated with the action that will modify the environment
        # it is done to keep track of the redispatching

        self.topo_vect = None
        self.shunt_topo_vect = None

        self.init_pp_backend = PandaPowerBackend()

        self.V = None
        self.max_it = 10
        self.tol = 1e-8  # tolerance for the solver

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

        self.thermal_limit_a = None

        self._iref_slack = None
        self._id_bus_added = None
        self._fact_mult_gen = -1
        self._what_object_where = None
        self._number_true_line = -1
        self._corresp_name_fun = {}
        self._get_vector_inj = {}
        self.dim_topo = -1

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

        self.__nb_powerline = self.init_pp_backend._grid.line.shape[0]
        self.__nb_bus_before = self.init_pp_backend.get_nb_active_bus()
        self._init_bus_load = 1.0 * self.init_pp_backend._grid.load["bus"].values
        self._init_bus_gen = 1.0 * self.init_pp_backend._grid.gen["bus"].values
        self._init_bus_lor = 1.0 * self.init_pp_backend._grid.line["from_bus"].values
        self._init_bus_lex = 1.0 * self.init_pp_backend._grid.line["to_bus"].values

        t_for = 1.0 * self.init_pp_backend._grid.trafo["hv_bus"].values
        t_fex = 1.0 * self.init_pp_backend._grid.trafo["lv_bus"].values
        self._init_bus_lor = np.concatenate((self._init_bus_lor, t_for)).astype(np.int)
        self._init_bus_lex = np.concatenate((self._init_bus_lex, t_fex)).astype(np.int)
        self._big_topo_to_obj = [(None, None) for _ in range(self.dim_topo)]

        nm_ = "load"
        for load_id, pos_big_topo  in enumerate(self.load_pos_topo_vect):
            self._big_topo_to_obj[pos_big_topo] = (load_id, nm_)
        nm_ = "gen"
        for gen_id, pos_big_topo  in enumerate(self.gen_pos_topo_vect):
            self._big_topo_to_obj[pos_big_topo] = (gen_id, nm_)
        nm_ = "lineor"
        for l_id, pos_big_topo in enumerate(self.line_or_pos_topo_vect):
            self._big_topo_to_obj[pos_big_topo] = (l_id, nm_)
        nm_ = "lineex"
        for l_id, pos_big_topo  in enumerate(self.line_ex_pos_topo_vect):
            self._big_topo_to_obj[pos_big_topo] = (l_id, nm_)

        self.prod_p = 1.0 * self.init_pp_backend._grid.gen["p_mw"].values
        self.next_prod_p = 1.0 * self.init_pp_backend._grid.gen["p_mw"].values

        # for shunts
        self.n_shunt = self.init_pp_backend.n_shunt
        self.shunt_to_subid = self.init_pp_backend.shunt_to_subid
        self.name_shunt = self.init_pp_backend.name_shunt
        self.shunts_data_available = self.init_pp_backend.shunts_data_available

        # number of object per bus, to activate, deactivate them
        self.nb_obj_per_bus = np.zeros(2 * self.__nb_bus_before, dtype=np.int)

        self.topo_vect = np.ones(self.dim_topo, dtype=np.int)
        if self.shunts_data_available:
            self.shunt_topo_vect = np.ones(self.n_shunt, dtype=np.int)

        self._count_object_per_bus()

    def _count_object_per_bus(self):
        # should be called only when self.topo_vect and self.shunt_topo_vect are set
        # todo factor that more properly to update it when it's modified, and not each time

        self.nb_obj_per_bus = np.zeros(2 * self.__nb_bus_before, dtype=np.int)

        arr_ = self.topo_vect[self.load_pos_topo_vect] - 1
        # TODO handle -1 here, eventually
        arr_ = self.load_to_subid + self.__nb_bus_before * arr_
        self.nb_obj_per_bus[arr_] += 1

        arr_ = self.topo_vect[self.gen_pos_topo_vect] - 1
        # TODO handle -1 here, eventually
        arr_ = self.gen_to_subid + self.__nb_bus_before * arr_
        self.nb_obj_per_bus[arr_] += 1

        arr_ = self.topo_vect[self.line_or_pos_topo_vect]
        is_connected = arr_ > 0  # powerline is disconnected
        arr_ = self.line_or_to_subid[is_connected] + self.__nb_bus_before * (arr_[is_connected] -1)
        self.nb_obj_per_bus[arr_] += 1

        arr_ = self.topo_vect[self.line_ex_pos_topo_vect]
        is_connected = arr_ > 0  # powerline is disconnected
        arr_ = self.line_ex_to_subid[is_connected] + self.__nb_bus_before * (arr_[is_connected] -1)
        self.nb_obj_per_bus[arr_] += 1

        if self.shunts_data_available:
            arr_ = self.shunt_topo_vect
            is_connected = arr_ > 0
            arr_ = self.shunt_to_subid[is_connected] + self.__nb_bus_before * (arr_[is_connected] - 1)
            self.nb_obj_per_bus[arr_] += 1

    def _deactivate_unused_bus(self):
        for bus_id, nb in enumerate(self.nb_obj_per_bus):
            if nb == 0:
                self._grid.deactivate_bus(bus_id)
            else:
                self._grid.reactivate_bus(bus_id)

    def close(self):
        self.init_pp_backend.close()
        self._grid = None

    def _convert_id_topo(self, id_big_topo):
        """
        convert an id of the big topo vector into:

        - the id of the object in its "only object" (eg if id_big_topo represents load 2, then it will be 2)
        - the type of object among: "load", "gen", "lineor" and "lineex"

        """
        return self._big_topo_to_obj[id_big_topo]

    def _switch_bus_me(self, tmp):
        """
        return 1 if tmp is 2 else 2 if tmp is one
        """
        return 2 * (1 - tmp) + 2

    def apply_action(self, action):
        # change the _injection if needed
        dict_injection, set_status, switch_status, set_topo_vect, switcth_topo_vect, redispatching, shunts = action()

        if "load_p" in dict_injection:
            tmp = dict_injection["load_p"]
            for i, val in enumerate(tmp):
                if np.isfinite(val):
                    self._grid.change_p_load(i, val)
        if "load_q" in dict_injection:
            tmp = dict_injection["load_q"]
            for i, val in enumerate(tmp):
                if np.isfinite(val):
                    self._grid.change_q_load(i, val)
        if "prod_p" in dict_injection:
            tmp = dict_injection["prod_p"]
            for i, val in enumerate(tmp):
                if np.isfinite(val):
                    self._grid.change_p_gen(i, val)
                    self.next_prod_p[i] = val
        if "prod_v" in dict_injection:
            tmp = dict_injection["prod_v"]
            for i, val in enumerate(tmp):
                if np.isfinite(val):
                    self._grid.change_v_gen(i, val / self.prod_pu_to_kv[i])

        if np.any(redispatching != 0.):
            for i, val in enumerate(redispatching):
                if np.isfinite(val):
                    if val != 0.:
                        self._grid.change_p_gen(i, val + self.next_prod_p[i])
                        self.next_prod_p[i] += val

        # shunts
        if shunts:
            arr_ = shunts["shunt_p"]
            for sh_id, new_p in enumerate(arr_):
                if np.isfinite(new_p):
                    self._grid.change_p_shunt(sh_id, new_p)
            arr_ = shunts["shunt_q"]
            for sh_id, new_q in enumerate(arr_):
                if np.isfinite(new_q):
                    self._grid.change_q_shunt(sh_id, new_q)
            arr_ = shunts["shunt_bus"]
            for sh_id, new_bus in enumerate(arr_):
                if new_bus == -1:
                    self._grid.deactivate_shunt(sh_id)
                    self.shunt_topo_vect[sh_id] = -1
                elif new_bus == 1:
                    self._grid.reactivate_shunt(sh_id)
                    self._grid.change_bus_shunt(sh_id, self.shunt_to_subid[sh_id])
                    self.shunt_topo_vect[sh_id] = 1
                elif new_bus == 2:
                    self._grid.reactivate_shunt(sh_id)
                    self._grid.change_bus_shunt(sh_id, self.shunt_to_subid[sh_id]+self.__nb_bus_before)
                    self.shunt_topo_vect[sh_id] = 2

        # topology
        possiblechange = set_topo_vect != 0
        if np.any(possiblechange) or np.any(switcth_topo_vect):
            actual_topo_full = 1.0 * self.get_topo_vect()
            # new topology vector
            for id_el, do_i_switch in enumerate(switcth_topo_vect):
                if do_i_switch:
                    new_bus_me = self._switch_bus_me(actual_topo_full[id_el])
                    self.topo_vect[id_el] = new_bus_me

            for id_el, new_bus in enumerate(set_topo_vect):
                if new_bus != 0.:
                    self.topo_vect[id_el] = new_bus

            if np.any(self.topo_vect != actual_topo_full):
                # i made at least a real change, so i implement it in the backend
                for id_el, new_bus in enumerate(self.topo_vect):
                    if new_bus < 0:
                        # I don't handle disconnection here
                        continue
                    # if self.topo_vect[id_el] == actual_topo_full[id_el]: continue
                    id_el_backend, type_obj = self._convert_id_topo(id_el)
                    if type_obj == "load":
                        new_bus_backend = self._klu_bus_from_grid2op_bus(new_bus, self._init_bus_load[id_el_backend])
                        self._grid.change_bus_load(id_el_backend, new_bus_backend)
                    elif type_obj == "gen":
                        new_bus_backend = self._klu_bus_from_grid2op_bus(new_bus, self._init_bus_gen[id_el_backend])
                        self._grid.change_bus_gen(id_el_backend, new_bus_backend)
                    elif type_obj == "lineor":
                        new_bus_backend = self._klu_bus_from_grid2op_bus(new_bus, self._init_bus_lor[id_el_backend])
                        if id_el_backend < self.__nb_powerline:
                            # it's a powerline
                            self._grid.change_bus_powerline_or(id_el_backend, new_bus_backend)
                        else:
                            # it's a trafo
                            self._grid.change_bus_trafo_hv(id_el_backend - self.__nb_powerline, new_bus_backend)
                    elif type_obj == "lineex":
                        new_bus_backend = self._klu_bus_from_grid2op_bus(new_bus, self._init_bus_lex[id_el_backend])
                        # if id_el_backend == 0: pdb.set_trace()
                        if id_el_backend < self.__nb_powerline:
                            # it's a powerline
                            self._grid.change_bus_powerline_ex(id_el_backend, new_bus_backend)
                        else:
                            # it's a trafo
                            self._grid.change_bus_trafo_lv(id_el_backend - self.__nb_powerline, new_bus_backend)

        # change line status if needed
        # note that it is a specification that lines status must override buses reconfiguration.
        if np.any(set_status != 0):
            # print("some set_status are non 0")
            for i, el in enumerate(set_status):
                # TODO performance optim here, it can be vectorized
                if el == -1:
                    self._disconnect_line(i)
                elif el == 1:
                    if i < self.__nb_powerline:
                        self._grid.reactivate_powerline(i)
                    else:
                        self._grid.reactivate_trafo(i - self.__nb_powerline)
        else:
            pass
            # print("all_set_status are 0")

        # switch line status if needed
        if np.any(switch_status):
            powerlines_current_status = self.get_line_status()
            for i, el in enumerate(switch_status):
                if el:
                    connected = powerlines_current_status[i]
                    if connected:
                        self._disconnect_line(i)
                    else:
                        # switch a connected powerline -> i reconnect it
                        # but first i need to check
                        bus_or = set_topo_vect[self.line_or_pos_topo_vect[i]]
                        bus_ex = set_topo_vect[self.line_ex_pos_topo_vect[i]]
                        if bus_ex == 0 or bus_or == 0:
                            raise InvalidLineStatus("Line {} was disconnected. The action switched its status, "
                                                    "without providing buses to connect it on both ends.".format(i))

                        # reconnection has then be handled in the topology
                        # and so is the update of topo_vect
                        if i < self.__nb_powerline:
                            self._grid.reactivate_powerline(i)
                        else:
                            self._grid.reactivate_trafo(i - self.__nb_powerline)

    def runpf(self, is_dc=False):
        try:
            self._count_object_per_bus()
            self._deactivate_unused_bus()
            if is_dc:
                raise NotImplementedError("Not fully implemented at the moment.")
                if self.V is None:
                    self.V = np.ones(self.nb_bus_total, dtype=np.complex_)
                self.V = self._grid.dc_pf(self.V, self.max_it, self.tol)
            else:
                if self.V is None:
                    # init from dc approx in this case
                    self.V = np.ones(self.nb_bus_total, dtype=np.complex_) * 1.04

                if self.initdc:
                    V = self._grid.dc_pf(self.V, self.max_it, self.tol)
                    if V.shape[0] == 0:
                        # V = self._grid.ac_pf(self.V, self.max_it, self.tol)
                        raise DivergingPowerFlow("divergence of powerflow (non connected grid)")
                    self.V = V

                V = self._grid.ac_pf(self.V, self.max_it, self.tol)
                if V.shape[0] == 0:
                    # V = self._grid.ac_pf(self.V, self.max_it, self.tol)
                    raise DivergingPowerFlow("divergence of powerflow")
                self.V = V
                self.V[self.V == 0.] = 1.
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
                self.next_prod_p[:] = self.prod_p
                res = True
                # TODO ! below !!! gen_q not handled!!!
                # self.prod_q = 1.0 * self.prod_p
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
            self.load_p = np.full(self.n_load, dtype=np.float, fill_value=np.NaN)
            self.load_q = self.load_p
            self.load_v = self.load_p
            self.prod_p = np.full(self.n_load, dtype=np.float, fill_value=np.NaN)
            self.next_prod_p[:] = self.prod_p
            self.prod_q = self.prod_p
            self.prod_v = self.prod_p
            res = False
        return res

    def copy(self):
        mygrid = self._grid
        self._grid = None
        res = copy.deepcopy(self)
        res._grid = init(self.init_pp_backend._grid)
        #TODO I need a c++ method that would just copy the state of the grid (bus connection, powerlines connected etc.)
        # TODO this could be done in a "get_action_to_set_me" and use to update obsenv for example!
        self._grid = mygrid
        return res

    def get_line_status(self):
        l_s = self._grid.get_lines_status()
        t_s = self._grid.get_trafo_status()
        return np.concatenate((l_s, t_s)).astype(np.bool)

    def get_line_flow(self):
        return self.a_or

    def _grid2op_bus_from_klu_bus(self, klu_bus):
        res = 0
        if klu_bus != 0:
            # object is connected
            res = 1 if klu_bus < self.__nb_bus_before else 2
        return res

    def _klu_bus_from_grid2op_bus(self, grid2op_bus, grid2op_bus_init):
        if grid2op_bus == 0:
            res = grid2op_bus_init
        elif grid2op_bus == 1:
            res = grid2op_bus_init
        elif grid2op_bus == 2:
            res = grid2op_bus_init + self.__nb_bus_before
        else:
            raise BackendError("grid2op bus must be 0 1 or 2")
        return int(res)

    def get_topo_vect(self):
        return self.topo_vect

    def generators_info(self):
        return self.prod_p, self.prod_q, self.prod_v

    def loads_info(self):
        return self.load_p, self.load_q, self.load_v

    def lines_or_info(self):
        return self.p_or, self.q_or, self.v_or, self.a_or

    def lines_ex_info(self):
        return self.p_ex, self.q_ex, self.v_ex, self.a_ex

    def shunt_info(self):
        tmp = self._grid.get_shunts_res()
        shunt_bus = [self._grid.get_bus_shunt(i) for i in range(self.n_shunt)]
        return (tmp[0], tmp[1], tmp[2], shunt_bus)

    def _disconnect_line(self, id_):
        self.topo_vect[self.line_ex_pos_topo_vect[id_]] = -1
        self.topo_vect[self.line_or_pos_topo_vect[id_]] = -1
        if id_ < self.__nb_powerline:
            self._grid.deactivate_powerline(id_)
        else:
            self._grid.deactivate_trafo(id_ - self.__nb_powerline)
