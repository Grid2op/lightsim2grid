# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import copy
import warnings
import numpy as np

from grid2op.Action import CompleteAction
from grid2op.Backend import Backend
from grid2op.Exceptions import InvalidLineStatus, BackendError, DivergingPowerFlow
from grid2op.Action._BackendAction import _BackendAction
from grid2op.dtypes import dt_float, dt_int, dt_bool

from lightsim2grid.initGridModel import init, SolverType


class LightSimBackend(Backend):
    def __init__(self, detailed_infos_for_cascading_failures=False):
        Backend.__init__(self,
                         detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)

        # lazy loading because it crashes...
        from grid2op.Backend import PandaPowerBackend
        from grid2op.Space import GridObjects  # lazy import
        self.__has_storage = hasattr(GridObjects, "n_storage")
        if not self.__has_storage:
            pass
            # warnings.warn("Please upgrade your grid2Op to >= 1.5.0. You are using a backward compatibility "
            #              "feature that will be removed in further lightsim2grid version.")

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

        self.storage_p = None
        self.storage_q = None
        self.storage_v = None

        self.thermal_limit_a = None

        self._iref_slack = None
        self._id_bus_added = None
        self._fact_mult_gen = -1
        self._what_object_where = None
        self._number_true_line = -1
        self._corresp_name_fun = {}
        self._get_vector_inj = {}
        self.dim_topo = -1
        self._init_action_to_set = None
        self._backend_action_class = None
        self.cst_1 = dt_float(1.0)
        self.__me_at_init = None
        self.__init_topo_vect = None

        # available solver in lightsim
        self.available_solvers = []
        self.comp_time = 0.  # computation time of just the powerflow
        self.__current_solver_type = None

        # hack for the storage unit:
        # in grid2op, for simplicity, I suppose that if a storage is alone on a busbar, and
        # that it produces / absorbs nothing, then that's fine
        # this behaviour in lightsim (c++ side) would be detected as a non connex grid and raise
        # a diverging powerflow
        # i "fake" to disconnect storage with these properties
        # TODO hummm we need to clarify that ! pandapower automatically disconnect this stuff  too ! This is super weird
        # TODO and should rather be handled in pandapower backend
        # backend SHOULD not do these kind of stuff
        self._idx_hack_storage = []

    def get_theta(self):
        """

        Returns
        -------
        line_or_theta: ``numpy.ndarray``
            For each orgin side of powerline, gives the voltage angle
        line_ex_theta: ``numpy.ndarray``
            For each extremity side of powerline, gives the voltage angle
        load_theta: ``numpy.ndarray``
            Gives the voltage angle to the bus at which each load is connected
        gen_theta: ``numpy.ndarray``
            Gives the voltage angle to the bus at which each generator is connected
        storage_theta: ``numpy.ndarray``
            Gives the voltage angle to the bus at which each storage unit is connected
        """
        line_or_theta = np.concatenate((self._grid.get_lineor_theta(), self._grid.get_trafohv_theta()))
        line_ex_theta = np.concatenate((self._grid.get_lineex_theta(), self._grid.get_trafolv_theta()))
        load_theta = self.cst_1 * self._grid.get_load_theta()
        gen_theta = self.cst_1 * self._grid.get_gen_theta()
        storage_theta = self.cst_1 * self._grid.get_storage_theta()
        return line_or_theta, line_ex_theta, load_theta, gen_theta, storage_theta

    def set_solver_type(self, solver_type):
        """
        Change the type of solver you want to use.

        Note that a powergrid should have been loaded for this function to work.

        This function does not modify :attr:`LightSimBackend.max_iter` nor :attr:`LightSimBackend.tol`. You might want
        to modify these values depending on the solver you are using.

        Notes
        ------
        By default, the fastest AC solver is used for your platform. This means that if KLU is available, then it is used
        otherwise it's SparseLU.

        This has to be set for every backend that you want to use. For example, you have to set it
        in the backend of the `_obs_env` of the observation and if you are using "grid2op.MultMixEnv` you
        have to set it in all mixes!

        Parameters
        ----------
        solver_type: lightsim2grid.SolverType
            The new type of solver you want to use. See backend.available_solvers for a list of available solver
            on your machine.
        """
        if not isinstance(solver_type, SolverType):
            raise BackendError(f"The solver type must be from type \"lightsim2grid.SolverType\" and not "
                               f"{type(solver_type)}")
        if solver_type not in self.available_solvers:
            raise BackendError(f"The solver type provided \"{solver_type}\" is not available on your system. Available"
                               f"solvers are {self.available_solvers}")
        self.__current_solver_type = copy.deepcopy(solver_type)
        self._grid.change_solver(self.__current_solver_type)

    def set_solver_max_iter(self, max_iter):
        """
        Set the maximum number of iteration the solver is allowed to perform.

        We do not recommend to modify the default value (10), unless you are using the GaussSeidel powerflow.
        This powerflow being slower, we do not recommend to use it.

        Recommendation:

        - for SolverType.SparseLU: 10
        - for SolverType.GaussSeidel: 10000
        - for SolverType.DC: this has no effect
        - for SolverType.SparseKLU: 10

        Parameters
        ----------
        max_iter: ``int``
            Maximum number of iteration the powerflow can run. It should be number >= 1

        Notes
        -------
        This has to be set for every backend that you want to use. For example, you have to set it
        in the backend of the `_obs_env` of the observation and if you are using "grid2op.MultMixEnv` you
        have to set it in all mixes!

        """
        try:
            max_iter = int(max_iter)
        except Exception as exc_:
            raise BackendError(f"Impossible to convert \"max_iter={max_iter}\" to an integer with exception \"{exc_}\"")
        if max_iter < 1:
            raise BackendError("max_iter should be a strictly positive integer (integer >= 1)")
        self.max_it = max_iter

    def set_tol(self, new_tol):
        """
        Set the tolerance of the powerflow. This means that the powerflow will stop when the Kirchhoff's Circuit Laws
        are met up to a tolerance of "new_tol".

        Decrease the tolerance might speed up the computation of the powerflow but will decrease the accuracy. We do
        not recommend to modify the default value of 1e-8.

        Parameters
        ----------
        new_tol: ``float``
            The new tolerance to use (should be a float > 0)

        Notes
        -------
        This has to be set for every backend that you want to use. For example, you have to set it
        in the backend of the `_obs_env` of the observation and if you are using "grid2op.MultMixEnv` you
        have to set it in all mixes!
        """
        try:
            new_tol = float(new_tol)
        except Exception as exc_:
            raise BackendError(f"Impossible to convert \"new_tol={new_tol}\" to an float with error \"{exc_}\"")
        if new_tol <= 0:
            raise BackendError("new_tol should be a strictly positive float (float > 0)")
        self.tol = new_tol
        self._idx_hack_storage = np.zeros(0, dtype=dt_int)

    def load_grid(self, path=None, filename=None):

        # if self.init_pp_backend is None:
        self.init_pp_backend.load_grid(path, filename)
        self.can_output_theta = True  # i can compute the "theta" and output it to grid2op

        self._grid = init(self.init_pp_backend._grid)
        self.available_solvers = self._grid.available_solvers()
        if SolverType.KLU in self.available_solvers:
            # use the faster KLU if available
            self._grid.change_solver(SolverType.KLU)
        if self.__current_solver_type is None:
            self.__current_solver_type = copy.deepcopy(self._grid.get_solver_type())

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

        # TODO storage check grid2op version and see if storage is available !
        if self.__has_storage:
            self.n_storage = self.init_pp_backend.n_storage
            self.storage_to_subid = self.init_pp_backend.storage_to_subid
            self.storage_pu_to_kv = self.init_pp_backend.storage_pu_to_kv
            self.name_storage = self.init_pp_backend.name_storage
            self.storage_to_sub_pos = self.init_pp_backend.storage_to_sub_pos
            self.storage_type = self.init_pp_backend.storage_type
            self.storage_Emin = self.init_pp_backend.storage_Emin
            self.storage_Emax = self.init_pp_backend.storage_Emax
            self.storage_max_p_prod = self.init_pp_backend.storage_max_p_prod
            self.storage_max_p_absorb = self.init_pp_backend.storage_max_p_absorb
            self.storage_marginal_cost = self.init_pp_backend.storage_marginal_cost
            self.storage_loss = self.init_pp_backend.storage_loss
            self.storage_discharging_efficiency = self.init_pp_backend.storage_discharging_efficiency
            self.storage_charging_efficiency = self.init_pp_backend.storage_charging_efficiency

        self.nb_bus_total = self.init_pp_backend._grid.bus.shape[0]

        self.thermal_limit_a = copy.deepcopy(self.init_pp_backend.thermal_limit_a)

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
        self._init_bus_lor = np.concatenate((self._init_bus_lor, t_for)).astype(int)
        self._init_bus_lex = np.concatenate((self._init_bus_lex, t_fex)).astype(int)
        self._init_bus_load = self._init_bus_load.astype(int)
        self._init_bus_gen = self._init_bus_gen.astype(int)

        tmp = self._init_bus_lor + self.__nb_bus_before
        self._init_bus_lor = np.concatenate((self._init_bus_lor.reshape(-1, 1),
                                             tmp.reshape(-1, 1)), axis=-1)
        tmp = self._init_bus_lex + self.__nb_bus_before
        self._init_bus_lex = np.concatenate((self._init_bus_lex.reshape(-1, 1),
                                             tmp.reshape(-1, 1)), axis=-1)
        tmp = self._init_bus_load + self.__nb_bus_before
        self._init_bus_load = np.concatenate((self._init_bus_load.reshape(-1, 1),
                                             tmp.reshape(-1, 1)), axis=-1)
        tmp = self._init_bus_gen + self.__nb_bus_before
        self._init_bus_gen = np.concatenate((self._init_bus_gen.reshape(-1, 1),
                                             tmp.reshape(-1, 1)), axis=-1)

        self._big_topo_to_obj = [(None, None) for _ in range(self.dim_topo)]

        self._compute_pos_big_topo()

        # set up the "lightsim grid" accordingly
        self._grid.set_n_sub(self.__nb_bus_before)
        self._grid.set_load_pos_topo_vect(self.load_pos_topo_vect)
        self._grid.set_gen_pos_topo_vect(self.gen_pos_topo_vect)
        self._grid.set_line_or_pos_topo_vect(self.line_or_pos_topo_vect[:self.__nb_powerline])
        self._grid.set_line_ex_pos_topo_vect(self.line_ex_pos_topo_vect[:self.__nb_powerline])
        self._grid.set_trafo_hv_pos_topo_vect(self.line_or_pos_topo_vect[self.__nb_powerline:])
        self._grid.set_trafo_lv_pos_topo_vect(self.line_ex_pos_topo_vect[self.__nb_powerline:])

        self._grid.set_load_to_subid(self.load_to_subid)
        self._grid.set_gen_to_subid(self.gen_to_subid)
        self._grid.set_line_or_to_subid(self.line_or_to_subid[:self.__nb_powerline])
        self._grid.set_line_ex_to_subid(self.line_ex_to_subid[:self.__nb_powerline])
        self._grid.set_trafo_hv_to_subid(self.line_or_to_subid[self.__nb_powerline:])
        self._grid.set_trafo_lv_to_subid(self.line_ex_to_subid[self.__nb_powerline:])

        # TODO storage check grid2op version and see if storage is available !
        if self.__has_storage:
            self._grid.set_storage_to_subid(self.storage_to_subid)
            self._grid.set_storage_pos_topo_vect(self.storage_pos_topo_vect)

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

        # TODO storage check grid2op version and see if storage is available !
        if self.__has_storage:
            nm_ = "storage"
            for l_id, pos_big_topo  in enumerate(self.storage_pos_topo_vect):
                self._big_topo_to_obj[pos_big_topo] = (l_id, nm_)

        self.prod_p = 1.0 * self.init_pp_backend._grid.gen["p_mw"].values
        self.next_prod_p = 1.0 * self.init_pp_backend._grid.gen["p_mw"].values

        # for shunts
        self.n_shunt = self.init_pp_backend.n_shunt
        self.shunt_to_subid = self.init_pp_backend.shunt_to_subid
        self.name_shunt = self.init_pp_backend.name_shunt

        if hasattr(self.init_pp_backend, "_sh_vnkv"):
            # attribute has been added in grid2op ~1.3 or 1.4
            self._sh_vnkv = self.init_pp_backend._sh_vnkv

        self.shunts_data_available = self.init_pp_backend.shunts_data_available

        # number of object per bus, to activate, deactivate them
        self.nb_obj_per_bus = np.zeros(2 * self.__nb_bus_before, dtype=dt_int)

        self.topo_vect = np.ones(self.dim_topo, dtype=dt_int)
        if self.shunts_data_available:
            self.shunt_topo_vect = np.ones(self.n_shunt, dtype=dt_int)

        self.p_or = np.full(self.n_line, dtype=dt_float, fill_value=np.NaN)
        self.q_or = np.full(self.n_line, dtype=dt_float, fill_value=np.NaN)
        self.v_or = np.full(self.n_line, dtype=dt_float, fill_value=np.NaN)
        self.a_or = np.full(self.n_line, dtype=dt_float, fill_value=np.NaN)
        self.p_ex = np.full(self.n_line, dtype=dt_float, fill_value=np.NaN)
        self.q_ex = np.full(self.n_line, dtype=dt_float, fill_value=np.NaN)
        self.v_ex = np.full(self.n_line, dtype=dt_float, fill_value=np.NaN)
        self.a_ex = np.full(self.n_line, dtype=dt_float, fill_value=np.NaN)

        self.load_p = np.full(self.n_load, dtype=dt_float, fill_value=np.NaN)
        self.load_q = np.full(self.n_load, dtype=dt_float, fill_value=np.NaN)
        self.load_v = np.full(self.n_load, dtype=dt_float, fill_value=np.NaN)

        self.prod_p = np.full(self.n_gen, dtype=dt_float, fill_value=np.NaN)
        self.prod_q = np.full(self.n_gen, dtype=dt_float, fill_value=np.NaN)
        self.prod_v = np.full(self.n_gen, dtype=dt_float, fill_value=np.NaN)

        # TODO storage check grid2op version and see if storage is available !
        if self.__has_storage:
            self.storage_p = np.full(self.n_storage, dtype=dt_float, fill_value=np.NaN)
            self.storage_q = np.full(self.n_storage, dtype=dt_float, fill_value=np.NaN)
            self.storage_v = np.full(self.n_storage, dtype=dt_float, fill_value=np.NaN)

        self._count_object_per_bus()
        self.__me_at_init = self._grid.copy()
        self.__init_topo_vect = np.ones(self.dim_topo, dtype=dt_int)
        self.__init_topo_vect[:] = self.topo_vect

    def assert_grid_correct_after_powerflow(self):
        """
        This method is called by the environment. It ensure that the backend remains consistent even after a powerflow
        has be run with :func:`Backend.runpf` method.

        :return: ``None``
        :raise: :class:`grid2op.Exceptions.EnvError` and possibly all of its derived class.
        """
        # test the results gives the proper size
        super().assert_grid_correct_after_powerflow()
        self.init_pp_backend.__class__ = self.init_pp_backend.init_grid(self)
        self._backend_action_class = _BackendAction.init_grid(self)
        self._init_action_to_set = self._backend_action_class()
        try:
            # feature added in grid2op 1.4 or 1.5
            _init_action_to_set = self.get_action_to_set()
        except TypeError:
            _init_action_to_set = self._get_action_to_set_deprecated()
        self._init_action_to_set += _init_action_to_set

    def _get_action_to_set_deprecated(self):
        warnings.warn("DEPRECATION: grid2op <=1.4 is not well supported with lightsim2grid. Lots of bugs have been"
                      "fixed since then. Please upgrade to grid2op >= 1.5",
                      DeprecationWarning)
        line_status = self.get_line_status()
        line_status = 2 * line_status - 1
        line_status = line_status.astype(dt_int)
        topo_vect = self.get_topo_vect()
        prod_p, _, prod_v = self.generators_info()
        load_p, load_q, _ = self.loads_info()
        complete_action_class = CompleteAction.init_grid(self)
        set_me = complete_action_class()
        set_me.update({"set_line_status": line_status,
                       "set_bus": topo_vect})
        return set_me

    def _count_object_per_bus(self):
        # should be called only when self.topo_vect and self.shunt_topo_vect are set
        # todo factor that more properly to update it when it's modified, and not each time

        self.nb_obj_per_bus = np.zeros(2 * self.__nb_bus_before, dtype=dt_int)

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

    def close(self):
        self.init_pp_backend.close()
        self._grid = None

    def apply_action(self, backendAction):
        """
        Specific implementation of the method to apply an action modifying a powergrid in the pandapower format.
        """
        active_bus, *_, topo__, shunts__ = backendAction()

        # update the injections
        self._grid.update_gens_p(backendAction.prod_p.changed,
                                 backendAction.prod_p.values)
        self._grid.update_gens_v(backendAction.prod_v.changed,
                                 backendAction.prod_v.values / self.prod_pu_to_kv)
        self._grid.update_loads_p(backendAction.load_p.changed,
                                  backendAction.load_p.values)
        self._grid.update_loads_q(backendAction.load_q.changed,
                                  backendAction.load_q.values)
        if self.__has_storage:
            # TODO
            # reactivate the storage that i deactivate because of the "hack". See
            # for stor_id in self._idx_hack_storage:
            #     self._grid.reactivate_storage(stor_id)
            self._grid.update_storages_p(backendAction.storage_power.changed,
                                         backendAction.storage_power.values)

        # handle shunts
        if self.shunts_data_available:
            shunt_p, shunt_q, shunt_bus = backendAction.shunt_p, backendAction.shunt_q, backendAction.shunt_bus
            # shunt topology
            # (need to be done before to avoid error like "impossible to set reactive value of a disconnected shunt")
            for sh_id, new_bus in shunt_bus:
                if new_bus == -1:
                    self._grid.deactivate_shunt(sh_id)
                else:
                    self._grid.reactivate_shunt(sh_id)
                    self._grid.change_bus_shunt(sh_id, self.shunt_to_subid[sh_id] * new_bus)

            for sh_id, new_p in shunt_p:
                self._grid.change_p_shunt(sh_id, new_p)
            for sh_id, new_q in shunt_q:
                self._grid.change_q_shunt(sh_id, new_q)

        # and now change the overall topology
        # TODO hack for storage units: if 0. production i pretend they are disconnected on the
        # TODO c++ side
        # this is to deal with the test that "if a storage unit is alone on a bus, but produces 0, then it's fine)
        # if self.__has_storage and self.n_storage > 0:
        #     chgt = copy.copy(backendAction.current_topo.changed)
        #     my_val = 1 * backendAction.current_topo.values
        #     self._idx_hack_storage = np.where((backendAction.storage_power.values == 0.))[0]
        #     idx_storage_topo = self.storage_pos_topo_vect[self._idx_hack_storage]
        #     changed[idx_storage_topo] = my_val[idx_storage_topo] != -1
        #     my_val[idx_storage_topo] = -1
        # else:
        #     self._idx_hack_storage = []
        #     chgt = backendAction.current_topo.changed
        #     my_val = backendAction.current_topo.values
        # self._grid.update_topo(changed, my_val)
        chgt = backendAction.current_topo.changed
        self._grid.update_topo(chgt, backendAction.current_topo.values)
        self.topo_vect[chgt] = backendAction.current_topo.values[chgt]
        # TODO c++ side: have a check to be sure that the set_***_pos_topo_vect and set_***_to_sub_id
        # TODO have been correctly called before calling the function self._grid.update_topo

    def runpf(self, is_dc=False):
        my_exc_ = None
        res = False
        try:
            if is_dc:
                msg_ = "LightSimBackend: the support of the DC approximation is fully supported at the moment"
                warnings.warn(msg_)
                raise RuntimeError(msg_)
                if self.V is None:
                    self.V = np.ones(self.nb_bus_total, dtype=np.complex_) * self._grid.get_init_vm_pu()
                V = self._grid.dc_pf(self.V, self.max_it, self.tol)
                if V.shape[0] == 0:
                    raise DivergingPowerFlow("divergence of powerflow (non connected grid)")
            else:
                if (self.V is None) or (self.V.shape[0] == 0):
                    # init from dc approx in this case
                    self.V = np.ones(self.nb_bus_total, dtype=np.complex_) * self._grid.get_init_vm_pu()

                if self.initdc:
                    self._grid.deactivate_result_computation()
                    V = self._grid.dc_pf(copy.deepcopy(self.V), self.max_it, self.tol)
                    self._grid.reactivate_result_computation()
                    if V.shape[0] == 0:
                        raise DivergingPowerFlow("divergence of powerflow (non connected grid)")
                    self.V[:] = V

                V = self._grid.ac_pf(self.V, self.max_it, self.tol)
                if V.shape[0] == 0:
                    # V = self._grid.ac_pf(self.V, self.max_it, self.tol)
                    raise DivergingPowerFlow("divergence of powerflow")

            self.comp_time += self._grid.get_computation_time()
            self.V[:] = V
            lpor, lqor, lvor, laor = self._grid.get_lineor_res()
            lpex, lqex, lvex, laex = self._grid.get_lineex_res()
            tpor, tqor, tvor, taor = self._grid.get_trafohv_res()
            tpex, tqex, tvex, taex = self._grid.get_trafolv_res()

            self.p_or[:] = np.concatenate((lpor, tpor))
            self.q_or[:] = np.concatenate((lqor, tqor))
            self.v_or[:] = np.concatenate((lvor, tvor))
            self.a_or[:] = 1000. * np.concatenate((laor, taor))

            self.p_ex[:] = np.concatenate((lpex, tpex))
            self.q_ex[:] = np.concatenate((lqex, tqex))
            self.v_ex[:] = np.concatenate((lvex, tvex))
            self.a_ex[:] = 1000. * np.concatenate((laex, taex))

            self.a_or[~np.isfinite(self.a_or)] = 0.
            self.v_or[~np.isfinite(self.v_or)] = 0.
            self.a_ex[~np.isfinite(self.a_ex)] = 0.
            self.v_ex[~np.isfinite(self.v_ex)] = 0.

            self.load_p[:], self.load_q[:], self.load_v[:] = self._grid.get_loads_res()
            self.prod_p[:], self.prod_q[:], self.prod_v[:] = self._grid.get_gen_res()
            if self.__has_storage:
                self.storage_p[:], self.storage_q[:], self.storage_v[:] = self._grid.get_storages_res()
            self.next_prod_p[:] = self.prod_p

            if np.any(~np.isfinite(self.load_v)) or np.any(self.load_v <= 0.):
                raise DivergingPowerFlow("One load is disconnected")
            if np.any(~np.isfinite(self.prod_v)) or np.any(self.prod_v <= 0.):
                raise DivergingPowerFlow("One generator is disconnected")
            # TODO storage case of divergence !

            res = True
        except Exception as exc_:
            # of the powerflow has not converged, results are Nan
            self._fill_nans()
            res = False
            my_exc_ = exc_

        # TODO grid2op compatibility ! (was a single returned element before storage were introduced)
        if self.__has_storage:
            res = res, my_exc_
        return res

    def _fill_nans(self):
        """fill the results vectors with nans"""
        self.p_or[:] = np.NaN
        self.q_or[:] = np.NaN
        self.v_or[:] = np.NaN
        self.a_or[:] = np.NaN
        self.p_ex[:] = np.NaN
        self.q_ex[:] = np.NaN
        self.v_ex[:] = np.NaN
        self.a_ex[:] = np.NaN
        self.load_p[:] = np.NaN
        self.load_q[:] = np.NaN
        self.load_v[:] = np.NaN
        self.prod_p[:] = np.NaN
        self.next_prod_p[:] = np.NaN
        self.prod_q[:] = np.NaN
        self.prod_v[:] = np.NaN
        self.topo_vect[:] = -1
        if self.__has_storage:
            self.storage_p[:] = np.NaN
            self.storage_q[:] = np.NaN
            self.storage_v[:] = np.NaN
        res = False

    def copy(self):
        # i can perform a regular copy, everything has been initialized
        mygrid = self._grid
        __me_at_init = self.__me_at_init
        inippbackend = self.init_pp_backend
        if __me_at_init is None:
            # __me_at_init is defined as being the copy of the grid,
            # if it's not defined then i can define it here.
            __me_at_init = self._grid.copy()

        self._grid = None
        self.__me_at_init = None
        self.init_pp_backend = None
        res = copy.deepcopy(self)

        res._grid = mygrid.copy()
        res.__me_at_init = __me_at_init.copy()
        res.init_pp_backend = inippbackend.copy()

        self._grid = mygrid
        self.init_pp_backend = inippbackend
        self.__me_at_init = __me_at_init
        return res

    def get_line_status(self):
        l_s = self._grid.get_lines_status()
        t_s = self._grid.get_trafo_status()
        return np.concatenate((l_s, t_s)).astype(dt_bool)

    def get_line_flow(self):
        return self.a_or

    def _grid2op_bus_from_klu_bus(self, klu_bus):
        res = 0
        if klu_bus != 0:
            # object is connected
            res = 1 if klu_bus < self.__nb_bus_before else 2
        return res

    def _klu_bus_from_grid2op_bus(self, grid2op_bus, grid2op_bus_init):
        return grid2op_bus_init[grid2op_bus - 1]

    def get_topo_vect(self):
        return self.topo_vect

    def generators_info(self):
        return self.cst_1 * self.prod_p, self.cst_1 * self.prod_q, self.cst_1 * self.prod_v

    def loads_info(self):
        return self.cst_1 * self.load_p, self.cst_1 * self.load_q, self.cst_1 * self.load_v

    def lines_or_info(self):
        return self.cst_1 * self.p_or, self.cst_1 * self.q_or, self.cst_1 * self.v_or, self.cst_1 * self.a_or

    def lines_ex_info(self):
        return self.cst_1 * self.p_ex, self.cst_1 * self.q_ex, self.cst_1 * self.v_ex, self.cst_1 * self.a_ex

    def storages_info(self):
        if not self.__has_storage:
            raise RuntimeError("Storage units are not supported with your grid2op version. Please upgrade to "
                               "grid2op >1.5")
        return self.cst_1 * self.storage_p, self.cst_1 * self.storage_q, self.cst_1 * self.storage_v

    def shunt_info(self):
        tmp = self._grid.get_shunts_res()
        shunt_bus = np.array([self._grid.get_bus_shunt(i) for i in range(self.n_shunt)], dtype=dt_int)
        res_bus = np.ones(shunt_bus.shape[0], dtype=dt_int)  # by default all shunts are on bus one
        res_bus[shunt_bus >= self.__nb_bus_before] = 2  # except the one that are connected to bus 2
        res_bus[shunt_bus == -1] = -1  # or the one that are disconnected
        return tmp[0].astype(dt_float), tmp[1].astype(dt_float), tmp[2].astype(dt_float), res_bus

    def _disconnect_line(self, id_):
        self.topo_vect[self.line_ex_pos_topo_vect[id_]] = -1
        self.topo_vect[self.line_or_pos_topo_vect[id_]] = -1
        if id_ < self.__nb_powerline:
            self._grid.deactivate_powerline(id_)
        else:
            self._grid.deactivate_trafo(id_ - self.__nb_powerline)

    def get_current_solver_type(self):
        return self.__current_solver_type

    def reset(self, grid_path, grid_filename=None):
        self.V = None
        self._fill_nans()
        self._grid = self.__me_at_init.copy()
        self._grid.change_solver(self.__current_solver_type)
        self.topo_vect[:] = self.__init_topo_vect
        self.comp_time = 0.
