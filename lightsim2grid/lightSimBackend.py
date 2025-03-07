# Copyright (c) 2020-2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import copy
from typing import Tuple, Optional, Any, Dict, Union
from packaging import version
try:
    from typing import Self
except ImportError:
    # python version is probably bellow 3.11
    from typing_extensions import Self
    
import warnings
import numpy as np
import pandas as pd
import time

import grid2op
from grid2op.Action import CompleteAction
from grid2op.Backend import Backend
from grid2op.Exceptions import BackendError, Grid2OpException
from grid2op.dtypes import dt_float, dt_int, dt_bool
try:
    from grid2op.Action._backendAction import _BackendAction
except ImportError as exc_:
    from grid2op.Action._BackendAction import _BackendAction

try:
    from grid2op.Space import DEFAULT_N_BUSBAR_PER_SUB
except ImportError:
    # for backward compatibility with grid2op <= 1.9.8
    DEFAULT_N_BUSBAR_PER_SUB = 2

try:
    from grid2op.Space import DEFAULT_ALLOW_DETACHMENT
except ImportError:
    # for backward compatibility with grid2op < 1.11.0
    DEFAULT_ALLOW_DETACHMENT = False
    
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal
    
from lightsim2grid.solver import SolverType


grid2op_min_cls_attr_ver = version.parse("1.6.4")
grid2op_min_shunt_cls_properly_handled = version.parse("1.9.7")


class LightSimBackend(Backend):
    """
    This is a specialization of the grid2op Backend class to use the lightsim2grid solver,
    coded in c++, aiming at speeding up the computations.
    """
    
    shunts_data_available = True
    glop_version = Backend.glop_version if hasattr(Backend, "glop_version") else grid2op.__version__
    if not hasattr(Backend, "n_busbar_per_sub"):
        # for legacy grid2op
        n_busbar_per_sub = DEFAULT_N_BUSBAR_PER_SUB
        
    if not hasattr(Backend, "detachment_is_allowed"):
        # for legacy grid2op (< 1.11.0)
        detachment_is_allowed = DEFAULT_ALLOW_DETACHMENT
        
    def __init__(self,
                 detailed_infos_for_cascading_failures: bool=False,
                 can_be_copied: bool=True,
                 max_iter: int=10,
                 tol: float=1e-8,
                 solver_type: Optional[SolverType]=None,
                 turned_off_pv : bool=True,  # are gen turned off (or with p=0) contributing to voltage or not
                 dist_slack_non_renew: bool=False,  # distribute the slack on non renewable turned on (and with P>0) generators
                 use_static_gen: bool=False, # add the static generators as generator gri2dop side
                 loader_method: Literal["pandapower", "pypowsybl"] = "pandapower",
                 loader_kwargs : Optional[dict] = None,
                 stop_if_load_disco : Optional[bool] = None,
                 stop_if_gen_disco : Optional[bool] = None,
                 stop_if_storage_disco : Optional[bool] = None,
                 ):
        #: ``int`` maximum number of iteration allowed for the solver
        #: if the solver has not converge after this, it will 
        #: send a "divergence error"
        self.max_it = max_iter
        
        #: ``float`` tolerance of the solver
        self.tol = tol  # tolerance for the solver
        
        self._check_suitable_solver_type(solver_type, check_in_avail_solver=False)
        self.__current_solver_type = solver_type
        
        #: does the "turned off" generators (including when p=0)
        #: are pv buses
        self._turned_off_pv = turned_off_pv
        
        #: distributed slack, on non renewable gen with P > 0
        self._dist_slack_non_renew = dist_slack_non_renew
        
        # add the static gen to the list of controlable gen in grid2Op
        self._use_static_gen = use_static_gen  # TODO implement it

        #: For now, you can initialize a "lightsim2grid" LightsimBackend
        #: either from pypowsybl or from pandapower.
        #: Use with `LightsimBackend(..., loader_method="pypowsybl")`
        #: or `LightsimBackend(..., loader_method="pandapower")` (default)
        self._loader_method = loader_method
        
        #: Which key-word arguments will be used to initialize the Gridmodel
        #: either from pandapower or lightsim2grid.
        #: It is not currently used for pandapower.
        #:
        #: For pypowsybl it can contain the following keys:
        #: 
        #:   - `n_busbar_per_sub` (``int``): number of independant buses for
        #:     each substation in the GridModel.
        #:   - `use_buses_for_sub` (``bool``): whether to use the buses (in the 
        #:     pypowsybl Network) to initialize the "substation" in the lightsim2grid
        #:     Gridmodel (if ``True``). If ``False`` it will use the `voltage_levels`
        #:     of the pypowsybl Network.
        #:   - `gen_slack_id` (``int``): which generator will be used for the slack
        #:   - `use_grid2op_default_names` (``bool``): whether to use the "default names"
        #:     assigne by grid2op or to read them from the the iidm grid.
        #:   - `reconnect_disco_gen` (``bool``): whether to reconnect the disconnected 
        #:     generators in the iidm grid. If set to ``True`` then the generators will be
        #:     reconnected with a p setpoint of 0. MW and a voltage setpoint of 1 pu
        #:   - `reconnect_disco_load` (``bool``): whether to reconnec the disconnected
        #:     load from in the iidm grid. If reconnected, load will have a target
        #:     p set to 0. and a target q set to 0.
        #: 
        self._loader_kwargs = loader_kwargs

        #: .. versionadded:: 0.8.0
        #:
        #: if set to ``True`` (default) then the backend will raise a 
        #: BackendError in case of disconnected load
        #:
        #: .. note:: When using grid2op >= 1.11.0 and lightsim2grid >= 0.10.0 this 
        #:    is automatically set-up with the call to grid2op.make, 
        #:    to match the behaviour expected by `allow_detachment`
        self._stop_if_load_disco = stop_if_load_disco
        
        #: .. versionadded:: 0.8.0
        #:
        #: if set to ``True`` (default) then the backend will raise a 
        #: BackendError in case of disconnected generator
        #:
        #: .. note:: When using grid2op >= 1.11.0 and lightsim2grid >= 0.10.0 this 
        #:    is automatically set-up with the call to grid2op.make, 
        #:    to match the behaviour expected by `allow_detachment`
        self._stop_if_gen_disco = stop_if_gen_disco
        
        #: .. versionadded:: 0.10.0
        #:
        #: if set to ``True`` (default) then the backend will raise a 
        #: BackendError in case of disconnected storage that are 
        #: asked to produce / absorb something
        #:
        #: .. note:: When using grid2op >= 1.11.0 and lightsim2grid >= 0.10.0 this 
        #:    is automatically set-up with the call to grid2op.make, 
        #:    to match the behaviour expected by `allow_detachment`
        self._stop_if_storage_disco = stop_if_storage_disco
                                        
        self._aux_init_super(detailed_infos_for_cascading_failures,
                             can_be_copied,
                             solver_type,
                             max_iter,
                             tol,
                             turned_off_pv,
                             dist_slack_non_renew,
                             use_static_gen,
                             loader_method,
                             loader_kwargs,
                             stop_if_load_disco,
                             stop_if_gen_disco,
                             stop_if_storage_disco)
        
        # backward compat: need to define it if not done by grid2op
        if not hasattr(self, "_can_be_copied"):
            self._can_be_copied = can_be_copied

        #: I can output the voltage angle
        self.can_output_theta = True
        
        #: .. versionadded:: 0.8.0
        #:
        #: Which type of grid format can be read by your backend.
        #: It is "json" if loaded from pandapower or
        #: "xiidm" if loaded from pypowsybl.
        self.supported_grid_format = None
        
        if loader_method == "pandapower":
            self.supported_grid_format = ("json", )  # new in 0.8.0
        elif loader_method == "pypowsybl":
            self.supported_grid_format = ("xiidm", )  # new in 0.8.0
        else:
            raise BackendError(f"Uknown loader_method : '{loader_method}'")
        
        # lazy loading because it crashes...
        from lightsim2grid._utils import _DoNotUseAnywherePandaPowerBackend
        from grid2op.Space import GridObjects  # lazy import
        self.__has_storage = hasattr(GridObjects, "n_storage")
        if not self.__has_storage:
            warnings.warn("Please upgrade your grid2op to >= 1.5.0. You are using a backward compatibility "
                          "feature that will be removed in further lightsim2grid version.")
        if version.parse(grid2op.__version__) < grid2op_min_shunt_cls_properly_handled:
            warnings.warn(f"You are using a legacy grid2op version. It is not possible to deactivate the shunts in lightsim2grid. "
                          f"Please upgrade to grid2op >= {grid2op_min_shunt_cls_properly_handled}")
            self.shunts_data_available = True  # needs to be self and not type(self) here
            type(self).shunts_data_available = True
            
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
        self._timer_preproc = 0.
        self._timer_postproc = 0.
        self._timer_solver = 0.
        self._timer_read_data_back = 0.
        self._timer_fetch_data_cpp = 0.
        self._timer_apply_act = 0.

        self.next_prod_p = None  # this vector is updated with the action that will modify the environment
        # it is done to keep track of the redispatching

        self.topo_vect = None
        self.shunt_topo_vect = None
        try:
            self.init_pp_backend = _DoNotUseAnywherePandaPowerBackend(with_numba=False)
        except TypeError as exc_:
            # oldest version of grid2op do not support the kwargs "with_numba"
            # (before 1.9.1)
            self.init_pp_backend = _DoNotUseAnywherePandaPowerBackend()
        
        self.V = None

        self.prod_pu_to_kv = None
        self.load_pu_to_kv = None
        self.lines_or_pu_to_kv = None
        self.lines_ex_pu_to_kv = None
        self.storage_pu_to_kv = None

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

        # shunts
        self.sh_p = None
        self.sh_q = None
        self.sh_v = None
        self.sh_theta = None
        self.sh_bus = None

        # voltage angle
        self.line_or_theta = None
        self.line_ex_theta = None
        self.load_theta = None
        self.gen_theta = None
        self.storage_theta = None

        self.thermal_limit_a = None

        self.dim_topo = -1
        self._init_action_to_set = None
        self._backend_action_class = None
        self.cst_1 = dt_float(1.0)
        self.__me_at_init = None
        self.__init_topo_vect = None

        # available solver in lightsim
        self.available_solvers = []
        
        # computation time of just the powerflow (when the grid is formatted 
        # by the gridmodel already)
        # it takes only into account the time spend in the powerflow algorithm
        self.comp_time = 0.  

        #: computation time of the powerflow
        #: it takes into account everything in the gridmodel, including the mapping 
        #: to the solver, building of Ybus and Sbus AND the time to solve the powerflow
        self.timer_gridmodel_xx_pf = 0.
        
        self._timer_postproc = 0.
        self._timer_preproc = 0.
        self._timer_solver = 0.
        self._timer_read_data_back = 0.
        self._timer_fetch_data_cpp = 0.
        self._timer_apply_act = 0.

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
        
        #: ..versionadded: 0.9.3
        #: sometimes some actions will make the grid fails
        #: but grid2op expect it to fail not on "apply_action"
        #: but rather when calling `runpf`
        #: this flags remembers it
        self._next_pf_fails : Optional[BackendError] = None
        
        # speed optimization
        self._lineor_res = None
        self._lineex_res = None
        self._load_res = None
        self._gen_res = None
        self._shunt_res = None
        self._trafo_hv_res = None
        self._trafo_lv_res = None
        self._storage_res = None
        self._reset_res_pointers()
        self._debug_Vdc = None   # use only for debug !
        self._orig_grid_pypowsybl = None
        
    def _aux_init_super(self, 
                        detailed_infos_for_cascading_failures,
                        can_be_copied,
                        solver_type,
                        max_iter,
                        tol,
                        turned_off_pv,
                        dist_slack_non_renew,
                        use_static_gen,
                        loader_method,
                        loader_kwargs,
                        stop_if_load_disco,
                        stop_if_gen_disco,
                        stop_if_storage_disco):
        try:
            # for grid2Op >= 1.7.1
            Backend.__init__(self,
                             detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures,
                             can_be_copied=can_be_copied,
                             solver_type=solver_type,
                             max_iter=max_iter,
                             tol=tol,
                             turned_off_pv=turned_off_pv,
                             dist_slack_non_renew=dist_slack_non_renew,
                             use_static_gen=use_static_gen,
                             loader_method=loader_method,
                             loader_kwargs=loader_kwargs,
                             stop_if_load_disco=stop_if_load_disco,
                             stop_if_gen_disco=stop_if_gen_disco,
                             stop_if_storage_disco=stop_if_storage_disco
                             )
        except TypeError as exc_:
            warnings.warn("Please use grid2op >= 1.7.1: with older grid2op versions, "
                          "you cannot set max_iter, tol nor solver_type arguments.")
            Backend.__init__(self,
                             detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
            
        if hasattr(type(self), "can_handle_more_than_2_busbar"):
            # do not forget to propagate this if needed
            self.can_handle_more_than_2_busbar()
            
        if hasattr(type(self), "can_handle_detachment"):
            # do not forget to propagate this if needed
            self.can_handle_detachment()
            
    def turnedoff_no_pv(self):
        self._turned_off_pv = False
        self._grid.turnedoff_no_pv()
        
    def turnedoff_pv(self):
        self._turned_off_pv = True
        self._grid.turnedoff_pv()

    def get_theta(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
        return self.cst_1 * self.line_or_theta, \
               self.cst_1 * self.line_ex_theta, \
               self.cst_1 * self.load_theta, \
               self.cst_1 * self.gen_theta, \
               self.cst_1 * self.storage_theta

    def get_solver_types(self) -> Union[SolverType, SolverType]:
        """Return the types of solver that are used in the form a tuple with 2 elements.
        
        The first one is the solver used for AC computation, the second one for DC computation (and also for
        initialization of AC computations)
        
        You can use it to check what is used:
        
        .. code-block:: python

            import grid2op
            import lightsim2grid
            from ligthsim2grid import LightSimBackend
            
            env_name = ...
            env = grid2op.make(env_name, backend=LightSimBackend())
            print(env.backend.get_solver_types())
            # >>> (<SolverType.KLUSingleSlack: 7>, <SolverType.KLUDC: 9>)  [can depend on your installation of lightsim2grid]
            
            env2 = grid2op.make(env_name, backend=LightSimBackend(solver_type=lightsim2grid.solver.SolverType.SparseLU))
            print(env2.backend.get_solver_types())
            # >>> (<SolverType.SparseLU: 0>, <SolverType.KLUDC: 9>)  [can depend on your installation of lightsim2grid]
            
        """
        return self._grid.get_solver_type(), self._grid.get_dc_solver_type()
        
    def set_solver_type(self, solver_type: SolverType) -> None:
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
        if solver_type is None:
            raise BackendError("Impossible to change the solver type to None. Please enter a valid solver type.")
        self._check_suitable_solver_type(solver_type)
        self.__current_solver_type = copy.deepcopy(solver_type)
        self._grid.change_solver(self.__current_solver_type)

    def _check_suitable_solver_type(self, solver_type, check_in_avail_solver=True):
        if solver_type is None:
            return
        
        if not isinstance(solver_type, SolverType):
            raise BackendError(f"The solver type must be from type \"lightsim2grid.SolverType\" and not "
                               f"{type(solver_type)}")
            
        if check_in_avail_solver and solver_type not in self.available_solvers:
            raise BackendError(f"The solver type provided \"{solver_type}\" is not available on your system. Available"
                               f"solvers are {self.available_solvers}")
            
    def set_solver_max_iter(self, max_iter: int) -> None:
        """
        Set the maximum number of iteration the solver is allowed to perform.

        We do not recommend to modify the default value (10), unless you are using the GaussSeidel powerflow.
        This powerflow being slower, we do not recommend to use it.

        Recommendation, for medium sized grid (**eg** based on the ieee 118):

        - for SolverType.SparseLU: 10
        - for SolverType.GaussSeidel: 10000
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

    def set_tol(self, new_tol: float) -> None:
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
            raise BackendError(f"Impossible to convert \"new_tol={new_tol}\" to an float with error \"{exc_}\"") from exc_
        if new_tol <= 0:
            raise BackendError("new_tol should be a strictly positive float (float > 0)")
        self.tol = new_tol
        self._idx_hack_storage = np.zeros(0, dtype=dt_int)

    def _handle_turnedoff_pv(self):
        if self._turned_off_pv:
            self._grid.turnedoff_pv()
        else:
            self._grid.turnedoff_no_pv()
    
    def _assign_right_solver(self):
        slack_weights = np.array([el.slack_weight for el in self._grid.get_generators()])
        nb_slack_nonzero = (np.abs(slack_weights) > 0.).sum()
        has_single_slack = nb_slack_nonzero == 1
        if has_single_slack and not self._dist_slack_non_renew:
            if SolverType.KLUSingleSlack in self.available_solvers:
                # use the faster KLU if available
                self._grid.change_solver(SolverType.KLUSingleSlack)
            else:
                self._grid.change_solver(SolverType.SparseLUSingleSlack)
        else:
            # grid has multiple slack      
            if SolverType.KLU in self.available_solvers:
                # use the faster KLU if available
                self._grid.change_solver(SolverType.KLU)
            else:
                self._grid.change_solver(SolverType.SparseLU)
    
    def _aux_set_correct_detach_flags_d_allowed(self):
        # user allowed detachment, I check the correct flags
        if self._stop_if_gen_disco is None:
            # user did not specify anything
            self._stop_if_gen_disco = False
        elif not self._stop_if_gen_disco:
            # force conversion to proper type
            self._stop_if_gen_disco = False
        elif self._stop_if_gen_disco:
            # erase default values and continue like the grid2op call specifies
            warnings.warn("Call to `grid2op.make(..., allow_detachement=True)` will erase the lightsim2grid kwargs `stop_if_gen_disco=True`")
            self._stop_if_gen_disco = False
            
        if self._stop_if_load_disco is None:
            # user did not specify anything
            self._stop_if_load_disco = False
        elif not self._stop_if_load_disco:
            # force conversion to the proper type
            self._stop_if_load_disco = False
        elif self._stop_if_load_disco:
            # erase default values and continue like the grid2op call specifies
            warnings.warn("Call to `grid2op.make(..., allow_detachement=True)` will erase the lightsim2grid kwargs `stop_if_load_disco=True`")
            self._stop_if_load_disco = False
            
        if self._stop_if_storage_disco is None:
            # user did not specify anything
            self._stop_if_storage_disco = False
        elif not self._stop_if_storage_disco:
            # force conversion to the proper type
            self._stop_if_storage_disco = False
        elif self._stop_if_storage_disco:
            # erase default values and continue like the grid2op call specifies
            warnings.warn("Call to `grid2op.make(..., allow_detachement=True)` will erase the lightsim2grid kwargs `stop_if_storage_disco=True`")
            self._stop_if_storage_disco = False
    
    def _aux_set_correct_detach_flags_d_not_allowed(self):# user did not allow detachment (or it's a legacy grid2op version), I check the correct flags
        if self._stop_if_gen_disco is None:
            # user did not specify anything
            self._stop_if_gen_disco = True
        elif self._stop_if_gen_disco:
            # force conversion to proper type
            self._stop_if_gen_disco = True
        elif not self._stop_if_gen_disco:
            # erase default values and continue like the grid2op call specifies
            warnings.warn("Call to `grid2op.make(..., allow_detachement=False)` will erase the lightsim2grid kwargs `stop_if_gen_disco=False`")
            self._stop_if_gen_disco = True
            
        if self._stop_if_load_disco is None:
            # user did not specify anything
            self._stop_if_load_disco = True
        elif self._stop_if_load_disco:
            # force conversion to proper type
            self._stop_if_load_disco = True
        elif not self._stop_if_load_disco:
            # erase default values and continue like the grid2op call specifies
            warnings.warn("Call to `grid2op.make(..., allow_detachement=False)` will erase the lightsim2grid kwargs `stop_if_load_disco=False`")
            self._stop_if_load_disco = True
            
        if self._stop_if_storage_disco is None:
            # user did not specify anything
            self._stop_if_storage_disco = True
        elif self._stop_if_storage_disco:
            # force conversion to proper type
            self._stop_if_storage_disco = True
        elif not self._stop_if_storage_disco:
            # erase default values and continue like the grid2op call specifies
            warnings.warn("Call to `grid2op.make(..., allow_detachement=False)` will erase the lightsim2grid kwargs `stop_if_storage_disco=False`")
            self._stop_if_storage_disco = True
                
    def _aux_set_correct_detach_flags(self):
        if self.detachment_is_allowed:
            self._aux_set_correct_detach_flags_d_allowed()
        else:
            self._aux_set_correct_detach_flags_d_not_allowed()
        
    def load_grid(self,
                  path : Union[os.PathLike, str],
                  filename : Optional[Union[os.PathLike, str]]=None) -> None: 
        cls = type(self)
        if hasattr(cls, "can_handle_more_than_2_busbar"):
            # grid2op version >= 1.10.0 then we use this
            self.can_handle_more_than_2_busbar()
            
        if hasattr(cls, "can_handle_detachment"):
            # grid2op version >= 1.11.0 then we use this
            self.can_handle_detachment()
            
        self._aux_set_correct_detach_flags()
            
        if self._loader_method == "pandapower":
            self._load_grid_pandapower(path, filename)
        elif self._loader_method == "pypowsybl":
            self._load_grid_pypowsybl(path, filename)
        else:
            raise BackendError(f"Impossible to initialize the backend with '{self._loader_method}'")
        self._grid.tell_solver_need_reset()
        self._reset_res_pointers()  # force the re reading of the accessors at next powerflow
        self.V = np.ones(self.nb_bus_total, dtype=np.complex_)
    
    def _should_not_have_to_do_this(self, path=None, filename=None):
        # included in grid2op now !
        # but before `make_complete_path` was introduced we need to still
        # be able to use lightsim2grid
        if path is None and filename is None:
            raise Grid2OpException(
                "You must provide at least one of path or file to load a powergrid."
            )
        if path is None:
            full_path = filename
        elif filename is None:
            full_path = path
        else:
            full_path = os.path.join(path, filename)
        if not os.path.exists(full_path):
            raise Grid2OpException('There is no powergrid at "{}"'.format(full_path))
        return full_path

    def _aux_get_substation_handling_from_loader_kwargs(self, loader_kwargs) -> Union[int, None]:
        res = None
        if "double_bus_per_sub" in loader_kwargs and loader_kwargs["double_bus_per_sub"]:
            res = DEFAULT_N_BUSBAR_PER_SUB
        if "n_busbar_per_sub" in loader_kwargs and loader_kwargs["n_busbar_per_sub"]:
            if res is not None:
                raise BackendError("When intializing a grid from pypowsybl, you cannot "
                                   "set both `double_bus_per_sub` and `n_busbar_per_sub` "
                                   "in the `loader_kwargs`. "
                                   "You can only set `n_busbar_per_sub` in this case.")
            res = int(loader_kwargs["n_busbar_per_sub"])
            if loader_kwargs["n_busbar_per_sub"] != res:
                raise BackendError("When initializing a grid from pypowsybl, the `n_busbar_per_sub` "
                                   "loader kwargs should be properly convertible to an int "
                                   "giving the default number of busbar sections per substation.")
        if res is not None:
            if self.n_busbar_per_sub != res:
                if self.n_busbar_per_sub != DEFAULT_N_BUSBAR_PER_SUB:
                    warnings.warn(f"You specified `n_busbar={self.n_busbar_per_sub}` when calling `grid2op.make(...)` "
                                  f"but also specified `n_busbar_per_sub` (or `double_bus_per_sub`) in the "
                                  f"`loader_kwargs` when building the `LightSimBackend` "
                                  f"(**eg** `LighSimBackend(..., loader_kwargs={{'n_busbar_per_sub': {res} }})`). "
                                  f"IN THIS CASE {self.n_busbar_per_sub} WILL BE USED (the things specified in the "
                                  f"`grid2op.make(...)`)")
                    res = self.n_busbar_per_sub
                elif res != DEFAULT_N_BUSBAR_PER_SUB:
                    warnings.warn(f"You specified `n_busbar_per_sub={res}` (or `double_bus_per_sub`) "
                                  f"in the `load_kwargs` when building the LightSimBackend "
                                  f"(**eg** `LighSimBackend(..., loader_kwargs={{'n_busbar_per_sub': {res} }})`) "
                                  f"without modifying the `n_busbar=...` when calling `grid2op.make(...)`. "
                                  f"LightSimBackend WILL OVERRIDE THE MAXIMUM NUMBER OF BUSES PER SUBSTATION "
                                  f"and set it to {res}. You can silence this warning by creating the env with "
                                  f"`grid2op.make(..., n_busbar={res})`")
                    self.n_busbar_per_sub = res
        return res
    
    def _get_subid_from_buses_legacy(self, buses_sub_id, el_sub_df):
        # buses_sub_id is the first element as returned by from_pypowsybl / init function
        # el_sub_df is an element dataframe returned by the same function
        tmp = pd.merge(el_sub_df.reset_index(), buses_sub_id, how="left", right_on="sub_id", left_on="sub_id")
        res = tmp.drop_duplicates(subset='id').set_index("id").sort_index()["bus_id"].values
        return res
    
    def _load_grid_pypowsybl(self, path=None, filename=None):
        from lightsim2grid.gridmodel.from_pypowsybl import init as init_pypow
        import pypowsybl.network as pypow_net
        loader_kwargs = {}
        if self._loader_kwargs is not None:
            loader_kwargs = self._loader_kwargs
        
        try:
            full_path = self.make_complete_path(path, filename)
        except AttributeError as exc_:
            warnings.warn("Please upgrade your grid2op version")
            full_path = self._should_not_have_to_do_this(path, filename)
            
        grid_tmp = pypow_net.load(full_path)
        self._orig_grid_pypowsybl = grid_tmp
        gen_slack_id = None
        if "gen_slack_id" in loader_kwargs:
            gen_slack_id = loader_kwargs["gen_slack_id"]
        
        df = grid_tmp.get_substations()
        buses_for_sub = True
        if "use_buses_for_sub" in loader_kwargs and loader_kwargs["use_buses_for_sub"]:
            df = grid_tmp.get_buses()
            buses_for_sub = False
            self.n_sub = df.shape[0]
            self.name_sub = ["sub_{}".format(i) for i, _ in enumerate(df.iterrows())]
            
        n_busbar_per_sub = self._aux_get_substation_handling_from_loader_kwargs(loader_kwargs)
        if n_busbar_per_sub is None:
            n_busbar_per_sub = self.n_busbar_per_sub
        self._grid, subs_id = init_pypow(grid_tmp,
                                         gen_slack_id=gen_slack_id,
                                         sort_index=True,
                                         return_sub_id=True,
                                         n_busbar_per_sub=n_busbar_per_sub,
                                         buses_for_sub=buses_for_sub,
                                         )
        (buses_sub_id, gen_sub, load_sub, (lor_sub, tor_sub), (lex_sub, tex_sub), 
         batt_sub, sh_sub, hvdc_sub_from_id, hvdc_sub_to_id) = subs_id
        if not buses_for_sub: 
            self.__nb_bus_before = self._grid.get_n_sub()
        else:
            self.__nb_bus_before = grid_tmp.get_buses().shape[0]
        self._aux_setup_right_after_grid_init()   
        
        # mandatory for the backend
        self.n_line = len(self._grid.get_lines()) + len(self._grid.get_trafos())
        self.n_gen = len(self._grid.get_generators())
        self.n_load = len(self._grid.get_loads())
        self.n_storage = len(self._grid.get_storages())
        if type(self).shunts_data_available:
            self.n_shunt = len(self._grid.get_shunts())
        else:
            self.n_shunt = None
        
        if not buses_for_sub:
            # consider that each "bus" in the powsybl grid is a substation
            # this is the "standard" behaviour for IEEE grid in grid2op
            # but can be considered "legacy" behaviour for more realistic grid
            this_load_sub = self._get_subid_from_buses_legacy(buses_sub_id, load_sub)
            this_gen_sub = self._get_subid_from_buses_legacy(buses_sub_id, gen_sub)
            this_lor_sub = self._get_subid_from_buses_legacy(buses_sub_id, lor_sub)
            this_tor_sub = self._get_subid_from_buses_legacy(buses_sub_id, tor_sub)
            this_lex_sub = self._get_subid_from_buses_legacy(buses_sub_id, lex_sub)
            this_tex_sub = self._get_subid_from_buses_legacy(buses_sub_id, tex_sub)
            this_batt_sub = self._get_subid_from_buses_legacy(buses_sub_id, batt_sub)
            this_sh_sub = self._get_subid_from_buses_legacy(buses_sub_id, sh_sub)
            
            self.load_to_subid = np.array(this_load_sub, dtype=dt_int)
            self.gen_to_subid = np.array(this_gen_sub, dtype=dt_int)
            self.line_or_to_subid = np.concatenate((this_lor_sub, this_tor_sub)).astype(dt_int)
            self.line_ex_to_subid = np.concatenate((this_lex_sub, this_tex_sub)).astype(dt_int)
            if self.__has_storage:
                self.storage_to_subid = np.array(this_batt_sub, dtype=dt_int)
            if self.n_shunt is not None:
                self.shunt_to_subid = np.array(this_sh_sub, dtype=dt_int)
        else:
            # consider effectively that each "voltage_levels" in powsybl grid
            # is a substation (in the underlying gridmodel)
            
            # TODO get back the sub id from the grid_tmp.get_voltage_levels()
            # need to work on that grid2op side: different make sure the labelling of the buses are correct !
 
            self.load_to_subid = np.array(load_sub.values.ravel(), dtype=dt_int)
            self.gen_to_subid = np.array(gen_sub.values.ravel(), dtype=dt_int)
            self.line_or_to_subid = np.concatenate((lor_sub.values.ravel(), tor_sub.values.ravel())).astype(dt_int)
            self.line_ex_to_subid = np.concatenate((lex_sub.values.ravel(), tex_sub.values.ravel())).astype(dt_int)
            if self.__has_storage:
                self.storage_to_subid = np.array(batt_sub.values.ravel(), dtype=dt_int)
            self.shunt_to_subid = np.array(sh_sub.values.ravel(), dtype=dt_int)
            if self.n_shunt is not None:
                self.n_sub = grid_tmp.get_voltage_levels().shape[0]
        
        # the names
        use_grid2op_default_names = True
        if "use_grid2op_default_names" in loader_kwargs and not loader_kwargs["use_grid2op_default_names"]:
            use_grid2op_default_names = False
            
        if use_grid2op_default_names:
            self.name_load = np.array([f"load_{el.bus_id}_{id_obj}" for id_obj, el in enumerate(self._grid.get_loads())]).astype(str)
            self.name_gen = np.array([f"gen_{el.bus_id}_{id_obj}" for id_obj, el in enumerate(self._grid.get_generators())]).astype(str)
            self.name_line = np.array([f"{el.bus_or_id}_{el.bus_ex_id}_{id_obj}"  for id_obj, el in enumerate(self._grid.get_lines())] +
                                    [f"{el.bus_hv_id}_{el.bus_lv_id}_{id_obj}"  for id_obj, el in enumerate(self._grid.get_trafos())]).astype(str)
            self.name_storage = np.array([f"storage_{el.bus_id}_{id_obj}"  for id_obj, el in enumerate(self._grid.get_storages())]).astype(str)
            self.name_shunt = np.array([f"shunt_{el.bus_id}_{id_obj}"  for id_obj, el in enumerate(self._grid.get_shunts())]).astype(str)
        else:
            self.name_load = np.array(load_sub.index.astype(str))
            self.name_gen = np.array(gen_sub.index.astype(str))
            self.name_line = np.concatenate((lor_sub.index.astype(str), tor_sub.index.astype(str)))
            self.name_storage = np.array(batt_sub.index.astype(str))
            self.name_shunt = np.array(sh_sub.index.astype(str))
            if not buses_for_sub:
                self.name_sub = np.array(buses_sub_id.index.astype(str))
            else:
                self.name_sub = np.array(grid_tmp.get_buses().loc[buses_sub_id.drop_duplicates("sub_id").index, "voltage_level_id"].values)

        # and now things needed by the backend (legacy)
        self.prod_pu_to_kv = 1.0 * self._grid.get_bus_vn_kv()[[el.bus_id for el in self._grid.get_generators()]]
        self.prod_pu_to_kv = self.prod_pu_to_kv.astype(dt_float)
        
        if "reconnect_disco_gen" in loader_kwargs and loader_kwargs["reconnect_disco_gen"]:
            for el in self._grid.get_generators():
                if not el.connected:
                    self._grid.reactivate_gen(el.id)
                    self._grid.change_bus_gen(el.id, self.gen_to_subid[el.id])
                    self._grid.change_p_gen(el.id, 0.)
                    self._grid.change_v_gen(el.id, self.prod_pu_to_kv[el.id])
                    
        if "reconnect_disco_load" in loader_kwargs and loader_kwargs["reconnect_disco_load"]:
            for el in self._grid.get_loads():
                if not el.connected:
                    self._grid.reactivate_load(el.id)
                    self._grid.change_bus_load(el.id, self.load_to_subid[el.id])
                    self._grid.change_p_load(el.id, 0.)
                    self._grid.change_q_load(el.id, 0.)
                    
        # complete the other vectors
        self._compute_pos_big_topo()
        
        self.__nb_powerline = len(self._grid.get_lines())
        
        # init this
        self.prod_p = np.array([el.target_p_mw for el in self._grid.get_generators()], dtype=dt_float)
        self.next_prod_p = np.array([el.target_p_mw for el in self._grid.get_generators()], dtype=dt_float)
        self.nb_bus_total = len(self._grid.get_bus_vn_kv())
        
        # and now things needed by the backend (legacy)
        self._big_topo_to_obj = [(None, None) for _ in range(type(self).dim_topo)]
        self.prod_pu_to_kv = 1.0 * self._grid.get_bus_vn_kv()[[el.bus_id for el in self._grid.get_generators()]]
        self.prod_pu_to_kv = self.prod_pu_to_kv.astype(dt_float)
        
        # TODO
        max_not_too_max = (np.finfo(dt_float).max * 0.5 - 1.)
        self.thermal_limit_a = max_not_too_max * np.ones(self.n_line, dtype=dt_float)
        bus_vn_kv = np.array(self._grid.get_bus_vn_kv())
        # shunt_bus_id = np.array([el.bus_id for el in self._grid.get_shunts()])
        self._sh_vnkv = bus_vn_kv[self.shunt_to_subid]
        self._aux_finish_setup_after_reading()
    
    def _aux_setup_right_after_grid_init(self):
        if  self._orig_grid_pypowsybl is None:
            self._grid.set_n_sub(self.__nb_bus_before)
        self._handle_turnedoff_pv()
            
        self.available_solvers = self._grid.available_solvers()
        if self.__current_solver_type is None:
            # previous default behaviour (< 0.7)
            # by default it builds the backend with the fastest solver
            # automatically found
            self._assign_right_solver()
            
            if SolverType.KLUDC in self.available_solvers:
                # use the faster KLU if available even for DC approximation
                self._grid.change_solver(SolverType.KLUDC)
                
            self.__current_solver_type = copy.deepcopy(self._grid.get_solver_type())
        else:
            # check that the solver type provided is installed with lightsim2grid
            self._check_suitable_solver_type(self.__current_solver_type)
            self._grid.change_solver(self.__current_solver_type)
                    
        # handle multiple busbar per substations
        if hasattr(type(self), "can_handle_more_than_2_busbar") and self._orig_grid_pypowsybl is None:
            # grid2op version >= 1.10.0 then we use this
            self._grid._max_nb_bus_per_sub = self.n_busbar_per_sub
            
        self._grid.tell_solver_need_reset()
    
    def init_from_loaded_pandapower(self, pp_net):
        if hasattr(type(self), "can_handle_more_than_2_busbar"):
            type(self.init_pp_backend).n_busbar_per_sub = self.n_busbar_per_sub
        self.init_pp_backend = pp_net.copy()
        self._aux_init_pandapower()
        
        # handles redispatching
        if type(pp_net).redispatching_unit_commitment_availble:
            self.redispatching_unit_commitment_availble = True
            for attr_nm in ["gen_type", "gen_pmin", "gen_pmax",
                            "gen_redispatchable", "gen_max_ramp_up",
                            "gen_max_ramp_down", "gen_min_uptime",
                            "gen_min_downtime", "gen_cost_per_MW",
                            "gen_startup_cost", "gen_shutdown_cost",
                            "gen_renewable"
                            ]:
                setattr(self, attr_nm, copy.deepcopy(getattr( type(pp_net), attr_nm)))
        
        # handles storages
        for attr_nm in ["storage_type",
                        "storage_Emax",
                        "storage_Emin",
                        "storage_max_p_prod" ,
                        "storage_max_p_absorb",
                        "storage_marginal_cost",
                        "storage_loss",
                        "storage_charging_efficiency",
                        "storage_discharging_efficiency",
                        ]:
            setattr(self, attr_nm, copy.deepcopy(getattr( type(pp_net), attr_nm)))
        
    def _load_grid_pandapower(self, path=None, filename=None):
        from lightsim2grid._utils import _DoNotUseAnywherePandaPowerBackend
        _DoNotUseAnywherePandaPowerBackend._clear_grid_dependant_class_attributes()
        if hasattr(type(self), "can_handle_more_than_2_busbar"):
            type(self.init_pp_backend).n_busbar_per_sub = self.n_busbar_per_sub
        type(self.init_pp_backend).set_env_name(type(self).env_name)
        if type(self).glop_version is not None:
            type(self.init_pp_backend).glop_version = type(self).glop_version
        self.init_pp_backend.load_grid(path, filename)
        self._aux_init_pandapower()
    
    def _aux_init_pandapower(self):
        from lightsim2grid.gridmodel import init_from_pandapower
        self._grid = init_from_pandapower(self.init_pp_backend._grid)    
        self.__nb_bus_before = self.init_pp_backend.get_nb_active_bus()  
        self._aux_setup_right_after_grid_init()        
        
        # deactive the buses that have been added
        for bus_id, bus_status in enumerate(self.init_pp_backend._grid.bus["in_service"].values):
            if bus_status:
                self._grid.reactivate_bus(bus_id)
            else:
                self._grid.deactivate_bus(bus_id)

        pp_cls = type(self.init_pp_backend)
        if pp_cls.n_line <= -1:
            warnings.warn("You are using a legacy (quite old now) grid2op version. Please upgrade it.")
            pp_cls = self.init_pp_backend
        self.n_line = pp_cls.n_line
        self.n_gen = pp_cls.n_gen
        self.n_load = pp_cls.n_load
        self.n_sub = pp_cls.n_sub
        self.sub_info = pp_cls.sub_info
        self.dim_topo = pp_cls.dim_topo
        self.load_to_subid = pp_cls.load_to_subid
        self.gen_to_subid = pp_cls.gen_to_subid
        self.line_or_to_subid = pp_cls.line_or_to_subid
        self.line_ex_to_subid = pp_cls.line_ex_to_subid
        self.load_to_sub_pos = pp_cls.load_to_sub_pos
        self.gen_to_sub_pos = pp_cls.gen_to_sub_pos
        self.line_or_to_sub_pos = pp_cls.line_or_to_sub_pos
        self.line_ex_to_sub_pos = pp_cls.line_ex_to_sub_pos
        
        self.name_gen = pp_cls.name_gen
        self.name_load = pp_cls.name_load
        self.name_line = pp_cls.name_line
        self.name_sub = pp_cls.name_sub

        self.prod_pu_to_kv = self.init_pp_backend.prod_pu_to_kv
        self.load_pu_to_kv = self.init_pp_backend.load_pu_to_kv
        self.lines_or_pu_to_kv = self.init_pp_backend.lines_or_pu_to_kv
        self.lines_ex_pu_to_kv = self.init_pp_backend.lines_ex_pu_to_kv

        # TODO storage check grid2op version and see if storage is available !
        if self.__has_storage:
            self.n_storage = pp_cls.n_storage
            self.storage_to_subid = pp_cls.storage_to_subid
            self.storage_pu_to_kv = self.init_pp_backend.storage_pu_to_kv
            self.name_storage = pp_cls.name_storage
            self.storage_to_sub_pos = pp_cls.storage_to_sub_pos
            self.storage_type = pp_cls.storage_type
            self.storage_Emin = pp_cls.storage_Emin
            self.storage_Emax = pp_cls.storage_Emax
            self.storage_max_p_prod = pp_cls.storage_max_p_prod
            self.storage_max_p_absorb = pp_cls.storage_max_p_absorb
            self.storage_marginal_cost = pp_cls.storage_marginal_cost
            self.storage_loss = pp_cls.storage_loss
            self.storage_discharging_efficiency = pp_cls.storage_discharging_efficiency
            self.storage_charging_efficiency = pp_cls.storage_charging_efficiency
        
        self.nb_bus_total = self.init_pp_backend._grid.bus.shape[0]

        self.thermal_limit_a = copy.deepcopy(self.init_pp_backend.thermal_limit_a)

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
        self.prod_p = 1.0 * self.init_pp_backend._grid.gen["p_mw"].values
        self.next_prod_p = 1.0 * self.init_pp_backend._grid.gen["p_mw"].values
        
        if type(self).shunts_data_available:
            if pp_cls.n_shunt is not None:
                # modern grid2op version
                self.n_shunt = pp_cls.n_shunt
                self.shunt_to_subid = pp_cls.shunt_to_subid
                self.name_shunt = pp_cls.name_shunt
            else:
                # legacy grid2op version...
                warnings.warn("You are using a legacy grid2op version, please upgrade grid2op.")
                self.n_shunt = self.init_pp_backend.n_shunt
                self.shunt_to_subid = self.init_pp_backend.shunt_to_subid
                self.name_shunt = self.init_pp_backend.name_shunt
        else:
            self.n_shunt = None
        self._compute_pos_big_topo()
        if hasattr(self.init_pp_backend, "_sh_vnkv"):
            # attribute has been added in grid2op ~1.3 or 1.4
            self._sh_vnkv = self.init_pp_backend._sh_vnkv
        
        self._aux_finish_setup_after_reading()

    def _aux_finish_setup_after_reading(self):
        # set up the "lightsim grid" accordingly
        cls = type(self)
        if cls.n_line <= -1:
            warnings.warn("You are using a legacy (quite old now) grid2op version. Please upgrade it.")
            cls = self
        self._grid.set_load_pos_topo_vect(cls.load_pos_topo_vect)
        self._grid.set_gen_pos_topo_vect(cls.gen_pos_topo_vect)
        self._grid.set_line_or_pos_topo_vect(cls.line_or_pos_topo_vect[:self.__nb_powerline])
        self._grid.set_line_ex_pos_topo_vect(cls.line_ex_pos_topo_vect[:self.__nb_powerline])
        self._grid.set_trafo_hv_pos_topo_vect(cls.line_or_pos_topo_vect[self.__nb_powerline:])
        self._grid.set_trafo_lv_pos_topo_vect(cls.line_ex_pos_topo_vect[self.__nb_powerline:])

        self._grid.set_load_to_subid(cls.load_to_subid)
        self._grid.set_gen_to_subid(cls.gen_to_subid)
        self._grid.set_line_or_to_subid(cls.line_or_to_subid[:self.__nb_powerline])
        self._grid.set_line_ex_to_subid(cls.line_ex_to_subid[:self.__nb_powerline])
        self._grid.set_trafo_hv_to_subid(cls.line_or_to_subid[self.__nb_powerline:])
        self._grid.set_trafo_lv_to_subid(cls.line_ex_to_subid[self.__nb_powerline:])

        # TODO storage check grid2op version and see if storage is available !
        if self.__has_storage:
            self._grid.set_storage_to_subid(cls.storage_to_subid)
            self._grid.set_storage_pos_topo_vect(cls.storage_pos_topo_vect)

        nm_ = "load"
        for load_id, pos_big_topo  in enumerate(cls.load_pos_topo_vect):
            self._big_topo_to_obj[pos_big_topo] = (load_id, nm_)
        nm_ = "gen"
        for gen_id, pos_big_topo  in enumerate(cls.gen_pos_topo_vect):
            self._big_topo_to_obj[pos_big_topo] = (gen_id, nm_)
        nm_ = "lineor"
        for l_id, pos_big_topo in enumerate(cls.line_or_pos_topo_vect):
            self._big_topo_to_obj[pos_big_topo] = (l_id, nm_)
        nm_ = "lineex"
        for l_id, pos_big_topo  in enumerate(cls.line_ex_pos_topo_vect):
            self._big_topo_to_obj[pos_big_topo] = (l_id, nm_)

        # TODO storage check grid2op version and see if storage is available !
        if self.__has_storage:
            nm_ = "storage"
            for l_id, pos_big_topo  in enumerate(cls.storage_pos_topo_vect):
                self._big_topo_to_obj[pos_big_topo] = (l_id, nm_)

        # number of object per bus, to activate, deactivate them
        self.nb_obj_per_bus = np.zeros(2 * self.__nb_bus_before, dtype=dt_int).reshape(-1)

        self.topo_vect = np.ones(cls.dim_topo, dtype=dt_int).reshape(-1)
        if cls.shunts_data_available:
            self.shunt_topo_vect = np.ones(cls.n_shunt, dtype=dt_int)
             # shunts
            self.sh_p = np.full(cls.n_shunt, dtype=dt_float, fill_value=np.NaN).reshape(-1)
            self.sh_q = np.full(cls.n_shunt, dtype=dt_float, fill_value=np.NaN).reshape(-1)
            self.sh_v = np.full(cls.n_shunt, dtype=dt_float, fill_value=np.NaN).reshape(-1)
            self.sh_theta = np.full(cls.n_shunt, dtype=dt_float, fill_value=np.NaN).reshape(-1)
            self.sh_bus = np.full(cls.n_shunt, dtype=dt_int, fill_value=-1).reshape(-1)

        self.p_or = np.full(cls.n_line, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.q_or = np.full(cls.n_line, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.v_or = np.full(cls.n_line, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.a_or = np.full(cls.n_line, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.p_ex = np.full(cls.n_line, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.q_ex = np.full(cls.n_line, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.v_ex = np.full(cls.n_line, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.a_ex = np.full(cls.n_line, dtype=dt_float, fill_value=np.NaN).reshape(-1)

        self.load_p = np.full(cls.n_load, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.load_q = np.full(cls.n_load, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.load_v = np.full(cls.n_load, dtype=dt_float, fill_value=np.NaN).reshape(-1)

        self.prod_p = np.full(cls.n_gen, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.prod_q = np.full(cls.n_gen, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.prod_v = np.full(cls.n_gen, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        
        self.line_or_theta = np.full(cls.n_line, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.line_ex_theta = np.full(cls.n_line, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.load_theta = np.full(cls.n_load, dtype=dt_float, fill_value=np.NaN).reshape(-1)
        self.gen_theta = np.full(cls.n_gen, dtype=dt_float, fill_value=np.NaN).reshape(-1)

        # TODO storage check grid2op version and see if storage is available !
        if self.__has_storage:
            self.storage_p = np.full(cls.n_storage, dtype=dt_float, fill_value=np.NaN).reshape(-1)
            self.storage_q = np.full(cls.n_storage, dtype=dt_float, fill_value=np.NaN).reshape(-1)
            self.storage_v = np.full(cls.n_storage, dtype=dt_float, fill_value=np.NaN).reshape(-1)
            self.storage_theta = np.full(cls.n_storage, dtype=dt_float, fill_value=np.NaN).reshape(-1)

        self._count_object_per_bus()
        
        # set the initial topology vector
        n_sub_cls_orig = cls.n_sub
        n_sub_ls_orig = LightSimBackend.n_sub
        try:
            cls.n_sub = self.n_sub
            LightSimBackend.n_sub = self.n_sub
            self.topo_vect[cls.load_pos_topo_vect] = cls.global_bus_to_local(np.array([el.bus_id for el in self._grid.get_loads()]),
                                                                            cls.load_to_subid)
            self.topo_vect[cls.gen_pos_topo_vect] = cls.global_bus_to_local(np.array([el.bus_id for el in self._grid.get_generators()]),
                                                                            cls.gen_to_subid)
            if self.__has_storage:
                self.topo_vect[cls.storage_pos_topo_vect] = cls.global_bus_to_local(np.array([el.bus_id for el in self._grid.get_storages()]),
                                                                                    cls.storage_to_subid)
            lor_glob_bus = np.concatenate((np.array([el.bus_or_id for el in self._grid.get_lines()]),
                                        np.array([el.bus_hv_id for el in self._grid.get_trafos()])))
            self.topo_vect[cls.line_or_pos_topo_vect] = cls.global_bus_to_local(lor_glob_bus,
                                                                                cls.line_or_to_subid)
            lex_glob_bus = np.concatenate((np.array([el.bus_ex_id for el in self._grid.get_lines()]),
                                        np.array([el.bus_lv_id for el in self._grid.get_trafos()])))
            self.topo_vect[cls.line_ex_pos_topo_vect] = cls.global_bus_to_local(lex_glob_bus,
                                                                                cls.line_ex_to_subid)
            
            self._grid.tell_solver_need_reset()
            self.__me_at_init = self._grid.copy()
            self.__init_topo_vect = np.ones(cls.dim_topo, dtype=dt_int)
            self.__init_topo_vect[:] = self.topo_vect
            if cls.shunts_data_available:
                self.sh_bus[:] = cls.global_bus_to_local(np.array([el.bus_id for el in self._grid.get_shunts()]),
                                                        cls.shunt_to_subid)
        finally:
            cls.n_sub = n_sub_cls_orig
            LightSimBackend.n_sub = n_sub_ls_orig
        
    def assert_grid_correct_after_powerflow(self) -> None:
        """
        This method is called by the environment. It ensure that the backend remains consistent even after a powerflow
        has be run with :func:`Backend.runpf` method.

        :return: ``None``
        :raise: :class:`grid2op.Exceptions.EnvError` and possibly all of its derived class.
        """
        # test the results gives the proper size
        super().assert_grid_correct_after_powerflow()
        self.init_pp_backend.__class__ = type(self.init_pp_backend).init_grid(type(self))
        from lightsim2grid._utils import _DoNotUseAnywherePandaPowerBackend  # lazy import
        _DoNotUseAnywherePandaPowerBackend._clear_grid_dependant_class_attributes()
        self._backend_action_class = _BackendAction.init_grid(type(self))
        self._init_action_to_set = self._backend_action_class()
        try:
            _init_action_to_set = self.get_action_to_set()
        except TypeError:
            # I am in legacy grid2op version...
            _init_action_to_set = _dont_use_get_action_to_set_legacy(self)
            
        self._init_action_to_set += _init_action_to_set
        if self.prod_pu_to_kv is not None:
            assert np.isfinite(self.prod_pu_to_kv).all()
        if self.load_pu_to_kv is not None:
            assert np.isfinite(self.load_pu_to_kv).all()
        if self.lines_or_pu_to_kv is not None:
            assert np.isfinite(self.lines_or_pu_to_kv).all()
        if self.lines_ex_pu_to_kv is not None:
            assert np.isfinite(self.lines_ex_pu_to_kv).all()
        if self.__has_storage and self.n_storage > 0 and self.storage_pu_to_kv is not None:
            assert np.isfinite(self.storage_pu_to_kv).all()

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

        if type(self).shunts_data_available:
            arr_ = self.shunt_topo_vect
            is_connected = arr_ > 0
            arr_ = self.shunt_to_subid[is_connected] + self.__nb_bus_before * (arr_[is_connected] - 1)
            self.nb_obj_per_bus[arr_] += 1

    def close(self) -> None:
        self.init_pp_backend.close()
        self._grid = None

    def apply_action(self, backendAction: Union["grid2op.Action._backendAction._BackendAction", None]) -> None:
        """
        Specific implementation of the method to apply an action modifying a powergrid in the pandapower format.
        """
        tick = time.perf_counter()
        self._next_pf_fails = None
        active_bus, *_, topo__, shunts__ = backendAction()

        # change the overall topology
        chgt = backendAction.current_topo.changed
        self._grid.update_topo(chgt, backendAction.current_topo.values)
        self.topo_vect[chgt] = backendAction.current_topo.values[chgt]
        
        # print(f" load p {backendAction.load_p.values[backendAction.load_p.changed]}")  # TODO DEBUG WINDOWS
        # print(f" prod_p p {backendAction.prod_p.values[backendAction.prod_p.changed]}")  # TODO DEBUG WINDOWS
        # update the injections
        try:
            self._grid.update_gens_p(backendAction.prod_p.changed,
                                     backendAction.prod_p.values)
            self._grid.update_gens_v(backendAction.prod_v.changed,
                                     backendAction.prod_v.values / self.prod_pu_to_kv)
            self._grid.update_loads_p(backendAction.load_p.changed,
                                      backendAction.load_p.values)
            self._grid.update_loads_q(backendAction.load_q.changed,
                                      backendAction.load_q.values)
        except RuntimeError as exc_:
            # see https://github.com/Grid2Op/lightsim2grid/issues/66 (even though it's not a "bug" and has not been replicated)
            raise BackendError(f"{exc_}") from exc_
        
        if self.__has_storage:
            # print(f"\t (in backend) storage_power {backendAction.storage_power.values[backendAction.storage_power.changed]}")  # TODO DEBUG WINDOWS
            try:
                self._grid.update_storages_p(backendAction.storage_power.changed,
                                             backendAction.storage_power.values)
            except RuntimeError as exc_:
                # modification of power of disconnected storage has no effect in lightsim2grid
                if self.detachment_is_allowed:
                    # a storage units is allowed to be disconnected in this case
                    pass
                else:
                    self._next_pf_fails = BackendError("Some storage units would be disconnected")

        # handle shunts
        if type(self).shunts_data_available:
            shunt_p, shunt_q, shunt_bus = backendAction.shunt_p, backendAction.shunt_q, backendAction.shunt_bus
            # shunt topology
            # (need to be done before to avoid error like "impossible to set reactive value of a disconnected shunt")
            for sh_id, new_bus in shunt_bus:
                if new_bus == -1:
                    self._grid.deactivate_shunt(sh_id)
                else:
                    self._grid.reactivate_shunt(sh_id)
                    if hasattr(type(self), "local_bus_to_global_int"):
                        self._grid.change_bus_shunt(sh_id, type(self).local_bus_to_global_int(new_bus, self.shunt_to_subid[sh_id]))
                    else:
                        self._grid.change_bus_shunt(sh_id, self.shunt_to_subid[sh_id] + (new_bus == 2) * type(self).n_sub)
            # remember the topology not to need to read it back from the grid
            self.sh_bus[shunt_bus.changed] = shunt_bus.values[shunt_bus.changed]
            for sh_id, new_p in shunt_p:
                self._grid.change_p_shunt(sh_id, new_p)
            for sh_id, new_q in shunt_q:
                self._grid.change_q_shunt(sh_id, new_q)
        
        self._handle_dist_slack()
        self._timer_apply_act += time.perf_counter() - tick
        
    def _handle_dist_slack(self):
        if self._dist_slack_non_renew:
            self._grid.update_slack_weights(type(self).gen_redispatchable)

    def _fetch_grid_data(self):
        beg_test = time.perf_counter()
        if self._lineor_res is None:
            self._lineor_res = self._grid.get_lineor_res_full()
        if self._lineex_res is None:
            self._lineex_res = self._grid.get_lineex_res_full()
        if self._load_res is None:
            self._load_res = self._grid.get_loads_res_full()
        if self._gen_res is None:
            self._gen_res = self._grid.get_gen_res_full()
        if self._trafo_hv_res is None:
            self._trafo_hv_res = self._grid.get_trafohv_res_full()
        if self._trafo_lv_res is None:
            self._trafo_lv_res = self._grid.get_trafolv_res_full()
        if self._storage_res is None:
            self._storage_res = self._grid.get_storages_res_full()
        if self._shunt_res is None:
            self._shunt_res = self._grid.get_shunts_res_full()
        self._timer_fetch_data_cpp += time.perf_counter() - beg_test
            
    def runpf(self, is_dc : bool=False) -> Tuple[bool, Union[Exception, None]]:
        my_exc_ = None
        res = False
        try:            
            if self._next_pf_fails is not None:
                raise self._next_pf_fails
            
            beg_preproc = time.perf_counter()
            if is_dc:
                # somehow, when asked to do a powerflow in DC, pandapower assign Vm to be
                # one everywhere...
                # But not when it initializes in DC mode... (see below)
                self.V = np.ones(self.nb_bus_total, dtype=np.complex_) #  * self._grid.get_init_vm_pu()
                tick = time.perf_counter()
                self._timer_preproc += tick - beg_preproc
                V = self._grid.dc_pf(self.V, self.max_it, self.tol)
                self._timer_solver += time.perf_counter() - tick
                if V.shape[0] == 0:
                    raise BackendError(f"Divergence of DC powerflow (non connected grid). Detailed error: {self._grid.get_dc_solver().get_error()}")
            else:
                if (self.V is None) or (self.V.shape[0] == 0):
                    # create the vector V as it is not created
                    self.V = np.ones(self.nb_bus_total, dtype=np.complex_) * self._grid.get_init_vm_pu()
                if self.initdc:
                    self._grid.deactivate_result_computation()
                    # if I init with dc values, it should depends on previous state
                    self.V[:] = 1. # self._grid.get_init_vm_pu()  # see issue 30
                    # apparently pandapower run a  "real" dc powerflow with vm_pu = 1
                    # when it initialize the AC powerflow, 
                    # print(f"\tLightSimBackend: self.V.shape = {self.V.shape}")  # TODO DEBUG WINDOWS
                    self._debug_Vdc = self._grid.dc_pf(copy.deepcopy(self.V), self.max_it, self.tol)
                    self._grid.reactivate_result_computation()
                    if self._debug_Vdc.shape[0] == 0:
                        raise BackendError(f"Divergence of DC powerflow (non connected grid) at the "
                                           f"initialization of AC powerflow. Detailed error: "
                                           f"{self._grid.get_dc_solver().get_error()}")
                    V_init = 1. * self._debug_Vdc
                else:
                    V_init = copy.deepcopy(self.V)
                tick = time.perf_counter()
                self._timer_preproc += tick - beg_preproc
                V = self._grid.ac_pf(V_init, self.max_it, self.tol)
                self._timer_solver += time.perf_counter() - tick
                if V.shape[0] == 0:
                    raise BackendError(f"Divergence of AC powerflow. Detailed error: {self._grid.get_solver().get_error()}")

            beg_postroc = time.perf_counter()
            if is_dc:
                self.comp_time += self._grid.get_dc_computation_time()
                self.timer_gridmodel_xx_pf += self._grid.timer_last_dc_pf
            else:
                self.comp_time += self._grid.get_computation_time()
                # NB get_computation_time returns "time_total_nr", which is
                # defined in the powerflow algorithm and not on the linear solver.
                # it takes into account everything needed to solve the powerflow
                # once everything is passed to the solver.
                # It does not take into account the time to format the data in the 
                # from the GridModel 
                
                self.timer_gridmodel_xx_pf += self._grid.timer_last_ac_pf
                # timer_gridmodel_xx_pf takes all the time within the gridmodel "ac_pf"
            
            beg_readback = time.perf_counter()
            self.V[:] = V
            self._fetch_grid_data()
                    
            (self.p_or[:self.__nb_powerline],
             self.q_or[:self.__nb_powerline],
             self.v_or[:self.__nb_powerline],
             self.a_or[:self.__nb_powerline],
             self.line_or_theta[:self.__nb_powerline]) = self._lineor_res
            (self.p_or[self.__nb_powerline:],
             self.q_or[self.__nb_powerline:],
             self.v_or[self.__nb_powerline:],
             self.a_or[self.__nb_powerline:],
             self.line_or_theta[self.__nb_powerline:]) = self._trafo_hv_res
            (self.p_ex[:self.__nb_powerline],
             self.q_ex[:self.__nb_powerline],
             self.v_ex[:self.__nb_powerline],
             self.a_ex[:self.__nb_powerline],
             self.line_ex_theta[:self.__nb_powerline]) = self._lineex_res
            (self.p_ex[self.__nb_powerline:],
             self.q_ex[self.__nb_powerline:],
             self.v_ex[self.__nb_powerline:],
             self.a_ex[self.__nb_powerline:],
             self.line_ex_theta[self.__nb_powerline:]) = self._trafo_lv_res

            self.a_or *= 1000.  # kA in lightsim, A expected in grid2op
            self.a_ex *= 1000.  # kA in lightsim, A expected in grid2op

            self.a_or[~np.isfinite(self.a_or)] = 0.
            self.v_or[~np.isfinite(self.v_or)] = 0.
            self.a_ex[~np.isfinite(self.a_ex)] = 0.
            self.v_ex[~np.isfinite(self.v_ex)] = 0.

            self.load_p[:], self.load_q[:], self.load_v[:], self.load_theta[:] = self._load_res
            self.prod_p[:], self.prod_q[:], self.prod_v[:], self.gen_theta[:] = self._gen_res

            if self.__has_storage:
                self.storage_p[:], self.storage_q[:], self.storage_v[:], self.storage_theta[:] = self._storage_res
                self.storage_v[self.storage_v == -1.] = 0.  # voltage is 0. for disconnected elements in grid2op
            self._timer_read_data_back +=  time.perf_counter() - beg_readback

            self.next_prod_p[:] = self.prod_p
            if self._stop_if_load_disco and ((~np.isfinite(self.load_v)).any() or (self.load_v <= 0.).any()):
                disco = (~np.isfinite(self.load_v)) | (self.load_v <= 0.)
                load_disco = disco.nonzero()[0]
                self._timer_postproc += time.perf_counter() - beg_postroc
                raise BackendError(f"At least one load is disconnected (check loads {load_disco})")
            if self._stop_if_gen_disco and ((~np.isfinite(self.prod_v)).any() or (self.prod_v <= 0.).any()):
                disco = (~np.isfinite(self.prod_v)) | (self.prod_v <= 0.)
                gen_disco = disco.nonzero()[0]
                self._timer_postproc += time.perf_counter() - beg_postroc
                raise BackendError(f"At least one generator is disconnected (check gen {gen_disco})")
            
            if self.__has_storage:
                sto_active = (np.abs(self.storage_p) > 0.)
                sto_act_disco = (((~np.isfinite(self.storage_v)) & sto_active).any() or 
                                 ((self.storage_v <= 0.) & sto_active).any()
                                )
                if self._stop_if_storage_disco and sto_act_disco:
                    disco = ((~np.isfinite(self.storage_v)) | (self.storage_v <= 0.)) & sto_active
                    sto_disco = disco.nonzero()[0]
                    self._timer_postproc += time.perf_counter() - beg_postroc
                    raise BackendError(f"At least one storage unit is disconnected (check gen {sto_disco})")
            # TODO storage case of divergence !

            if type(self).shunts_data_available:
                self._set_shunt_info()

            if (np.abs(self.line_or_theta) >= 1e6).any() or (np.abs(self.line_ex_theta) >= 1e6).any():
                raise BackendError(f"Some theta are above 1e6 which should not be happening !")
            res = True
            self._grid.unset_changes()
            self._timer_postproc += time.perf_counter() - beg_postroc
        except Exception as exc_:
            # of the powerflow has not converged, results are Nan
            self._grid.unset_changes()
            self._fill_nans()
            res = False
            my_exc_ = exc_
            if not isinstance(my_exc_, BackendError):
                my_exc_ = BackendError(f"Converted the error of type {type(my_exc_)}, message was: {my_exc_}")
            if is_dc:
                # set back the solver to its previous state
                self._grid.change_solver(self.__current_solver_type)
            self._grid.tell_solver_need_reset()

        # TODO grid2op compatibility ! (was a single returned element before storage were introduced)
        if self.__has_storage:
            res = res, my_exc_
        return res

    def _fill_nans(self):
        """fill the results vectors with nans"""
        self._next_pf_fails = None
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
        self.line_or_theta[:] = np.NaN
        self.line_ex_theta[:] = np.NaN
        self.load_theta[:] = np.NaN
        self.gen_theta[:] = np.NaN

        if type(self).shunts_data_available:
            self.sh_p[:] = np.NaN
            self.sh_q[:] = np.NaN
            self.sh_v[:] = np.NaN
            self.sh_theta[:] = np.NaN
            self.sh_bus[:] = -1

        if self.__has_storage:
            self.storage_p[:] = np.NaN
            self.storage_q[:] = np.NaN
            self.storage_v[:] = np.NaN
            self.storage_theta[:] = np.NaN
        self.V[:] = self._grid.get_init_vm_pu()  # reset the V to its "original" value (see issue 30)
        self._reset_res_pointers()
    
    def _reset_res_pointers(self):
        self._lineor_res  = None
        self._lineex_res = None
        self._load_res = None
        self._gen_res = None
        self._trafo_hv_res = None
        self._trafo_lv_res = None
        self._storage_res = None
        self._shunt_res = None
        
    def __deepcopy__(self, memo):
        result = self.copy()
        memo[id(self)] = result
        return result

    def copy(self) -> Self:
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

        ####################
        # res = copy.deepcopy(self)  # super slow
        res = type(self).__new__(type(self))
        # make sure to init the "base class"
        # in particular with "new" attributes in future grid2op Backend
        res._aux_init_super(self.detailed_infos_for_cascading_failures,
                            self._can_be_copied,
                            self.__current_solver_type,
                            self.max_it,
                            self.tol,
                            self._turned_off_pv,
                            self._dist_slack_non_renew,
                            self._use_static_gen,
                            self._loader_method,
                            self._loader_kwargs,
                            self._stop_if_load_disco,
                            self._stop_if_gen_disco,
                            self._stop_if_storage_disco)
        
        # for backward compat (attribute was not necessarily present in early grid2op)
        if not hasattr(res, "_can_be_copied"):
            res._can_be_copied = self._can_be_copied
            
        res.comp_time = self.comp_time
        res.timer_gridmodel_xx_pf = self.timer_gridmodel_xx_pf

        # copy the regular attribute
        res.__has_storage = self.__has_storage
        res.__current_solver_type = self.__current_solver_type  # forced here because of special `__`
        res.__nb_powerline = self.__nb_powerline
        res.__nb_bus_before = self.__nb_bus_before
        res.cst_1 = dt_float(1.0)
        li_regular_attr = ["comp_time", "can_output_theta", "_is_loaded",
                           "nb_bus_total", "initdc",
                           "_big_topo_to_obj", "dim_topo",
                           "_idx_hack_storage",
                           "_timer_preproc", "_timer_postproc", "_timer_solver",
                           "_timer_read_data_back", "_timer_apply_act",
                           "supported_grid_format", 
                           "max_it", "tol", "_turned_off_pv", "_dist_slack_non_renew",
                           "_use_static_gen", "_loader_method", "_loader_kwargs",
                           "_stop_if_load_disco", "_stop_if_gen_disco", "_stop_if_storage_disco",
                           "_timer_fetch_data_cpp", "_next_pf_fails"
                           ]
        for attr_nm in li_regular_attr:
            if hasattr(self, attr_nm):
                # this test is needed for backward compatibility with older grid2op version
                setattr(res, attr_nm, copy.deepcopy(getattr(self, attr_nm)))

        # copy the numpy array
        li_attr_npy = ["thermal_limit_a", "_sh_vnkv", "_init_bus_load", "_init_bus_gen",
                       "_init_bus_lor", "_init_bus_lex", "nb_obj_per_bus", "next_prod_p", "topo_vect",
                       "shunt_topo_vect", "V", "prod_pu_to_kv", "load_pu_to_kv", "lines_or_pu_to_kv",
                       "lines_ex_pu_to_kv", "storage_pu_to_kv",
                       "p_or", "q_or", "v_or", "a_or",
                       "p_ex", "q_ex", "v_ex", "a_ex",
                       "load_p", "load_q", "load_v",
                       "prod_p", "prod_q", "prod_v",
                       "storage_p", "storage_q", "storage_v",
                       "sh_p", "sh_q", "sh_v", "sh_bus", "sh_theta",
                       "line_or_theta", "line_ex_theta", "load_theta", "gen_theta", "storage_theta",   
                       "_debug_Vdc"                
                       ]
        for attr_nm in li_attr_npy:
            if hasattr(self, attr_nm):
                # this test is needed for backward compatibility with older grid2op version
                setattr(res, attr_nm, copy.deepcopy(getattr(self, attr_nm)))

        # copy class attribute for older grid2op version (did not use the class attribute)
        if version.parse(grid2op.__version__) <= grid2op_min_cls_attr_ver:
            cls_attr = ["env_name", 
                        "n_line", "n_gen", "n_load", "n_sub",
                        "name_gen", "name_load", "name_line", "name_sub", "name_storage",
                        "sub_info", "dim_topo", 
                        "load_to_subid", "gen_to_subid", "line_or_to_subid",
                        "line_ex_to_subid", "load_to_sub_pos", "line_or_to_sub_pos",
                        "line_ex_to_sub_pos", "load_pos_topo_vect", "gen_pos_topo_vect",
                        "line_or_pos_topo_vect", "line_ex_pos_topo_vect", "storage_pos_topo_vect",

                        "shunts_data_available", "shunt_to_subid", "n_shunt", "name_shunt",

                        "n_storage", "storage_to_subid", "storage_to_sub_pos", "storage_type", 
                        "storage_Emax", "storage_Emin", "storage_max_p_prod", "storage_max_p_absorb",
                        "storage_marginal_cost", "storage_loss", "storage_charging_efficiency", 
                        "storage_discharging_efficiency",

                        "alarms_area_names", "alarms_lines_area", "alarms_area_lines"
                        ] + type(self)._li_attr_disp

            for attr_nm in cls_attr:
                if hasattr(self, attr_nm):
                    # this test is needed for backward compatibility with other grid2op version
                    setattr(res, attr_nm, copy.deepcopy(getattr(self, attr_nm)))
        ###############

        # handle the most complicated
        res._grid = mygrid.copy()
        res.__me_at_init = __me_at_init.copy()  # this is const
        res.init_pp_backend = inippbackend  # this is const
        res._init_action_to_set = copy.deepcopy(self._init_action_to_set)
        res._backend_action_class = self._backend_action_class  # this is const
        res.__init_topo_vect = self.__init_topo_vect
        res.available_solvers = self.available_solvers
        res._orig_grid_pypowsybl = self._orig_grid_pypowsybl 
        res._reset_res_pointers()
        res._fetch_grid_data()
        # assign back "self" attributes
        self._grid = mygrid
        self.init_pp_backend = inippbackend
        self.__me_at_init = __me_at_init
        return res

    def get_line_status(self) -> np.ndarray:
        l_s = self._grid.get_lines_status()
        t_s = self._grid.get_trafo_status()
        return np.concatenate((l_s, t_s)).astype(dt_bool)

    def get_line_flow(self) -> np.ndarray:
        return self.a_or
    
    def get_topo_vect(self) -> np.ndarray:
        return 1 * self.topo_vect

    def generators_info(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self.cst_1 * self.prod_p, self.cst_1 * self.prod_q, self.cst_1 * self.prod_v

    def loads_info(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self.cst_1 * self.load_p, self.cst_1 * self.load_q, self.cst_1 * self.load_v

    def lines_or_info(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return self.cst_1 * self.p_or, self.cst_1 * self.q_or, self.cst_1 * self.v_or, self.cst_1 * self.a_or

    def lines_ex_info(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return self.cst_1 * self.p_ex, self.cst_1 * self.q_ex, self.cst_1 * self.v_ex, self.cst_1 * self.a_ex

    def storages_info(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if not self.__has_storage:
            raise RuntimeError("Storage units are not supported with your grid2op version. Please upgrade to "
                               "grid2op >1.5")
        return self.cst_1 * self.storage_p, self.cst_1 * self.storage_q, self.cst_1 * self.storage_v

    def shunt_info(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return self.cst_1 * self.sh_p, self.cst_1 * self.sh_q, self.cst_1 * self.sh_v, self.sh_bus

    def _compute_shunt_bus_with_compat(self, shunt_bus):
        cls = type(self)
        if hasattr(cls, "global_bus_to_local"):
            self.sh_bus[:] = cls.global_bus_to_local(shunt_bus, cls.shunt_to_subid)
        else:
            res = (1 * shunt_bus).astype(dt_int)  # make a copy
            if hasattr(cls, "n_busbar_per_sub"):
                n_busbar_per_sub = cls.n_busbar_per_sub
            else:
                # backward compat when this was not defined:
                n_busbar_per_sub = DEFAULT_N_BUSBAR_PER_SUB
            for i in range(n_busbar_per_sub):
                res[(i * self.n_sub <= shunt_bus) & (shunt_bus < (i+1) * self.n_sub)] = i + 1
            res[shunt_bus == -1] = -1
            self.sh_bus[:] = res
        
    def _set_shunt_info(self):
        tick = time.perf_counter()
        self.sh_p[:], self.sh_q[:], self.sh_v[:], self.sh_theta[:]  = self._shunt_res
        self.sh_v[self.sh_v == -1.] = 0.  # in grid2op disco element have voltage of 0. and -1.
        self._timer_read_data_back += time.perf_counter() - tick
        # self._timer_fetch_data_cpp += time.perf_counter() - tick
        
    def _disconnect_line(self, id_):
        self.topo_vect[self.line_ex_pos_topo_vect[id_]] = -1
        self.topo_vect[self.line_or_pos_topo_vect[id_]] = -1
        if id_ < self.__nb_powerline:
            self._grid.deactivate_powerline(id_)
        else:
            self._grid.deactivate_trafo(id_ - self.__nb_powerline)

    def get_current_solver_type(self) -> SolverType:
        return self.__current_solver_type

    def reset(self,
              path : Union[os.PathLike, str],
              grid_filename : Optional[Union[os.PathLike, str]]=None) -> None:
        self._fill_nans()
        self._grid = self.__me_at_init.copy()
        self._grid.unset_changes()
        self._grid.change_solver(self.__current_solver_type)
        self._handle_turnedoff_pv()
        self.comp_time = 0.
        self.timer_gridmodel_xx_pf = 0.
        self._timer_postproc = 0.
        self._timer_preproc = 0.
        self._timer_solver = 0.
        self._timer_read_data_back = 0.
        self._timer_fetch_data_cpp = 0.
        self._timer_apply_act = 0.
        self._grid.tell_solver_need_reset()
        self._reset_res_pointers()
        if type(self).shunts_data_available:
            self.sh_bus[:] = 1  # TODO self._compute_shunt_bus_with_compat(self._grid.get_all_shunt_buses())
        self.topo_vect[:] = self.__init_topo_vect  # TODO#


def _dont_use_global_bus_to_local_legacy(cls, global_bus: np.ndarray, to_sub_id: np.ndarray) -> np.ndarray:
    res = (1 * global_bus).astype(dt_int)  # make a copy
    assert cls.n_busbar_per_sub >= 1, f"cls.n_busbar_per_sub should be >=1, found {cls.n_busbar_per_sub}"
    assert cls.n_sub >= 1, f"cls.n_sub should be >=1, found {cls.n_sub}"
    for i in range(cls.n_busbar_per_sub):
        res[(i * cls.n_sub <= global_bus) & (global_bus < (i+1) * cls.n_sub)] = i + 1
    res[global_bus == -1] = -1
    return res
        
        
def _dont_use_get_action_to_set_legacy(self, *arg, **kwargs):
    warnings.warn("DEPRECATION: grid2op <=1.4 is not well supported with lightsim2grid. Lots of bugs have been"
                    "fixed since then. Please upgrade to grid2op >= 1.10",
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


if not hasattr(Backend, "global_bus_to_local"):
    # for legacy grid2op
    setattr(LightSimBackend, "global_bus_to_local", classmethod(_dont_use_global_bus_to_local_legacy))
