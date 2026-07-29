# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["ContingencyAnalysisCPP", "LimitViolation", "ViolationElementType",
           "LimitViolationType", "PreContingencyResult",
           "ContingencyResult", "SecurityAnalysisResult"]

import copy
from dataclasses import dataclass
from typing import List, Optional, Tuple
import numpy as np
from collections.abc import Iterable

from lightsim2grid.algorithm import AlgorithmType
from .lightsim2grid_cpp import (ContingencyAnalysisCPP, LimitViolation,
                                 ViolationElementType, LimitViolationType)

try:
    from lightsim2grid.lightSimBackend import LightSimBackend
    __all__.append("ContingencyAnalysis")
    GRID2OP_INSTALLED = True
except ImportError as exc_:  # noqa: F841
    # grid2op is not installed
    GRID2OP_INSTALLED = False


@dataclass
class PreContingencyResult:
    """Limit violations for the pre-contingency ("n", no disconnection) case.

    .. note::
        ``converged`` is always True here: `ContingencyAnalysisCPP.compute` raises a
        ``RuntimeError`` if the pre-contingency powerflow itself does not converge (every
        contingency is solved relative to this base case, so a diverging base case makes the
        whole analysis meaningless).
    """
    converged: bool
    limit_violations: List[LimitViolation]


@dataclass
class ContingencyResult:
    """Limit violations for a single simulated contingency.

    .. note::
        If ``converged`` is False, ``limit_violations`` contains exactly one ``LimitViolation``
        with ``element_type == ViolationElementType.GRID`` and ``violation_type`` either
        ``LimitViolationType.NOT_SIMULATED`` (a pre-check skipped this contingency, eg it splits
        the grid, without ever invoking the solver) or ``LimitViolationType.DIVERGENCE`` (the
        solver ran but did not converge).
    """
    element_ids: List[int]  #: branch ids (lines then trafos) disconnected by this contingency
    element_names: List[str]  #: names (`env.name_line`) of the elements disconnected by this contingency
    contingency_name: Optional[str]  #: user-supplied name, see `add_single_contingency`
    converged: bool
    limit_violations: List[LimitViolation]


@dataclass
class SecurityAnalysisResult:
    """Result of `ContingencyAnalysis.run` / `run_ac` / `run_dc`, modeled after pypowsybl's
    security analysis result."""
    pre_contingency_result: PreContingencyResult
    post_contingency_results: List[ContingencyResult]


class __ContingencyAnalysis(object):
    """
    This class allows to perform a "security analysis" from a given grid state.

    For now, you cannot change the grid state, and it only computes the security analysis with 
    current flows at origin of powerlines.
    
    Feel free to post a feature request if you want to extend it.

    This class is used in 4 phases:
    
    0) you create it from a grid2op environment (the grid topology will not be modified from this environment)
    1) you add some contingencies to simulate
    2) you start the simulation
    3) you read back the results
    
    
    Examples
    --------
    An example is given here
    
    .. code-block:: python

        import grid2op
        from lightsim2grid import ContingencyAnalysis
        from lightsim2grid import LightSimBackend
        env_name = ...
        env = grid2op.make(env_name, backend=LightSimBackend())

        0) you create
        security_analysis = ContingencyAnalysis(env)
        
        1) you add some contingencies to simulate
        security_analysis.add_multiple_contingencies(...) # or security_analysis.add_single_contingency(...)
        
        2) you start the simulation (done automatically)
        3) you read back the results
        res_p, res_a, res_v = security_analysis.get_flows()

        # in this results, then
        # res_a[row_id] will be the flows, on all powerline corresponding to the `row_id` contingency.
        # you can retrieve it with `security_analysis.contingency_order[row_id]`

    Notes
    ------

    Sometimes, the behaviour might differ from grid2op. For example, if simulating a contingency
    leads to a non connected grid, then this function will return "Nan" for the flows and 0. for
    the voltages.

    In grid2op, it would be, in this case, 0. for the flows and 0. for the voltages.

    By default, a contingency that splits the grid in multiple connected components is not
    simulated (its voltages are left at 0.). If you set the `handle_disconnected_grid` attribute
    to ``True``, such contingencies are instead simulated on their largest connected component:
    the buses of the other component(s) are "masked" and their voltage is reported as 0. This is
    done without triggering any extra matrix re-factorization (the symbolic factorization of the
    solver is reused). It is supported by the Newton-Raphson family (AC) and by the DC solver; a
    non Newton-Raphson AC algorithm (*eg* Gauss-Seidel or Fast-Decoupled) is rejected.

    """
    STR_TYPES = (str, np.str_)  # np.str deprecated in numpy 1.20 and earlier versions not supported anyway
        
    def __init__(self, grid2op_env, compute_limit_violations: bool = False):
        if not GRID2OP_INSTALLED:
            raise RuntimeError("Impossible to use the python wrapper `ContingencyAnalysis` "
                               "when grid2op is not installed. Please fall back to the "
                               "c++ version (available in python) with:\n"
                               "\tfrom lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP\n"
                               "and refer to the appropriate documentation.")
        from grid2op.Environment import Environment # type: ignore
        self.__is_closed = False
        if isinstance(grid2op_env, Environment):    
            if not isinstance(grid2op_env.backend, LightSimBackend):
                raise RuntimeError("This class only works with LightSimBackend")
            self._ls_backend = grid2op_env.backend.copy()
        elif isinstance(grid2op_env, LightSimBackend):
            if hasattr(grid2op_env, "_is_loaded") and not grid2op_env._is_loaded:
                raise RuntimeError("Impossible to init a `ContingencyAnalysis` "
                                   "with a backend that has not been initialized.")
            self._ls_backend = grid2op_env.copy()
        else:
            raise RuntimeError("`ContingencyAnalysis` can only be created "
                               "with a grid2op `Environment` or a `LightSimBackend`")
        self.computer = ContingencyAnalysisCPP(self._ls_backend._grid, bool(compute_limit_violations))
        self._contingency_order = {}  # key: contingency (as tuple), value: order in which it is entered
        self._all_contingencies = []
        self._contingency_names = {}  # key: contingency (as tuple), value: user-supplied name (or None)
        self.__computed = False
        self._vs = None
        self._ampss = None
        self._mws = None

        self.available_default_algorithms = self.computer.available_default_algorithms()
        if AlgorithmType.NR_KLU in self.available_default_algorithms:
            # use the faster KLU if available
            self.computer.change_algorithm(AlgorithmType.NR_KLU)

    @property
    def all_contingencies(self):
        return copy.deepcopy(self._all_contingencies)

    @all_contingencies.setter
    def all_contingencies(self, val):
        raise RuntimeError("Impossible to add new topologies like this. Please use `add_single_contingency` "
                           "or `add_multiple_contingencies`.")

    @property
    def init_from_n_powerflow(self):
        return self.computer.init_from_n_powerflow
    
    @init_from_n_powerflow.setter
    def init_from_n_powerflow(self, val: bool):
        if bool(val) != val:
            raise ValueError("The `init_from_n_powerflow` attribute must be a boolean.")
        self.computer.init_from_n_powerflow = bool(val)

    @property
    def handle_disconnected_grid(self):
        return self.computer.handle_disconnected_grid

    @handle_disconnected_grid.setter
    def handle_disconnected_grid(self, val: bool):
        if bool(val) != val:
            raise ValueError("The `handle_disconnected_grid` attribute must be a boolean.")
        self.computer.handle_disconnected_grid = bool(val)

    @property
    def compute_limit_violations(self):
        """Whether limit violations are computed inline, per contingency, during `run` /
        `run_ac` / `run_dc` (see also `get_violations` on the underlying `computer`). Default:
        false. Computing violations means an extra per-element current / voltage check in
        every contingency's solve, so leave this off if you only need `get_flows`. Can only be
        set at construction time (`ContingencyAnalysis(env, compute_limit_violations=True)`) or
        via this setter; changing it clears any previously-computed results.
        """
        return self.computer.compute_limit_violations

    @compute_limit_violations.setter
    def compute_limit_violations(self, val: bool):
        if bool(val) != val:
            raise ValueError("The `compute_limit_violations` attribute must be a boolean.")
        val = bool(val)
        if val == self.computer.compute_limit_violations:
            return  # no-op, matches the C++ side (which also no-ops and does not clear)
        self.computer.compute_limit_violations = val  # this clears the C++-side results / contingencies
        self.clear()  # keep the python-side bookkeeping (contingency order / names) in sync

    @property
    def nb_thread(self):
        """Number of OS threads used to solve the contingencies (default: 1).

        With ``nb_thread == 1`` the behaviour is identical to the legacy
        sequential computation. With ``nb_thread > 1`` the contingencies are
        split across that many threads (each with its own solver and admittance
        matrix copy); the results do not depend on the number of threads.
        """
        return self.computer.nb_thread

    @nb_thread.setter
    def nb_thread(self, val: int):
        if int(val) != val:
            raise ValueError("The `nb_thread` attribute must be an integer.")
        self.computer.nb_thread = int(val)

    # TODO implement that !
    def __update_grid(self, backend_act):
        raise NotImplementedError("TODO !")
        self.clear(with_contlist=False)
        self._ls_backend.apply_action(backend_act)
        # run the powerflow
        self._ls_backend.runpf()
        # update the computer
        # self.computer = ContingencyAnalysisCPP(self._ls_backend._grid)
        # self.computer.update_grid(...)  # not implemented
        
    def clear(self, with_contlist=True):
        """
        Clear the list of contingencies to simulate
        """
        
        if self.__is_closed:
            raise RuntimeError("This is closed, you cannot use it.")

        self.__computed = False
        self._vs = None
        self._ampss = None
        self._mws = None
        if with_contlist:
            self.computer.clear()
            self._contingency_order = {}
            self._all_contingencies = []
            self._contingency_names = {}
        else:
            self.computer.clear_results_only()

    def _single_cont_to_li_int(self, single_cont):
        li_disc = []
        if isinstance(single_cont, int):
            single_cont = [single_cont]

        for stuff in single_cont:
            if isinstance(stuff, type(self).STR_TYPES):
                stuff = (type(self._ls_backend).name_line == stuff).nonzero()
                stuff = stuff[0]
                if stuff.size == 0:
                    # name is not found
                    raise RuntimeError(f"Impossible to find a powerline named \"{stuff}\" in the environment")
                stuff = int(stuff[0])
            else:
                stuff = int(stuff)
            li_disc.append(stuff)
        return li_disc

    def add_single_contingency(self, *args, name=None):
        """
        This function allows to add a single contingency specified by either the powerlines names
        (which should match env.name_line) or by their ID.

        The contingency added can be a "n-1" which will simulate a single powerline disconnection
        or a "n-k" which will simulate the disconnection of multiple powerlines.

        It does not accept any positional keyword arguments, but accepts the keyword-only
        `name` argument: an optional, user-supplied string used to identify this contingency
        in the result of `run` / `run_ac` / `run_dc` (`ContingencyResult.contingency_name`). If
        not given, `contingency_name` is `None` for this contingency. If this exact contingency
        was already registered (same set of disconnected elements), `name` is ignored (the
        first registration wins).

        Examples
        --------

        .. code-block:: python

            import grid2op
            from lightsim2grid import ContingencyAnalysis
            from lightsim2grid import LightSimBackend
            env_name = ...
            env = grid2op.make(env_name, backend=LightSimBackend())

            security_analysis = ContingencyAnalysis(env)
            # the single (n-1) contingency "disconnect powerline 0" is added
            security_analysis.add_single_contingency(0)

            # add the single (n-1) contingency "disconnect line 1
            security_analysis.add_single_contingency(env.name_line[1])

            # add a single contingency that disconnect powerline 2 and 3 at the same time
            security_analysis.add_single_contingency(env.name_line[2], 3)

        Notes
        -----
        If it raises an error for a given contingency, the object might be not properly initialized.
        In this case, we recommend you to clear it (using the `clear()` method and to attempt to 
        add contingencies again.)

        """
        if self.__is_closed:
            raise RuntimeError("This is closed, you cannot use it.")
        
        li_disc = self._single_cont_to_li_int(args)
        li_disc_tup = tuple(li_disc)
        if li_disc_tup not in self._contingency_order:
            # this is really the first time this contingency is seen
            try:
                self.computer.add_nk(li_disc)
                my_id = len(self._contingency_order)
                self._contingency_order[li_disc_tup] = my_id
                self._all_contingencies.append(li_disc_tup)
                self._contingency_names[li_disc_tup] = name
            except Exception as exc_:
                raise RuntimeError(f"Impossible to add the contingency {args}. The most likely cause "
                                   f"is that you try to disconnect a powerline that is not present "
                                   f"on the grid") from exc_

    def add_multiple_contingencies(self, *args):
        """
        This function will add multiple contingencies at the same time.

        This code is equivalent to:

        .. code-block:: python

            for single_cont in args:
                self.add_single_contingency(single_cont)

        It does not accept any keword arguments.

        Examples
        --------

        .. code-block:: python

            import grid2op
            from lightsim2grid import ContingencyAnalysis
            from lightsim2grid import LightSimBackend
            env_name = ...
            env = grid2op.make(env_name, backend=LightSimBackend())

            security_analysis = ContingencyAnalysis(env)

            # add a single contingency that disconnect powerline 2 and 3 at the same time
            security_analysis.add_single_contingency(env.name_line[2], 3)

            # add a multiple contingencies the first one disconnect powerline 2 and 
            # and the second one disconnect powerline 3
            security_analysis.add_multiple_contingencies(env.name_line[2], 3)
        """     
        if self.__is_closed:
            raise RuntimeError("This is closed, you cannot use it.")
           
        for single_cont in args:
            if isinstance(single_cont, Iterable) and not isinstance(single_cont, type(self).STR_TYPES):
                # this is a contingency consisting in cutting multiple powerlines
                self.add_single_contingency(*single_cont)
            else:
                # this is likely an int or a string representing a contingency
                self.add_single_contingency(single_cont)

    def add_all_n1_contingencies(self):
        """
        This method registers as the contingencies that will be computed all the contingencies that disconnects 1 powerline

        This is equivalent to:

        .. code-block:: python

            for single_cont_id in range(env.n_line):
                self.add_single_contingency(single_cont_id)
        """
        if self.__is_closed:
            raise RuntimeError("This is closed, you cannot use it.")
        
        for single_cont_id in range(type(self._ls_backend).n_line):
            self.add_single_contingency(single_cont_id)

    def get_flows(self, *args):
        """
        Retrieve the flows after each contingencies has been simulated.

        Each row of the resulting flow matrix will correspond to a contingency simulated in the arguments.

        You can require only the result on some contingencies with the `args` argument, but in each case, all the results will
        be computed. If you don't specify anything, the results will be returned for all contingencies (which we recommend to do)
        
        Examples
        --------

        .. code-block:: python

            import grid2op
            from lightsim2grid import ContingencyAnalysis
            from lightsim2grid import LightSimBackend
            env_name = ...
            env = grid2op.make(env_name, backend=LightSimBackend())

            security_analysis = ContingencyAnalysis(env)
            security_analysis.add_multiple_contingencies(...) # or security_analysis.add_single_contingency(...)
            res_p, res_a, res_v = security_analysis.get_flows()

            # in this results, then
            # res_a[row_id] will be the flows, on all powerline corresponding to the `row_id` contingency.
            # you can retrieve it with `security_analysis.contingency_order[row_id]`
        """
        if self.__is_closed:
            raise RuntimeError("This is closed, you cannot use it.")
        
        all_defaults = self.computer.my_defaults()
        if len(args) == 0:
            # default: i consider all contingencies
            orders_ = np.zeros(len(all_defaults), dtype=int)
            for id_cpp, cont_ in enumerate(all_defaults):
                tup_ = tuple(cont_)
                orders_[self._contingency_order[tup_]] = id_cpp
        else:
            # a list of interesting contingencies has been provided
            orders_ = np.zeros(len(args), dtype=int)
            all_defaults = [tuple(cont) for cont in all_defaults]
            for id_me, cont_ in enumerate(args):
                cont_li = self._single_cont_to_li_int(cont_)
                tup_ = tuple(cont_li)
                if tup_ not in self._contingency_order:
                    raise RuntimeError(f"Contingency {cont_} is not simulated by this class. Have you called "
                                       f"`add_single_contingency` or `add_multiple_contingencies` ?")
                id_cpp = all_defaults.index(tup_)
                orders_[id_me] = id_cpp

        if not self.__computed:
            self.compute_V()
            self.compute_A()
            self.compute_P()
        
        return self._mws[orders_], self._ampss[orders_], self._vs[orders_]

    def run(self) -> SecurityAnalysisResult:
        """
        Run this contingency analysis and report, for the pre-contingency ("n") case and for
        each registered contingency, the list of limit violations (bus voltage out of
        [vmin_kv, vmax_kv], line/trafo current above limit_a1_ka / limit_a2_ka -- see
        `LSGrid.set_bus_voltage_limits` / `set_line_current_limit_side1` / `set_line_current_limit_side2`
        / `set_trafo_current_limit_side1` / `set_trafo_current_limit_side2`).

        This requires `compute_limit_violations=True` (either passed at construction time or
        set via ``this_instance.compute_limit_violations = True``), else a `RuntimeError` is
        raised. Prefer `run_ac` / `run_dc` if you want to also select the algorithm family.

        A `RuntimeError` is also raised if the pre-contingency ("n") powerflow itself does not
        converge -- every contingency is solved relative to this base case, so a diverging base
        case makes the whole analysis meaningless.

        The returned object mimics pypowsybl's security analysis result:

        .. code-block:: python

            res = security_analysis.run()
            for v in res.pre_contingency_result.limit_violations:
                ...
            for cont in res.post_contingency_results:  # a list, ordered like `add_single_contingency` calls
                cont.element_ids       # branch ids disconnected by this contingency (always present)
                cont.element_names     # names (env.name_line) of these same elements (always present)
                cont.contingency_name  # optional, user-supplied via add_single_contingency(..., name=...)
                cont.converged
                cont.limit_violations

        .. note::
            A `converged == False` post-contingency entry has exactly one `LimitViolation` in
            `limit_violations`, with `element_type == ViolationElementType.GRID` and
            `violation_type` either `LimitViolationType.NOT_SIMULATED` (a pre-check skipped it,
            eg it splits the grid) or `LimitViolationType.DIVERGENCE` (the solver ran but did not
            converge) -- unlike a `converged == True` entry with an empty list, which genuinely
            means no violation was found.
        """
        if self.__is_closed:
            raise RuntimeError("This is closed, you cannot use it.")
        if not self.computer.compute_limit_violations:
            raise RuntimeError("`run` (and `run_ac` / `run_dc`) require `compute_limit_violations=True`, "
                               "set either at construction time (`ContingencyAnalysis(env, "
                               "compute_limit_violations=True)`) or via "
                               "`this_instance.compute_limit_violations = True`.")
        if not self.__computed:
            self.compute_V()

        all_defaults = self.computer.my_defaults()
        orders_ = np.zeros(len(all_defaults), dtype=int)
        for id_cpp, cont_ in enumerate(all_defaults):
            tup_ = tuple(cont_)
            orders_[self._contingency_order[tup_]] = id_cpp

        converged = self.computer.converged()
        violations = self.computer.get_violations()

        pre_contingency_result = PreContingencyResult(
            converged=self.computer.converged_n(),
            limit_violations=list(self.computer.get_violations_n()),
        )
        post_contingency_results = [
            ContingencyResult(
                element_ids=list(all_defaults[id_cpp]),
                element_names=[str(self._ls_backend.name_line[el_id]) for el_id in all_defaults[id_cpp]],
                contingency_name=self._contingency_names.get(self._all_contingencies[id_me]),
                converged=bool(converged[id_cpp]),
                limit_violations=list(violations[id_cpp]),
            )
            for id_me, id_cpp in enumerate(orders_)
        ]
        return SecurityAnalysisResult(pre_contingency_result, post_contingency_results)

    def _change_algorithm_family(self, want_dc: bool):
        """internal: switch to a default AC / DC algorithm, but only if the current one is not
        already of the requested family, so this does not needlessly `clear()` (and thus lose
        the registered contingencies) when it is already the case."""
        current_is_dc = "DC_" in self.computer.get_algo_type().name
        if current_is_dc == want_dc:
            return
        preferred = (["DC_KLU", "DC_SparseLU"] if want_dc else ["NR_KLU", "NR_SparseLU"])
        available = {a.name: a for a in self.computer.available_default_algorithms()}
        for name in preferred:
            if name in available:
                self.change_algorithm(available[name])
                return
        raise RuntimeError(f"Impossible to find a default {'DC' if want_dc else 'AC'} algorithm "
                           f"among {list(available.keys())}.")

    def run_ac(self) -> SecurityAnalysisResult:
        """Like `run`, but first makes sure an AC algorithm is selected (switches to NR_KLU /
        NR_SparseLU if the current algorithm is a DC one -- which clears any previously
        registered contingency, exactly like `change_algorithm`; does nothing if already AC)."""
        self._change_algorithm_family(want_dc=False)
        return self.run()

    def run_dc(self) -> SecurityAnalysisResult:
        """Like `run`, but first makes sure the DC algorithm is selected (switches to DC_KLU /
        DC_SparseLU if the current algorithm is an AC one -- which clears any previously
        registered contingency, exactly like `change_algorithm`; does nothing if already DC)."""
        self._change_algorithm_family(want_dc=True)
        return self.run()

    def compute_V(self):
        """
        This function allows to retrieve the complex voltage at each bus of the grid for each contingency.

        .. warning:: Order of the results

            The order in which the results are returned is NOT necessarily the order in which the contingencies have
            been entered. Please use `get_flows()` method for easier reading back of the results

        """
        if self.__is_closed:
            raise RuntimeError("This is closed, you cannot use it.")
        if self._ls_backend._debug_Vdc is not None:
            v_init = self._ls_backend._debug_Vdc.copy()
        else:
            v_init = self._ls_backend.V.copy()
        self.computer.compute(v_init,
                              self._ls_backend.max_it,
                              self._ls_backend.tol)
        self._vs = self.computer.get_voltages().copy()
        self.__computed = True
        return self._vs

    def compute_A(self):
        """
        This function returns the current flows (in Amps, A) at the origin / high voltage side

        .. warning:: Order of the results

            The order in which the results are returned is NOT necessarily the order in which the contingencies have
            been entered. Please use `get_flows()` method for easier reading back of the results !

        """
        if self.__is_closed:
            raise RuntimeError("This is closed, you cannot use it.")
        
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute_V has been sucessfully called")
        self.computer.compute_flows()
        self._ampss = 1e3 * self.computer.get_flows()
        return self._ampss

    def compute_P(self):
        """
        This function returns the active power flows (in MW) at the origin / high voltage side

        .. warning:: Order of the results

            The order in which the results are returned is NOT necessarily the order in which the contingencies have
            been entered. Please use `get_flows()` method for easier reading back of the results !

        """
        if self.__is_closed:
            raise RuntimeError("This is closed, you cannot use it.")
        
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute_V has been sucessfully called")
        self.computer.compute_power_flows()
        self._mws = self.computer.get_power_flows().copy()
        return self._mws

    def close(self):
        """permanently close the object"""
        if self.__is_closed:
            return
        self.clear()
        self.computer.close()
        self._ls_backend.close()
        self.__is_closed = True
        
    def change_algorithm(self, solver_type):
        if self.__is_closed:
            raise RuntimeError("This is closed, you cannot use it.")
        self.computer.change_algorithm(solver_type)
        self.clear()

    def change_solver(self, solver_type):
        # kept as a backward-compatible alias of `change_algorithm`
        self.change_algorithm(solver_type)
        
        
if GRID2OP_INSTALLED:
    ContingencyAnalysis = __ContingencyAnalysis
