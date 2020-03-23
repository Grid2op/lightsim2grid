from time import time
import numpy as np
from scipy import sparse

from pandapower.auxiliary import _add_auxiliary_elements
from pandapower.pd2ppc import _pd2ppc
from pandapower.results import reset_results, verify_results
from pandapower.pf.ppci_variables import _get_pf_variables_from_ppci, _store_results_from_pf_in_ppci
from pandapower.powerflow import _ppci_to_net
from pandapower.pf.run_newton_raphson_pf import ppci_to_pfsoln, _get_Y_bus, _get_Sbus, _store_internal
from pandapower.pf.run_newton_raphson_pf import _get_numba_functions, _run_dc_pf
from pandapower.run import _passed_runpp_parameters, _init_runpp_options, _check_bus_index_and_print_warning_if_high
from pandapower.run import _check_gen_index_and_print_warning_if_high
# from pandapower.run_newton_raphson_pf import ppci_to_pfsoln, _get_Y_bus, _get_Sbus, _store_internal, _get_numba_functions

from pyklu_cpp import KLUSolver, DataModel

import pdb

# columns of ppci id to their meaning
ID2Colname = {
    "bus": ["bus_i", "type", "Pd", "Qd", "Gs", "Bs", "area", "Vm", "Va", "baseKV", "zone", "Vmax", "Vmin"],
    "branch": ["fbus", "tbus", "r", "x", "b", "rateA", "rateB", "rateC", "ratio", "angle", "status", "angmin", "angmax"],
    "gen": ["bus", "Pg", "Qg", "Qmax", "Qmin", "Vg", "mBase", "status", "Pmax", "Pmin", "Pc1", "Pc2", "Qc1min",
            "Qc1max", "Qc2min", "Qc2max", "ramp_agc", "ramp_10", "ramp_30", "ramp_q", "apf"],
}
# meaning of the ppc columns to their id
ColID2Names = {key: {col: i for i, col in enumerate(val)} for key, val in ID2Colname.items()}

# baseR = np.square(base_kv) / net.sn_mva
# branch[f:t, BR_R] = line["r_ohm_per_km"].values * length_km / baseR / parallel
# branch[f:t, BR_X] = line["x_ohm_per_km"].values * length_km / baseR / parallel


# TODO optim: when i reuse the powerflow, i just update Sbus and reuse V and all the other stuff
# especially, i don't re create a solver etc.
# TODO just i want to test


def newtonpf(Ybus, V, Sbus, pv, pq, ppci, options):

    max_it = options["max_iteration"]
    tol = options['tolerance_mva']
    # initialize the solver
    solver = KLUSolver()
    Ybus = sparse.csc_matrix(Ybus)
    # do the newton raphson algorithm
    solver.solve(Ybus, V, Sbus, pv, pq, max_it, tol)

    # extract the results
    Va = solver.get_Va()
    Vm = solver.get_Vm()
    V = Vm * np.exp(1j * Va)
    J = solver.get_J()
    converged = solver.converged()
    iterations = solver.get_nb_iter()
    return V, converged, iterations, J


class KLU4Pandapower():
    def __init__(self):
        self.solver = KLUSolver()
        self.ppci = None
        self.V = None
        self.Ybus = None
        self.pv = None
        self.pq = None

        # TODO
        self.ppci = None
        self.Ybus = None
        self.Yf = None
        self.Yt = None
        self.bus = None
        self.gen = None
        self.branch = None
        self.ref = None
        self.pv = None
        self.pq = None
        self.ref_gens = None
        self.baseMVA = None

    def runpp(self, net, max_iteration=10,  need_reset=True, **kwargs):
        # ---------- pp.run.runpp() -----------------
        t0_start = time()

        t0_options = time()
        passed_parameters = _passed_runpp_parameters(locals())
        _init_runpp_options(net, algorithm="nr", calculate_voltage_angles="auto",
                            init="auto", max_iteration=max_iteration, tolerance_mva=1e-8,
                            trafo_model="t", trafo_loading="current",
                            enforce_q_lims=False, check_connectivity=False,
                            voltage_depend_loads=True,
                            consider_line_temperature=False,
                            passed_parameters=passed_parameters, numba=True, **kwargs)
        _check_bus_index_and_print_warning_if_high(net)
        _check_gen_index_and_print_warning_if_high(net)
        et_options = time() - t0_options

        # ---------- pp.powerflow._powerflow() -----------------
        """
        Gets called by runpp or rundcpp with different arguments.
        """
        # get infos from options
        t0_early_init = time()
        init_results = net["_options"]["init_results"]
        ac = net["_options"]["ac"]
        algorithm = net["_options"]["algorithm"]

        net["converged"] = False
        net["OPF_converged"] = False
        _add_auxiliary_elements(net)

        if not ac or init_results:
            verify_results(net)
        else:
            reset_results(net)

        # TODO remove this when zip loads are integrated for all PF algorithms
        if algorithm not in ['nr', 'bfsw']:
            net["_options"]["voltage_depend_loads"] = False


        _add_auxiliary_elements(net)
        # convert pandapower net to ppc
        ppc, self.ppci = _pd2ppc(net)

        # store variables
        net["_ppc"] = ppc

        if not "VERBOSE" in kwargs:
            kwargs["VERBOSE"] = 0

        # ----- run the powerflow -----
        options = net["_options"]
        et_early_init = time() - t0_early_init

        # ---------- pp.powerflow._run_pf_algorithm() ----------------
        # ---------- pp.pf.run_newton_raphson_pf.run_newton_raphson_pf() ----------------
        t0 = time()
        t0_init = t0
        et_init_dc = 0.
        if need_reset:
            if isinstance(options["init_va_degree"], str) and options["init_va_degree"] == "dc":
                self.ppci = _run_dc_pf(self.ppci)
                et_init_dc = time() - t0
            if options["enforce_q_lims"]:
                raise NotImplementedError("enforce_q_lims not yet implemented")

            t0_init = time()
        # ---------- pp.pf.run_newton_raphson_pf._run_ac_pf_without_qlims_enforced ----------
        # ppci, success, iterations = _run_ac_pf_without_qlims_enforced(ppci, options)
            makeYbus, pfsoln = _get_numba_functions(self.ppci, options)
            self.baseMVA, self.bus, self.gen, self.branch, self.ref, self.pv, self.pq, _, _, V0, self.ref_gens = _get_pf_variables_from_ppci(self.ppci)
            self.ppci, self.Ybus, self.Yf, self.Yt = _get_Y_bus(self.ppci, options, makeYbus, self.baseMVA, self.bus, self.branch)


            # TODO i have a problem here for the order of the bus / id of bus
            tmp_bus_ind = np.argsort(net.bus.index)
            model = DataModel()
            model.set_sn_mva(net.sn_mva)
            model.set_f_hz(net.f_hz)

            # init_but should be called first among all the rest
            model.init_bus(net.bus.iloc[tmp_bus_ind]["vn_kv"].values, net.line.shape[0], net.line.shape[0])

            # init the shunts
            model.init_powerlines(net.line["r_ohm_per_km"].values * net.line["length_km"].values,
                                  net.line["x_ohm_per_km"].values * net.line["length_km"].values,
                                  net.line["c_nf_per_km"].values * net.line["length_km"].values,
                                  # net.line["g_us_per_km"].values * net.line["length_km"].values,
                                  net.line["from_bus"].values,
                                  net.line["to_bus"].values
                                  )

            # init the shunts
            model.init_shunt(net.shunt["p_mw"].values,
                             net.shunt["q_mvar"].values,
                             net.shunt["bus"].values
                             )
            # init trafo
            if net.trafo.shape[0]:
                trafo_r, trafo_x, trafo_b = model.get_trafo_param(net.trafo["vn_hv_kv"].values,
                                                                  net.trafo["vn_lv_kv"].values,
                                                                  net.trafo["vk_percent"].values,
                                                                  net.trafo["vkr_percent"].values,
                                                                  net.trafo["sn_mva"].values,
                                                                  net.trafo["pfe_kw"].values,
                                                                  net.trafo["i0_percent"].values,
                                                                  net.trafo["lv_bus"].values
                                                                  )

                # trafo_branch = ppc["branch"][net.line.shape[0]:, :]

                tap_step_pct = net.trafo["tap_step_percent"].values
                tap_step_pct[~np.isfinite(tap_step_pct)] = 0.

                tap_pos = net.trafo["tap_pos"].values
                tap_pos[~np.isfinite(tap_pos)] = 0.

                is_tap_hv_side = net.trafo["tap_side"].values == "hv"
                is_tap_hv_side[~np.isfinite(tap_pos)] = True
                model.init_trafo(trafo_r,
                                 trafo_x,
                                 trafo_b,
                                 tap_step_pct,
                                 tap_pos,
                                 is_tap_hv_side,
                                 net.trafo["hv_bus"].values,
                                 net.trafo["lv_bus"].values)

            model.init_loads(net.load["p_mw"].values,
                             net.load["q_mvar"].values,
                             net.load["bus"].values
                             )
            model.init_generators(net.gen["p_mw"].values,
                                  net.gen["vm_pu"].values,
                                  net.gen["bus"].values
                                  )
            # TODO better way here!
            model.add_slackbus(net.ext_grid["bus"].values)


            model.init_Ybus()
            Ybus = model.get_Ybus()

            # be careful, the order is not the same between this and pandapower, you need to change it
            # Ybus_proper_oder = Ybus[np.array([net.bus.index]).T, np.array([net.bus.index])]
            # self.Ybus_proper_oder = self.Ybus
            Ybus_proper_oder = Ybus
            self.Ybus_proper_oder = self.Ybus[np.array([tmp_bus_ind]).T, np.array([tmp_bus_ind])]
            tmp = np.abs(Ybus_proper_oder - self.Ybus_proper_oder)  # > 1e-7
        else:
            pass
            # TODO update self.ppci with new values of generation - load such that  Sbus is properly udpated


        # compute complex bus power injections [generation - load]
        Sbus = _get_Sbus(self.ppci, False)
        Sbus_me = model.get_Sbus()

        # pdb.set_trace()
        Sbus_me_r = np.real(Sbus_me)
        dcYbus = model.dc_pf(Sbus_me_r)
        pdb.set_trace()

        # run the newton power  flow
        # ------------------- pp.pypower.newtonpf ---------------------
        max_it = options["max_iteration"]
        tol = options['tolerance_mva']
        self.Ybus = sparse.csc_matrix(self.Ybus)
        et_init = time() - t0_init

        t0__ = time()
        if need_reset:
            # reset the solver
            self.solver.reset()
            self.V = V0
        else:
            # reuse previous voltages
            pass
        self.solver.solve(self.Ybus, self.V, Sbus, self.pv, self.pq, max_it, tol)
        et__ = time() - t0__

        t0_ = time()
        Va = self.solver.get_Va()
        Vm = self.solver.get_Vm()
        self.V = Vm * np.exp(1j * Va)
        J = self.solver.get_J()
        success = self.solver.converged()
        iterations = self.solver.get_nb_iter()
        # timer_Fx_, timer_solve_, timer_initialize_, timer_check_, timer_dSbus_, timer_fillJ_, timer_total_nr_
        timers = self.solver.get_timers()
        et_ = time() - t0_
        # ---------------------- pp.pypower.newtonpf ---------------------

        self.ppci = _store_internal(self.ppci, {"J": J, "Vm_it": None, "Va_it": None, "bus": self.bus, "gen": self.gen,
                                                "branch": self.branch,
                                                "baseMVA": self.baseMVA, "V": self.V, "pv": self.pv, "pq": self.pq,
                                                "ref": self.ref, "Sbus": Sbus,
                                                "ref_gens": self.ref_gens, "Ybus": self.Ybus, "Yf": self.Yf, "Yt": self.Yt,
                                                "timers": timers, "time_get_res": et_, "time_solve": et__,
                                                "time_init": et_init,
                                                "time_init_dc": et_init_dc,
                                                "time_early_init": et_early_init,
                                                "time_options": et_options}
                               )
        t0_ppci_to_pfsoln = time()
        # update data matrices with solution store in ppci
        # ---------- pp.pf.run_newton_raphson_pf._run_ac_pf_without_qlims_enforced ----------
        self.bus, self.gen, self.branch = ppci_to_pfsoln(self.ppci, options)
        te_ppci_to_pfsoln = time() - t0_ppci_to_pfsoln
            # these are the values from pypower / matlab
        t0_store_res = time()
        et = t0_store_res - t0
        result = _store_results_from_pf_in_ppci(self.ppci, self.bus, self.gen, self.branch, success, iterations, et)
        t0_to_net = time()
        et_store_res = t0_to_net - t0_store_res
        # ---------- pp.pf.run_newton_raphson_pf.run_newton_raphson_pf() ----------------
        # ---------- pp.powerflow._run_pf_algorithm() ----------------

        # read the results (=ppci with results) to net
        _ppci_to_net(result, net)
        et_to_net = time() - t0_to_net
        # ---------- pp.powerflow._powerflow() ----------------
        # ---------- pp.run.runpp() -----------------

        # added
        et_start = time() - t0_start
        self.ppci = _store_internal(self.ppci, {"time_store_res": et_store_res, "time_to_net": et_to_net,
                                                "time_all": et_start,
                                                "time_ppci_to_pfsoln": te_ppci_to_pfsoln})