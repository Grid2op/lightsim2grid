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

from pyklu_cpp import KLUSolver

import pdb

# TODO optim: when i reuse the powerflow, i just update Sbus and reuse V and all the other stuff
# especially, i don't re create a solver etc.
# TODO just i want to test

def runpp(net, max_iteration=10, **kwargs):
    # ---------- pp.run.runpp() -----------------
    passed_parameters = _passed_runpp_parameters(locals())
    _init_runpp_options(net, algorithm="nr", calculate_voltage_angles="auto",
                        init="auto", max_iteration=max_iteration, tolerance_mva=1e-8,
                        trafo_model="t", trafo_loading="current",
                        enforce_q_lims=False, check_connectivity=True,
                        voltage_depend_loads=True,
                        consider_line_temperature=False,
                        passed_parameters=passed_parameters, numba=True, **kwargs)
    _check_bus_index_and_print_warning_if_high(net)
    _check_gen_index_and_print_warning_if_high(net)

    # ---------- pp.powerflow._powerflow() -----------------
    """
    Gets called by runpp or rundcpp with different arguments.
    """
    # get infos from options
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
    ppc, ppci = _pd2ppc(net)

    # store variables
    net["_ppc"] = ppc

    if not "VERBOSE" in kwargs:
        kwargs["VERBOSE"] = 0

    # ----- run the powerflow -----
    options = net["_options"]
    # ---------- pp.powerflow._run_pf_algorithm() ----------------
    # ---------- pp.pf.run_newton_raphson_pf.run_newton_raphson_pf() ----------------
    t0 = time()
    if isinstance(options["init_va_degree"], str) and options["init_va_degree"] == "dc":
        ppci = _run_dc_pf(ppci)
    if options["enforce_q_lims"]:
        raise NotImplementedError("enforce_q_lims not yet implemented")

    # ---------- pp.pf.run_newton_raphson_pf._run_ac_pf_without_qlims_enforced ----------
    # ppci, success, iterations = _run_ac_pf_without_qlims_enforced(ppci, options)
    makeYbus, pfsoln = _get_numba_functions(ppci, options)
    baseMVA, bus, gen, branch, ref, pv, pq, _, _, V0, ref_gens = _get_pf_variables_from_ppci(ppci)

    ppci, Ybus, Yf, Yt = _get_Y_bus(ppci, options, makeYbus, baseMVA, bus, branch)

    # compute complex bus power injections [generation - load]
    Sbus = _get_Sbus(ppci, False)

    # run the newton power  flow
    # ------------------- pp.pypower.newtonpf ---------------------
    max_it = options["max_iteration"]
    tol = options['tolerance_mva']
    Ybus = sparse.csc_matrix(Ybus)

    # t0 = time()
    solver = KLUSolver()
    solver.solve(Ybus, V0, Sbus, pv, pq, max_it, tol)
    # et = time() - t0

    t0_ = time()
    Va = solver.get_Va()
    Vm = solver.get_Vm()
    V = Vm * np.exp(1j * Va)
    J = solver.get_J()
    success = solver.converged()
    iterations = solver.get_nb_iter()
    et_ = time() - t0_
    # timer_Fx_, timer_solve_, timer_initialize_, timer_check_, timer_dSbus_, timer_fillJ_, timer_total_nr_
    timers = solver.get_timers()
    # ---------------------- pp.pypower.newtonpf ---------------------

    ppci = _store_internal(ppci, {"J": J, "Vm_it": None, "Va_it": None, "bus": bus, "gen": gen, "branch": branch,
                                  "baseMVA": baseMVA, "V": V, "pv": pv, "pq": pq, "ref": ref, "Sbus": Sbus,
                                  "ref_gens": ref_gens, "Ybus": Ybus, "Yf": Yf, "Yt": Yt,
                                  "timers": timers, "time_get_res": et_})

    # update data matrices with solution store in ppci
    # ---------- pp.pf.run_newton_raphson_pf._run_ac_pf_without_qlims_enforced ----------
    bus, gen, branch = ppci_to_pfsoln(ppci, options)
    et = time() - t0
    result = _store_results_from_pf_in_ppci(ppci, bus, gen, branch, success, iterations, et)
    # ---------- pp.pf.run_newton_raphson_pf.run_newton_raphson_pf() ----------------
    # ---------- pp.powerflow._run_pf_algorithm() ----------------

    # read the results (=ppci with results) to net
    _ppci_to_net(result, net)
    # ---------- pp.powerflow._powerflow() ----------------
    # ---------- pp.run.runpp() -----------------