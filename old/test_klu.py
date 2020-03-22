import numpy as np
import time

import pandapower as pp
import pandapower.networks as pn
from pyklu.compute_powerflow import KLU4Pandapower

import pdb
tol = 1e-6

grid1 = pn.case118()
grid2 = pn.case118()

# grid1 = pn.case39()
# grid2 = pn.case39()

# grid1 = pn.case6515rte()
# grid2 = pn.case6515rte()

# grid1 = pn.case9241pegase()
# grid2 = pn.case9241pegase()
# grid1 = pn.case1888rte()
# grid2 = pn.case1888rte()

# without trafo, ordering issue nor shunt, no parrallel line
grid1 = pn.case5()
grid2 = pn.case5()

# without trafo, ordering issue, nor shunt, no parrallel line
grid1 = pn.case6ww()
grid2 = pn.case6ww()

# without trafo, but ordering issue and shunt, no parrallel line
grid1 = pn.case30()
grid2 = pn.case30()

# with trafo, ordering issue, and shunt, no parrallel line
grid1 = pn.case14()
grid2 = pn.case14()


nb_iteration = 1  # number of powerflow run
nb_max_newton_it = 10  # maximum number of iteration for the solver

### code begins here
powerflow_time_pp = 0
powerflow_time_cpp = 0
sucees_cpp = False
sucees_pp = False
timers_cpp = np.zeros(7)
nb_it_pp = -1
nb_it_cpp = -1
total_time_cpp = 0.
total_time_pp = 0.

time_get_res = 0.  # time to convert the results from c++ to python
time_init = 0.  # time to initialize Ybus from ppc
time_solve = 0.  # time to do the newton raphson
time_store_res = 0.  # time to do the newton raphson
time_to_net = 0.  # time to do the newton raphson
time_init_dc = 0.  # time to do the initialization
time_early_init = 0.  # time to do the initialization of the parameters
time_options = 0.  # time to parse the options
time_all = 0.  # time to do all computations
time_ppci_to_pfsoln = 0.  # time to parse pcci back to proper shape

cpp_solver = KLU4Pandapower()
# remove the first call with numba, that compiles the code
pp.runpp(grid1, max_iteration=nb_max_newton_it, numba=True)
start_time_pp = time.time()
for i in range(nb_iteration):
    start_time_pp = time.time()
    try:
        pp.runpp(grid1, max_iteration=nb_max_newton_it, numba=True)
    except pp.powerflow.LoadflowNotConverged:
        pass
    total_time_pp += time.time() - start_time_pp
    powerflow_time_pp += grid1._ppc['et']
    nb_it_pp = grid1._ppc['iterations']
    sucees_pp = grid1._ppc['success']
end_time_pp = time.time()

for i in range(nb_iteration):
    start_time_cpp = time.time()
    try:
        cpp_solver.runpp(grid2, max_iteration=nb_max_newton_it,
                         need_reset=True   # reset the KLU solver to an original state, need to be done each time the Ymatrix is changed (might be slow)
                         )
    except pp.powerflow.LoadflowNotConverged:
        if i == 0:
            print("I didn't converge for c++")
        pass
    total_time_cpp += time.time() - start_time_cpp

    powerflow_time_cpp += grid2._ppc['et']
    sucees_cpp = grid2._ppc['success']
    timers_cpp += np.array(grid2._ppc["internal"]['timers'])
    nb_it_cpp = grid2._ppc['iterations']
    time_get_res += grid2._ppc["internal"]['time_get_res']
    time_solve += grid2._ppc["internal"]['time_solve']
    time_init += grid2._ppc["internal"]['time_init']
    time_to_net += grid2._ppc["internal"]['time_to_net']
    time_store_res += grid2._ppc["internal"]['time_store_res']
    time_init_dc += grid2._ppc["internal"]['time_init_dc']
    time_early_init += grid2._ppc["internal"]['time_early_init']
    time_options += grid2._ppc["internal"]['time_options']
    time_all += grid2._ppc["internal"]['time_all']
    time_ppci_to_pfsoln += grid2._ppc["internal"]['time_ppci_to_pfsoln']


print("Is there a solution:\n\tpandapower: {}\n\tcustom implementation: {}".format(sucees_pp, sucees_pp))
print("Total time:\n\tpandapower: {:.3f}s\n\tc++ custom implementation: {:.3f}s [results for {} repeat(s)]".format(
    total_time_pp, total_time_cpp, nb_iteration))

print("Number of newton iteration (per cases): \n\tpandapower: {}\n\tc++ custom implementation: {}".format(
    nb_it_pp, nb_it_cpp))

print("Powerflow time: \n\tpandapower: {:.3f}s\n\tc++ custom implementation: {:.3f}s [results for {} repeat(s)]".format(
    powerflow_time_pp, powerflow_time_cpp, nb_iteration))

print("\nDetailled timers for c++ implementations are:")
for time_, nm_var in zip(timers_cpp,
                         ["timer_Fx_", "timer_solve_", "timer_initialize_", "timer_check_",
                          "timer_dSbus_", "timer_fillJ_", "timer_total_nr_"]):
    print("\t {}: {:.4f}s".format(nm_var, time_))
print("\nConvertion timers:")
print("\t time to configure the options: {:.3f}".format(time_options))
print("\t time to do the preliminary steps: {:.3f}".format(time_early_init))
print("\t time to do the powerflow: {:.3f}".format(powerflow_time_cpp))
print("\t\t time to do the dc approx to init: {:.3f}".format(time_init_dc))
print("\t\t time to init the newton raphson: {:.3f}".format(time_init))
print("\t\t time to solve the newton raphson (python side): {:.3f}".format(time_solve))
print("\t\t time to convert results c++ -> python: {:.3f}".format(time_get_res))
print("\t\t time to get back results (ppci_to_pfsoln): {:.3f}".format(time_ppci_to_pfsoln))
print("\t time to store results: {:.3f}".format(time_store_res))
print("\t time to convert to net: {:.3f}".format(time_to_net))
print("time for all KLUPP: {:.3f}".format(time_all))
# print(time_options + time_early_init + powerflow_time_cpp + time_store_res + time_to_net)
# check the results are the same for buses
print()
passed = True
for colname in ["vm_pu", "va_degree", "p_mw", "q_mvar"]:
        if np.sum(np.abs(grid1.res_bus[colname] - grid2.res_bus[colname])) > tol:
            print("mismatch for \"{}\" in res_bus".format(colname))
            passed = False
# check the results are the same for generators
for colname in ["vm_pu", "va_degree", "p_mw", "q_mvar"]:
        if np.sum(np.abs(grid1.res_gen[colname] - grid2.res_gen[colname])) > tol:
            print("mismatch for \"{}\" in res_gen".format(colname))
            passed = False

# check the results are the same for generators
for colname in ["p_mw", "q_mvar"]:
        if np.sum(np.abs(grid1.res_load[colname] - grid2.res_load[colname])) > tol:
            print("mismatch for \"{}\" in res_load".format(colname))
            passed = False

# same for line
for colname in ["p_from_mw", "q_from_mvar", "p_to_mw", "q_to_mvar",
                "vm_from_pu", "va_from_degree", "vm_to_pu", "va_to_degree"]:
        if np.sum(np.abs(grid1.res_line[colname] - grid2.res_line[colname])) > tol:
            print("mismatch for \"{}\" in res_line".format(colname))
            passed = False

# and trafo
for colname in ["p_hv_mw", "q_hv_mvar", "p_lv_mw", "q_lv_mvar",
                "vm_hv_pu", "va_hv_degree", "vm_lv_pu", "va_lv_degree"]:
        if np.sum(np.abs(grid1.res_trafo[colname] - grid2.res_trafo[colname])) > tol:
            print("mismatch for \"{}\" in res_line".format(colname))
            passed = False

if passed:
    print("All test succesfully passed!")
else:
    print("Some results are differents! FAIL")


