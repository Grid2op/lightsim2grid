import numpy as np
import time

import pandapower as pp
import pandapower.networks as pn
from pyklu.compute_powerflow import runpp
import pdb
tol = 1e-5

grid1 = pn.case118()
grid2 = pn.case118()

grid1 = pn.case6515rte()
grid2 = pn.case6515rte()
nb_iteration = 100  # number of powerflow run
nb_max_newton_it = 10  # maximum number of iteration for the solver

### code begins here
powerflow_time_pp = 0
powerflow_time_cpp = 0
sucees_cpp = False
sucees_pp = False
timers_cpp = np.zeros(7)

# remove the first call with numba, that compiles the code
pp.runpp(grid1, max_iteration=nb_max_newton_it, numba=True)
start_time_pp = time.time()
for i in range(nb_iteration):
    pp.runpp(grid1, max_iteration=nb_max_newton_it, numba=True)
    powerflow_time_pp += grid1._ppc['et']
    sucees_pp = grid1._ppc['success']
end_time_pp = time.time()

start_time_cpp = time.time()
for i in range(nb_iteration):
    runpp(grid2, max_iteration=nb_max_newton_it)
    powerflow_time_cpp += grid2._ppc['et']
    sucees_cpp = grid2._ppc['success']
    timers_cpp += np.array(grid2._ppc["internal"]['timers'])

end_time_cpp = time.time()
print("Is there a solution:\n\tpandapower: {}\n\tcustom implementation: {}".format(sucees_pp, sucees_pp))
print("Total time:\n\tpandapower: {:.3f}s\n\tc++ custom implementation: {:.3f}s [results for {} iteration(s)]".format(
    end_time_pp-start_time_pp, end_time_cpp-start_time_cpp, nb_iteration))

print("Powerflow time: \n\tpandapower: {:.3f}s\n\tc++ custom implementation: {:.3f}s [results for {} iteration(s)]".format(
    powerflow_time_pp, powerflow_time_cpp, nb_iteration))

print("Detail timers for c++ implementations are:")
for time_, nm_var in zip(timers_cpp,
                         ["timer_Fx_", "timer_solve_", "timer_initialize_", "timer_check_",
                          "timer_dSbus_", "timer_fillJ_", "timer_total_nr_"]):
    print("\t {}: {:.4f}s".format(nm_var, time_))

# check the results are the same for buses
for colname in ["vm_pu", "va_degree", "p_mw", "q_mvar"]:
        if np.sum(np.abs(grid1.res_bus[colname] - grid2.res_bus[colname])) > tol:
            print("mismatch for \"{}\" in res_bus".format(colname))
# check the results are the same for generators
for colname in ["vm_pu", "va_degree", "p_mw", "q_mvar"]:
        if np.sum(np.abs(grid1.res_gen[colname] - grid2.res_gen[colname])) > tol:
            print("mismatch for \"{}\" in res_gen".format(colname))

# check the results are the same for generators
for colname in ["p_mw", "q_mvar"]:
        if np.sum(np.abs(grid1.res_load[colname] - grid2.res_load[colname])) > tol:
            print("mismatch for \"{}\" in res_load".format(colname))

# same for line
for colname in ["p_from_mw", "q_from_mvar", "p_to_mw", "q_to_mvar",
                "vm_from_pu", "va_from_degree", "vm_to_pu", "va_to_degree"]:
        if np.sum(np.abs(grid1.res_line[colname] - grid2.res_line[colname])) > tol:
            print("mismatch for \"{}\" in res_line".format(colname))

# and trafo
for colname in ["p_hv_mw", "q_hv_mvar", "p_lv_mw", "q_lv_mvar",
                "vm_hv_pu", "va_hv_degree", "vm_lv_pu", "va_lv_degree"]:
        if np.sum(np.abs(grid1.res_trafo[colname] - grid2.res_trafo[colname])) > tol:
            print("mismatch for \"{}\" in res_line".format(colname))


