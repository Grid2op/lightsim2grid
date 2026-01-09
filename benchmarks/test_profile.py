#!/home/donnotben/Documents/lightsim2grid/venv_ls/bin/python3

#!/bin/python3
import pickle
import time

from kiwisolver import Solver
import grid2op
from lightsim2grid import LightSimBackend
from lightsim2grid.solver import SolverType
import warnings
import pandapower as pp
import numpy as np

# usage:
# perf record ./test_profile.py
# perf report
# env_name = "l2rpn_neurips_2020_track2_small"

from  benchmark_grid_size import (
    get_loads_gens,
    make_grid2op_env
)
prng = np.random.default_rng(42)
    
CASE_NAME = "case9241pegase.json"
NB_TS = 10


def my_grid2op_env(case_name, nb_ts, prng):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        # env_tmp = grid2op.make(env_name)
        
        # load the case file
        case = pp.from_json(case_name)
        pp.runpp(case)  # for slack
        load_p_init = 1.0 * case.load["p_mw"].values
        load_q_init = 1.0 * case.load["q_mvar"].values
        gen_p_init = 1.0 * case.gen["p_mw"].values
        sgen_p_init = 1.0 * case.sgen["p_mw"].values
        load_p, load_q, gen_p, sgen_p = get_loads_gens(load_p_init, load_q_init, gen_p_init, sgen_p_init, prng)

        slack_gens =  np.zeros((nb_ts, case.ext_grid.shape[0]))
        if "res_ext_grid" in case:
            slack_gens += np.tile(case.res_ext_grid["p_mw"].values.reshape(1,-1), (nb_ts, 1))
        gen_p_g2op = np.concatenate((gen_p, slack_gens), axis=1)  
        env = make_grid2op_env(case,
                            case_name,
                            load_p,
                            load_q,
                            gen_p_g2op,
                            sgen_p)
            
        # param = env_tmp.parameters
        # param.NO_OVERFLOW_DISCONNECTION = True
        # env = grid2op.make(env_name, backend=LightSimBackend(), param=param)
        env.reset(seed=0, options={"time serie id": 0})
        env.backend._timer_preproc = 0
        env.backend._timer_solver = 0
        
        with open(f"gridmodel_{case_name}.pickle", "wb") as f:
            pickle.dump(obj=env.backend._grid, file=f)
        return env


def make_steps_glop(env, nb=NB_TS, reset_algo=False):
    for i in range(nb):
        if reset_algo:
            env.backend._grid.tell_solver_need_reset()
        _ = env.step(env.action_space())
    print(f"Total time: {env.backend._timer_preproc + env.backend._timer_solver}")


def main_glop():
    env = my_grid2op_env(CASE_NAME, NB_TS, prng)
    make_steps_glop(env, NB_TS, reset_algo=True)
    

def main_gridmodel(case_name=CASE_NAME, nb_ts=NB_TS, reset_algo=True, solver_used=SolverType.KLU):
    time_ls = 0.
    ls_timer_Fx = 0.
    ls_timer_solve = 0.
    ls_timer_initialize = 0.
    ls_timer_check = 0.
    ls_timer_dSbus = 0.
    ls_timer_fillJ = 0.
    ls_timer_Va_Vm = 0.
    ls_timer_pre_proc = 0.
    ls_timer_total_nr = 0.
    
    with open(f"gridmodel_{case_name}.pickle", "rb") as f:
        ls_grid = pickle.load(f)
    ls_grid.change_solver(solver_used)
    v_init = ls_grid.dc_pf(np.ones(ls_grid.get_bus_vn_kv().shape[0], dtype=complex) * 1.04, 1, 0.1)
    for _ in range(nb_ts):
        if reset_algo:
            ls_grid.tell_solver_need_reset()
        beg_ls = time.perf_counter()
        V = ls_grid.ac_pf(v_init, 10, 1e-6)
        end_ls = time.perf_counter()
        time_ls += end_ls - beg_ls
        if V.shape[0] == 0:
            raise RuntimeError("Divergence")
        
        (timer_Fx_, timer_solve_, timer_initialize_, 
         timer_check_, timer_dSbus_, timer_fillJ_, 
         timer_Va_Vm_, timer_pre_proc_, timer_total_nr_
         ) = ls_grid.get_solver().get_timers_jacobian()
        ls_grid.unset_changes()  # tell lightsim2grid that the state of the grid is consistent
        ls_timer_Fx += timer_Fx_
        ls_timer_solve += timer_solve_
        ls_timer_initialize += timer_initialize_
        ls_timer_check += timer_check_
        ls_timer_dSbus += timer_dSbus_
        ls_timer_fillJ += timer_fillJ_
        ls_timer_Va_Vm += timer_Va_Vm_
        ls_timer_pre_proc += timer_pre_proc_
        ls_timer_total_nr += timer_total_nr_
    print(f"Solver used: {solver_used}")
    print(f"Nb iter: {nb_ts}")
    print(f"Do reset each step: {reset_algo}")
    print("--------------------------------------")
    print("Detailed lightsim2grid timings: ")
    print(f"Total time {time_ls * 1000.:.2e}ms => {time_ls / nb_ts * 1e3:.2e} ms/pf | {nb_ts / time_ls:.0f} pf /s")
    print(f"Total time spent in the solver: {1e3 * ls_timer_total_nr:.2e} ms ({100. * ls_timer_total_nr / time_ls:.0f}% of total time spent in lightsim2grid)")
    print(f"\t Time to pre process Ybus, Sbus etc.: {1e3 * ls_timer_pre_proc:.2e} ms ({100. * ls_timer_pre_proc / ls_timer_total_nr:.0f} % of time in solver)")
    print(f"\t Time to initialize linear solver {1e3 * ls_timer_initialize:.2e} ms ({100. * ls_timer_initialize / ls_timer_total_nr:.0f} % of time in solver)")
    print(f"\t Time to compute dS/dV {1e3 * ls_timer_dSbus : .2e} ms ({100. * ls_timer_dSbus / ls_timer_total_nr:.0f} % of time in solver)")
    print(f"\t Time to fill the Jacobian {1e3 * ls_timer_fillJ:.2e} ms ({100. * ls_timer_fillJ / ls_timer_total_nr:.0f} % of time in solver)")
    print(f"\t Time to solve the Jacobian linear system: {1e3 * ls_timer_solve:.2e} ms ({100. * ls_timer_solve / ls_timer_total_nr:.0f} % of time in solver)")
    print(f"\t Time to update Va and Vm {1e3*ls_timer_Va_Vm:.2e} ms ({100. * ls_timer_Va_Vm / ls_timer_total_nr:.0f} % of time in solver)")
    print(f"\t Time to evaluate p,q mismmatch at each bus {1e3*ls_timer_Fx:.2e} ms ({100. * ls_timer_Fx / ls_timer_total_nr:.0f} % of time in solver)")
    print(f"\t Time to evaluate cvg criteria {1e3*ls_timer_check:.2e} ms ({100. * ls_timer_check / ls_timer_total_nr:.0f} % of time in solver)")
    print("--------------------------------------\n")
    
    
if __name__ == "__main__":
    # my_grid2op_env(case_name, nb_ts, prng)
    main_gridmodel(CASE_NAME, NB_TS, reset_algo=True, solver_used=SolverType.KLUSingleSlack)
