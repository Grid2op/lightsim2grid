import warnings
import pandapower as pp
import numpy as np        
from grid2op import make, Parameters
from grid2op.Chronics import FromNPY
from lightsim2grid import LightSimBackend
import tempfile
import os

try:
    from tabulate import tabulate
    TABULATE_AVAIL = True
except ImportError:
    print("The tabulate package is not installed. Some output might not work properly")
    TABULATE_AVAIL = False


case_names = [
            #   "case14.json",
              "case118.json",
            #   "case_illinois200.json",
            #   "case300.json",
            #   "case1354pegase.json",
              "case1888rte.json",
            #   "GBnetwork.json",  # 2224 buses
              "case2848rte.json",
            #   "case2869pegase.json",
            #   "case3120sp.json",
              "case6495rte.json",
              "case6515rte.json",
              "case9241pegase.json"
              ]

case_name = "case6495rte.json"
case_name = "case14.json"
    
def make_grid2op_env(pp_case, case_name, load_p, load_q, gen_p, sgen_p):
    param = Parameters.Parameters()
    param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})
        
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env_lightsim = make("blank",
                            param=param,
                            test=True,
                            backend=LightSimBackend(),
                            chronics_class=FromNPY,
                            data_feeding_kwargs={"load_p": load_p,
                                                 "load_q": load_q,
                                                 "prod_p": gen_p
                                                },
                            grid_path=case_name,
                            _add_to_name=f"{case_name}",
                            )
    return env_lightsim

if __name__ == "__main__":

    import pandapower.networks as pn
    for case_name in case_names:
        tmp_nm =  os.path.splitext(case_name)[0]
        print(f"====================== {tmp_nm} ======================")
        case = getattr(pn,tmp_nm)()
        pp.runpp(case)  # for slack
        
        load_p_init = 1.0 * case.load["p_mw"].values
        load_q_init = 1.0 * case.load["q_mvar"].values
        gen_p_init = 1.0 * case.gen["p_mw"].values
        sgen_p_init = 1.0 * case.sgen["p_mw"].values
        
        nb_ts = 1
        # add slack !
        slack_gens =  np.zeros((nb_ts, case.ext_grid.shape[0]))
        if "res_ext_grid" in case:
            slack_gens += np.tile(case.res_ext_grid["p_mw"].values.reshape(1,-1), (nb_ts, 1))
        gen_p_g2op = np.concatenate((gen_p_init.reshape(1,-1), slack_gens), axis=1)  
        
        with tempfile.TemporaryDirectory() as tmpdirname:
            pp.to_json(case, os.path.join(tmpdirname, case_name))
            with open(os.path.join(tmpdirname, "config.py"), "w") as f:
                f.write("config = {}")
            
            env = make_grid2op_env(None,
                                os.path.join(tmpdirname, case_name),
                                load_p=load_p_init.reshape(1,-1),
                                load_q=load_q_init.reshape(1,-1),
                                gen_p=gen_p_g2op.reshape(1,-1),
                                sgen_p=None)
            
        env.backend._grid.tell_solver_need_reset()
        _ = env.step(env.action_space())
        ls_solver =  env.backend._grid.get_solver()
        nb_iter_solver = ls_solver.get_nb_iter()
        timers = ls_solver.get_timers_jacobian()
        (timer_Fx, timer_solve, timer_init, timer_check, 
        timer_compute_dS, timer_fillJ, timer_compVa_Vm, timer_preproc, timer_total) = timers
        print(f"Total time for the powerflow (=pre proc + NR + post proc): {env.backend._grid.timer_last_ac_pf:.2e}s")
        print(f"Total time spent in the Newton Raphson: {timer_total:.2e}s")
        print(f"Time to pre process input data: {timer_preproc:.2e}s")
        print(f"Time to intialize linear solver: {timer_init:.2e}s")
        print(f"Then for all iterations (cumulated time over all {nb_iter_solver} iterations)")
        print(f"\ttotal time to compute dS/dVm and dS/dVa: {timer_compute_dS:.2e}s")
        print(f"\ttotal time fill jacobian matrix (from dS/dVm and dS/dVa): {timer_fillJ:.2e}s")
        print(f"\ttotal time to solve J.x = b: {timer_solve:.2e}s")
        print(f"\ttotal time to compute V, Va and Vm: {timer_compVa_Vm:.2e}s")
        print(f"\ttotal time to compute p, q mismatch at buses: {timer_Fx:.2e}s")
        print(f"\ttotal time to check if p,q mismatch at buses are within tolerance: {timer_check:.2e}s")
        print(f"====================== {' '*len(tmp_nm)} ======================")
        