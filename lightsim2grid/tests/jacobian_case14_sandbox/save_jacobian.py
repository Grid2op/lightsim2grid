import pickle

import grid2op
import numpy as np

from lightsim2grid import LightSimBackend
try:
    from lightsim2grid.algorithm import NRSing_KLU
except ImportError:
    # lightsim2grid version < 0.14
    from lightsim2grid.solver import KLUSolver as NRSing_KLU
    
env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
obs = env.reset(seed=0, options={"time serie id": 0})

n_iter = 0
Ybus = env.backend._grid.get_Ybus_solver().copy()
res = {
"init_state":
    {
        "tol": 1e-8,
        "Sbus": env.backend._grid.get_Sbus_solver().copy(),
        "Ybus": Ybus,
        "v_init": np.ones(Ybus.shape[0], dtype=complex),
        "pv": env.backend._grid.get_pv_solver().copy(),
        "pq":  env.backend._grid.get_pq_solver().copy(),
        "slack_weights": env.backend._grid.get_slack_weights_solver().copy(),
        "slack_ids": env.backend._grid.get_slack_ids_solver().copy(),
    }
}
for n_iter in range(5):
    solver = NRSing_KLU()
    solver.compute_pf(
        res["init_state"]["Ybus"],
        res["init_state"]["v_init"],
        res["init_state"]["Sbus"],
        res["init_state"]["slack_ids"],
        res["init_state"]["slack_weights"],
        res["init_state"]["pv"],
        res["init_state"]["pq"],
        n_iter,
        res["init_state"]["tol"])
    J = solver.get_J()
    V = solver.get_V()
    res[f"{n_iter}"] = {
        "J": J.copy(),
        "V": V.copy()
    }
assert solver.converged()

with open("saved_jacobian_single.pkl", "wb") as f:
    pickle.dump(res, f)
