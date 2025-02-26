import time
from lightsim2grid_cpp import LightEnv
from lightsim2grid_cpp import Protections
import numpy as np
import grid2op

from lightsim2grid import LightSimBackend
# env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
env = grid2op.make("l2rpn_idf_2023", test=True, backend=LightSimBackend())

light_env = LightEnv(env.backend._grid)

# raises error due to lack of load_p
# obs = light_env.reset()

nb_ts = env.chronics_handler.real_data.data.load_p.shape[0]
light_env.assign_time_series(
    env.chronics_handler.real_data.data.load_p,
    env.chronics_handler.real_data.data.load_q,
    env.chronics_handler.real_data.data.prod_p,
    np.full((nb_ts, env.n_gen), fill_value=np.nan),  # gen_v
    np.full((nb_ts, env.n_storage), fill_value=np.nan),  # storage p
    np.full((nb_ts, env.n_shunt), fill_value=np.nan),  # shunt_p
    np.full((nb_ts, env.n_shunt), fill_value=np.nan),  # shunt_q
    np.full((nb_ts, 0), fill_value=np.nan),  # static gen p 
    np.full((nb_ts, 0), fill_value=np.nan),  # static gen q
)

protections = Protections()
protections.set_thermal_limit_or(10. * env.get_thermal_limit().astype(dtype=float))
protections.set_thermal_limit_ex(env.get_thermal_limit().astype(dtype=float) / 
                                 env.backend.lines_ex_pu_to_kv *
                                 env.backend.lines_or_pu_to_kv * 
                                 10.)
protections.set_max_line_time_step_overflow(np.array([3 for _ in range(env.n_line)], dtype=np.int32))
light_env.protections = protections
print("starting now")
obs = light_env.reset()
done = False
while not done:
    obs, reward, done, truncated, info = light_env.step(0)
# import pdb
# pdb.set_trace()
# light_env.protections.line_disconnected
total_time_ms = 1e3*(light_env.step_time + light_env.reset_time)
nb_stuff = 1 + light_env.current_step # 1+ to take into account the env.reset call
step_per_sec = int(1e3 * nb_stuff / total_time_ms / 100) * 100
pow_per_sec = int(nb_stuff / light_env.protections.powerflow_time / 100) * 100
print(f"Total computation time: {total_time_ms:.2f}ms for {nb_stuff} steps "
      f"({total_time_ms/nb_stuff:.2e} ms / step, {step_per_sec:2.0f} step / s)")
print("Among which: ")
print(f"\t - {1e3 * light_env.reset_time:.2f}ms spent for the initial reset")
print(f"\t - {1e3 * light_env.step_time:.2f}ms spent for the steps")
print(f"\t - {1e3 * light_env.obs_time:.2f}ms spent to retrieve the observation from the protections")
print(f"\t - {1e3 * light_env.update_gridmodel_time:.2f}ms spent to update the gridmodel (new injections)")
print(f"\t - {1e3 * light_env.protections.total_time:.2f}ms spent in the protections")
print(f"\t - {1e3 * light_env.protections.powerflow_time:.2f}ms spent in the protections to perform the powerflows (DC + AC)")
print(f"\t   this is {1e3 * light_env.protections.powerflow_time / nb_stuff:.3f}ms/powerflow or {pow_per_sec} powerflow/s")
