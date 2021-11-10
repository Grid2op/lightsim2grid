# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import time
import numpy as np

import grid2op
from grid2op.Parameters import Parameters
from grid2op.Action import BaseAction
from grid2op.Chronics import ChangeNothing
import warnings
from lightsim2grid import LightSimBackend, SecurityAnalysis

env_name = "l2rpn_neurips_2020_track2_small"
test = False

# Create the grid2op environment
param = Parameters()
param.NO_OVERFLOW_DISCONNECTION = True
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    multi_mix_env = grid2op.make(env_name,
                                 backend=LightSimBackend(),
                                 # ignore the protection, that are NOT simulated
                                 # by the TimeSerie module !
                                 param=param,
                                 test=test
                                 )

key_env = max([el for el in multi_mix_env.keys()])
env = multi_mix_env[key_env]

# Run the environment on a scenario using the TimeSerie module
security_analysis = SecurityAnalysis(env)
security_analysis.add_all_n1_contingencies()
a_or, voltages = security_analysis.get_flows()
# the 3 lines above are the only lines you need to do to perform a security analysis !

computer = security_analysis.computer
print(f"For environment: {env_name} ({computer.nb_solved()} n-1 simulated)")
print(f"Total time spent in \"computer\" to solve everything: {1e3*computer.total_time():.1f}ms "
      f"({computer.nb_solved() / computer.total_time():.0f} pf / s), "
      f"{1000.*computer.total_time() / computer.nb_solved():.2f} ms / pf)")
print(f"\t - time to compute the coefficients to simulate line disconnection: {1e3*computer.preprocessing_time():.2f}ms")
print(f"\t - time to pre process Ybus: {1e3*computer.modif_Ybus_time():.2f}ms")
print(f"\t - time to perform powerflows: {1e3*computer.solver_time():.2f}ms "
      f"({computer.nb_solved() / computer.solver_time():.0f} pf / s, "
      f"{1000.*computer.solver_time() / computer.nb_solved():.2f} ms / pf)")
print(f"In addition, it took {1e3*computer.amps_computation_time():.2f} ms to retrieve the current "
      f"from the complex voltages (in total "
      f"{computer.nb_solved() / ( computer.total_time() + computer.amps_computation_time()):.1f} "
      "pf /s, "
      f"{1000.*( computer.total_time() + computer.amps_computation_time()) / computer.nb_solved():.2f} ms / pf)")

#### Comparison with running grid2op
obs = env.get_obs()
beg_ = time.perf_counter()
for cont_id in range(env.n_line):
    sim_obs, sim_r, sim_d, sim_info = obs.simulate(env.action_space({"set_line_status": [(cont_id, -1)]}),
                                                   time_step=0)
end_ = time.perf_counter()
total_time = end_ - beg_
print("Comparison with raw grid2op timings")
print(f"It took grid2op: {end_ - beg_:.2f}s to perform the same computation")
print(f"This is a {(end_ - beg_) / ( computer.total_time() + computer.amps_computation_time()) :.1f} "
      f"speed up from SecurityAnalysis over raw grid2op (using obs.simulate)")


#### Check that the results match
for cont_id in range(env.n_line):
    action = env.action_space({"set_line_status": [(cont_id, -1)]})
    sim_obs, sim_r, sim_d, sim_info = obs.simulate(action,
                                                   time_step=0)
    if not sim_d:
        if np.max(np.abs(sim_obs.a_or - a_or[cont_id,:])) > 1e-4:
            raise RuntimeError(f"wrong amps for contingency {cont_id}")
    else:
        assert np.all(~np.isfinite(a_or[cont_id,:])), f"amps should be Nan for cont {cont_id}"
