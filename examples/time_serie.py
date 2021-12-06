#!/bin/python3

# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

# sudo sh -c 'echo 1 >/proc/sys/kernel/perf_event_paranoid'
# perf record ./time_serie.py
# 
import time
import numpy as np

import grid2op
from grid2op.Parameters import Parameters
from grid2op.Chronics import GridStateFromFileWithForecastsWithoutMaintenance
import warnings
from lightsim2grid import LightSimBackend, TimeSerie

env_name = "l2rpn_neurips_2020_track2"  # use "l2rpn_neurips_2020_track2_small" if `test=True`
test = True
scenario_id = 1

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
                                 test=test,
                                 # ignore the maintenance, that are NOT simulated
                                 # by the TimeSerie module !
                                 data_feeding_kwargs={"gridvalueClass": GridStateFromFileWithForecastsWithoutMaintenance}
                                 )
key_env = max([el for el in multi_mix_env.keys()])
env = multi_mix_env[key_env]

# Run the environment on a scenario using the TimeSerie module
time_serie = TimeSerie(env)
time_serie.compute_V(scenario_id=scenario_id)
a_or = time_serie.compute_A()
computer = time_serie.computer

print(f"For environment: {env_name}")
print(f"Total time spent in \"computer\" to solve everything: {computer.total_time():.2f}s "
      f"({computer.nb_solved() / computer.total_time():.0f} pf / s), "
      f"{1000.*computer.total_time() / computer.nb_solved():.2f} ms / pf)")
print(f"\t - time to pre process the injections: {computer.preprocessing_time():.2f}s")
print(f"\t - time to perform powerflows: {computer.solver_time():.2f}s "
      f"({computer.nb_solved() / computer.solver_time():.0f} pf / s, "
      f"{1000.*computer.solver_time() / computer.nb_solved():.2f} ms / pf)")
print(f"In addition, it took {computer.amps_computation_time():.2f} s to retrieve the current "
      f"from the complex voltages (in total "
      f"{computer.nb_solved() / ( computer.total_time() + computer.amps_computation_time()):.1f} "
      "pf /s, "
      f"{1000.*( computer.total_time() + computer.amps_computation_time()) / computer.nb_solved():.2f} ms / pf)")

#### Comparison with running grid2op
env.set_id(scenario_id)
_ = env.reset()
beg_ = time.perf_counter()
done = False
while not done:
      *_, done, info = env.step(env.action_space())
end_ = time.perf_counter()
print("Comparison with raw grid2op timings")
print(f"It took grid2op: {end_ - beg_:.2f}s to perform the same computation")
print(f"This is a {(end_ - beg_) / ( computer.total_time() + computer.amps_computation_time()) :.1f} speed up from TimeSerie over raw grid2op")

#### Check that the results matches
env.set_id(scenario_id)
obs = env.reset()
assert np.max(np.abs(obs.a_or - a_or[0])) <= 1e-4, "error for first step (just after reset)"
ts = 1
done = False
while not done:
      obs, reward, done, info = env.step(env.action_space())
      assert np.max(np.abs(obs.a_or - a_or[ts])) <= 1e-4, f"wrong amps at step {ts}"
      ts += 1