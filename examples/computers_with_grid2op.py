# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

# ADVANCED USAGE
# This files explains how to use the Computers cpp class, for easier use
# please consult the documentation of TimeSeries or the
# time_serie.py file !

import grid2op
from grid2op.Parameters import Parameters
import warnings
import numpy as np
from lightsim2grid import LightSimBackend
from lightsim2grid_cpp import Computers

env_name = "l2rpn_neurips_2020_track2_small"
test = False
param = Parameters()
param.NO_OVERFLOW_DISCONNECTION = True
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    env = grid2op.make(env_name, backend=LightSimBackend(), param=param, test=test)

nb_bus = env.n_sub
obs = env.reset()
grid = env.backend._grid
Vinit = env.backend.V
prod_p = 1.0 * env.chronics_handler.real_data.data.prod_p
load_p = 1.0 * env.chronics_handler.real_data.data.load_p
load_q = 1.0 * env.chronics_handler.real_data.data.load_q
nb_sim = prod_p.shape[0]

# now perform the computation
computer = Computers(grid)
# print("start the computation")
status = computer.compute_Vs(prod_p,
                            np.zeros((nb_sim, 0)),  # no static generators for now !
                            load_p,
                            load_q,
                            Vinit,
                            env.backend.max_it,
                            env.backend.tol)
if status != 1:
    raise RuntimeError(f"Some error occurred, the powerflow has diverged after {computer.nb_solved()} step(s)")

assert computer.nb_solved() == nb_sim, f"Error should have made {nb_sim} powerflows, but did {computer.nb_solved()}"
Vs = computer.get_voltages()

# check i can call the method to get the buses
sbuses = computer.get_sbuses()
ampss = computer.compute_flows()

# results should be different
assert np.any(np.abs(np.min(np.abs(Vs), axis=0) - np.max(np.abs(Vs), axis=0)) > 0.01)

# I got the same voltages as a normal pf
for it_num in range(100):
    obs, *_ = env.step(env.action_space())
    if np.max(np.abs(env.backend.V[:nb_bus] - Vs[1 + it_num, :nb_bus]))  > 1e-6:
        raise RuntimeError(f"error at it {it_num} for voltages")
    if np.max(np.abs(ampss[1 + it_num] - obs.a_or / 1000.)) > 1e-6:
        raise RuntimeError(f"error at it {it_num} for amps")

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
