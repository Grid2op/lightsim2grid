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
import warnings
from lightsim2grid import LightSimBackend, ContingencyAnalysis

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
                                 # by the ContingencyAnalysis module !
                                 param=param,
                                 test=test
                                 )
    multi_mix_env_pp = grid2op.make(env_name,
                                    # ignore the protection, that are NOT simulated
                                    # by the ContingencyAnalysis module !
                                    param=param,
                                    test=test
                                    )

key_env = max([el for el in multi_mix_env.keys()])
env = multi_mix_env[key_env]
env_pp = multi_mix_env_pp[key_env]

# Run the environment on a scenario using the ContingencyAnalysis module
contingency_analysis = ContingencyAnalysis(env)
contingency_analysis.add_all_n1_contingencies()
# (optional, usually faster) warm-start each contingency's powerflow with the "n" case
# voltage solution instead of a flat start
contingency_analysis.init_from_n_powerflow = True
p_or, a_or, voltages = contingency_analysis.get_flows()
# the 3 lines above (ignoring the optional one) are all you need to perform a contingency analysis !

computer = contingency_analysis.computer
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
total_time_glop_ls = end_ - beg_

obs_pp = env_pp.get_obs()
beg_ = time.perf_counter()
for cont_id in range(env.n_line):
    sim_obs, sim_r, sim_d, sim_info = obs_pp.simulate(env.action_space({"set_line_status": [(cont_id, -1)]}),
                                                      time_step=0)
end_ = time.perf_counter()
total_time_glop_pp = end_ - beg_

full_time_sa = computer.total_time() + computer.amps_computation_time()
print()
print("Comparison with raw grid2op timings")
print(f"It took grid2op (with lightsim2grid, using obs.simulate): {total_time_glop_ls:.2f}s to perform the same computation")
print(f"\t This is a {(total_time_glop_ls) / (full_time_sa) :.1f} "
      f"speed up from ContingencyAnalysis over raw grid2op (using obs.simulate and lightsim2grid)")
print(f"It took grid2op (with pandapower, using obs.simulate): {total_time_glop_pp:.2f}s to perform the same computation")
print(f"\t This is a {(total_time_glop_pp) / (full_time_sa) :.1f} "
      f"speed up from ContingencyAnalysis over raw grid2op (using obs.simulate and pandapower)")


#### Check that the results match
for cont_id in range(env.n_line):
    action = env.action_space({"set_line_status": [(cont_id, -1)]})
    sim_obs, sim_r, sim_d, sim_info = obs.simulate(action,
                                                   time_step=0)
    if not sim_d:
        # simulation converged
        assert np.all(np.isfinite(a_or[cont_id,:])), f"amps should not be Nan for cont {cont_id} (converged)"
        if np.max(np.abs(sim_obs.a_or - a_or[cont_id,:])) > 1e-4:
            raise RuntimeError(f"wrong amps for contingency {cont_id}")
        if np.max(np.abs(sim_obs.p_or - p_or[cont_id,:])) > 1e-4:
            raise RuntimeError(f"wrong active power for contingency {cont_id}")
    else:
        # simulation diverged
        assert np.all(~np.isfinite(a_or[cont_id,:])), f"amps should be Nan for cont {cont_id} (diverged)"
        assert np.all(np.abs(p_or[cont_id,:]) <= 1e-6), f"active power should be 0. for cont {cont_id}"

print("All results match !")


#### Advanced usage: multi-threading, disconnected grids and limit violations
# `compute_limit_violations` must be set (either here or at construction time) to use
# `run` / `run_ac` / `run_dc`; it makes the C++ side check every bus voltage and
# branch/trafo current against the limits set on the grid, per contingency. Changing it
# clears any previously registered contingency, so the n-1 list must be re-added after.
contingency_analysis.compute_limit_violations = True
# split the contingencies across 2 OS threads instead of solving them sequentially;
# results do not depend on `nb_thread`
contingency_analysis.nb_thread = 2
# a contingency that splits the grid into several components would otherwise make that
# component's powerflow diverge; with this set to True, only the elements that end up
# disconnected from the slack are reported as "not simulated" instead
contingency_analysis.handle_disconnected_grid = True
contingency_analysis.add_all_n1_contingencies()  # re-add: compute_limit_violations cleared them above

# `run_ac` (or `run_dc`) also makes sure the corresponding algorithm family is selected
result = contingency_analysis.run_ac()
print(f"\nPre-contingency ('n' case) converged: {result.pre_contingency_result.converged}, "
      f"{len(result.pre_contingency_result.limit_violations)} limit violation(s)")
nb_diverged = sum(not cont.converged for cont in result.post_contingency_results)
print(f"Out of {len(result.post_contingency_results)} contingencies, {nb_diverged} did not converge")
