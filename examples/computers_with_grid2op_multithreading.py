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

import time
import warnings
import numpy as np
import sys
import concurrent.futures  # thread
from multiprocessing import Pool  # multiprocessing
import asyncio  # asyncio

import grid2op
from grid2op.Parameters import Parameters
from lightsim2grid import LightSimBackend
from lightsim2grid.timeSerie import Computers

NB_THREAD = 4
ENV_NAME = "l2rpn_neurips_2020_track2_small"

def get_env(env_name):
    param = Parameters()
    param.NO_OVERFLOW_DISCONNECTION = True
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env = grid2op.make(env_name, backend=LightSimBackend(), param=param)
    return env

def get_injs(env):
    nb_bus = env.n_sub
    prods_p = []
    loads_p = []
    loads_q = []
    for _ in range(NB_THREAD):
        obs = env.reset()
        grid = env.backend._grid
        Vinit = env.backend.V
        prods_p.append(1.0 * env.chronics_handler.real_data.data.prod_p)
        loads_p.append(1.0 * env.chronics_handler.real_data.data.load_p)
        loads_q.append(1.0 * env.chronics_handler.real_data.data.load_q)
    return prods_p, loads_p, loads_q

def get_flows(grid, Vinit, prod_p, load_p, load_q, max_it=10, tol=1e-8):
    # now perform the computation
    computer = Computers(grid)
    # print("start the computation")
    status = computer.compute_Vs(prod_p,
                                 np.zeros((prod_p.shape[0], 0)),  # no static generators for now !
                                 load_p,
                                 load_q,
                                 Vinit,
                                 max_it,
                                 tol)
    if status != 1:
        raise RuntimeError(f"Some error occurred, the powerflow has diverged after {computer.nb_solved()} step(s)")
    nb_sim = prod_p.shape[0]

    assert computer.nb_solved() == nb_sim, f"Error should have made {nb_sim} powerflows, but did {computer.nb_solved()}"

    # check i can call the method to get the buses
    ampss = computer.compute_flows()
    res = 1.0 * ampss
    return res

### asyncio related
async def main(ls_grid, Vinit, prods_p, loads_p, loads_q):
    tasks = []
    loop = asyncio.get_running_loop()

    # inform the computation
    for th_id in range(NB_THREAD):  # see https://docs.python.org/3/library/asyncio-eventloop.html#asyncio.loop.run_in_executor
        task = loop.run_in_executor(None,
                                   get_flows,
                                   ls_grid, Vinit, prods_p[th_id], loads_p[th_id], loads_q[th_id])
        # task = asyncio.create_task(get_flows_async(ls_grid, Vinit, prods_p[th_id], loads_p[th_id], loads_q[th_id]))
        tasks.append(task)
    # run them
    for task in tasks:
        await task
    return [el.result() for el in tasks]
    

if __name__ == "__main__":
    env = get_env(ENV_NAME)
    prods_p, loads_p, loads_q = get_injs(env)
    Vinit = env.backend.V
    ls_grid = env.backend._grid

    # using threads
    beg_par_ = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=NB_THREAD) as executor:
        flows_future = {executor.submit(get_flows, ls_grid, Vinit, prods_p[th_id], loads_p[th_id], loads_q[th_id]): th_id
                        for th_id in range(NB_THREAD)}
        data_par = [None for _ in range(NB_THREAD)]
        for future in concurrent.futures.as_completed(flows_future):
            th_id = flows_future[future]
            data_par[th_id] = future.result()
    end_par_ = time.perf_counter()
    print(f"Execution time in parrallel using threads: {end_par_ - beg_par_:.2f}s")

    # using multiprocessing
    beg_mp_ = time.perf_counter()
    with Pool(processes=NB_THREAD) as pool:
        data_mp = pool.starmap(get_flows,
                               [(ls_grid, Vinit, prods_p[th_id], loads_p[th_id], loads_q[th_id]) for th_id in range(NB_THREAD)]
                              )
    end_mp_ = time.perf_counter()
    print(f"Execution time in parrallel using multiprocessing : {end_mp_ - beg_mp_:.2f}s")

    # using asyncio
    beg_asyncio_ = time.perf_counter()
    amps_asyncio = asyncio.run(main(ls_grid, Vinit, prods_p, loads_p, loads_q))
    end_asyncio_ = time.perf_counter()
    print(f"Execution time in parrallel using asyncio : {end_asyncio_ - beg_asyncio_:.2f}s")

    # sequential
    seq_prod_p = np.vstack(prods_p)
    seq_load_p = np.vstack(loads_p)
    seq_load_q = np.vstack(loads_q)
    beg_seq_ = time.perf_counter()
    ref_as = get_flows(ls_grid, Vinit, seq_prod_p, seq_load_p, seq_load_q)
    end_seq_ = time.perf_counter()
    print(f"Execution time (sequential) : {end_seq_ - beg_seq_:.2f}s")
    
    print(f"Using {NB_THREAD} \"cores\" the speed up where (sequential = 1)")
    print(f"Using multithreading: {(end_seq_ - beg_seq_) / (end_par_ - beg_par_):.2f}")
    print(f"Using multi processing: {(end_seq_ - beg_seq_) / (end_mp_ - beg_mp_):.2f}")
    print(f"Using asyncio: {(end_seq_ - beg_seq_) / (end_asyncio_ - beg_asyncio_):.2f}")

    # checking the results
    res_threading = np.vstack(data_par)
    res_mp = np.vstack(data_mp)
    res_asyncio = np.vstack(amps_asyncio)
    assert np.max(np.abs(res_threading - ref_as)) <= 1e-6, "error for multi threading"
    assert np.max(np.abs(res_mp - ref_as)) <= 1e-6, "error for multi processing"
    assert np.max(np.abs(res_asyncio - ref_as)) <= 1e-6, "error for asyncio"
