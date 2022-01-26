# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import warnings
from numpy.core.fromnumeric import size
import pandapower as pp
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from grid2op import make, Parameters
from grid2op.Chronics import ChangeNothing
from lightsim2grid import LightSimBackend, TimeSerie
from tqdm import tqdm
import os

VERBOSE = False

case_names = ["case118.json",
              "case_illinois200.json",
              "case300.json",
              "case1354pegase.json",
              "case1888rte.json",
              "GBnetwork.json",  # 2224 buses
              "case2848rte.json",
              "case2869pegase.json",
              "case3120sp.json",
              "case6495rte.json",
              "case9241pegase.json"
              ]

def make_grid2op_env(pp_case, casse_name):
    param = Parameters.Parameters()
    param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})
        
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env_lightsim = make("blank",
                            param=param, test=True,
                            backend=LightSimBackend(),
                            data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                            grid_path=case_name,
                            _add_to_name=f"{case_name}")
    return env_lightsim

def get_loads_gens(load_p_init, load_q_init, gen_p_init, sgen_p_init):
    # scale loads

    # use some French time series data for loads
    # see https://github.com/BDonnot/data_generation for where to find this file
    coeffs = {"sources": {
    "country": "France",
    "year": "2012",
    "web": "http://clients.rte-france.com/lang/fr/visiteurs/vie/vie_stats_conso_inst.jsp"
    },
    "month": {
    "jan": 1.21,
    "feb": 1.40,
    "mar": 1.05,
    "apr": 1.01,
    "may": 0.86,
    "jun": 0.84,
    "jul": 0.84,
    "aug": 0.79,
    "sep": 0.85,
    "oct": 0.94,
    "nov": 1.01,
    "dec": 1.20
    },
    "day": {
    "mon": 1.01,
    "tue": 1.05,
    "wed": 1.05,
    "thu": 1.05,
    "fri": 1.03,
    "sat": 0.93,
    "sun": 0.88
    },
    "hour": {
    "00:00": 1.00,
    "01:00": 0.93,
    "02:00": 0.91,
    "03:00": 0.86,
    "04:00": 0.84,
    "05:00": 0.85,
    "06:00": 0.90,
    "07:00": 0.97,
    "08:00": 1.03,
    "09:00": 1.06,
    "10:00": 1.08,
    "11:00": 1.09,
    "12:00": 1.09,
    "13:00": 1.09,
    "14:00": 1.06,
    "15:00": 1.03,
    "16:00": 1.00,
    "17:00": 1.00,
    "18:00": 1.04,
    "19:00": 1.09,
    "20:00": 1.05,
    "21:00": 1.01,
    "22:00": 0.99,
    "23:00": 1.03
    }
    }
    vals = list(coeffs["hour"].values())
    x_final = np.arange(12*len(vals))

    # interpolate them at 5 minutes resolution (instead of 1h)
    vals.append(vals[0])
    vals = np.array(vals) * coeffs["month"]["mar"] * coeffs["day"]["mon"]
    x_interp = 12 * np.arange(len(vals))
    coeffs = interp1d(x=x_interp, y=vals, kind="cubic")
    all_vals = coeffs(x_final)

    # compute the "smooth" loads matrix
    load_p_smooth = all_vals.reshape(-1, 1) * load_p_init.reshape(1, -1)
    load_q_smooth = all_vals.reshape(-1, 1) * load_q_init.reshape(1, -1)

    # add a bit of noise to it to get the "final" loads matrix
    load_p = load_p_smooth * np.random.lognormal(mean=0., sigma=0.003, size=load_p_smooth.shape)
    load_q = load_q_smooth * np.random.lognormal(mean=0., sigma=0.003, size=load_q_smooth.shape)

    # scale generators accordingly
    gen_p = load_p.sum(axis=1).reshape(-1, 1) / load_p_init.sum() * gen_p_init.reshape(1, -1)
    sgen_p = load_p.sum(axis=1).reshape(-1, 1) / load_p_init.sum() * sgen_p_init.reshape(1, -1)

    return load_p, load_q, gen_p, sgen_p


if __name__ == "__main__":
    np.random.seed(42)
    times = []
    speeds = []
    sizes = []
    for case_name in tqdm(case_names):

        if not os.path.exists(case_name):
            import pandapower.networks as pn
            case = getattr(pn, os.path.splitext(case_name)[0])()
            pp.to_json(case, case_name)

        # load the case file
        case = pp.from_json(case_name)
        # extract reference data
        load_p_init = 1.0 * case.load["p_mw"].values
        load_q_init = 1.0 * case.load["q_mvar"].values
        gen_p_init = 1.0 * case.gen["p_mw"].values
        sgen_p_init = 1.0 * case.sgen["p_mw"].values

        res_time = 1.
        res_unit = "s"
        if len(load_p_init) <= 1000:
            # report results in ms if there are less than 1000 loads
            # only affects "verbose" printing
            res_time = 1e3
            res_unit = "ms"

        # simulate the data
        load_p, load_q, gen_p, sgen_p = get_loads_gens(load_p_init, load_q_init, gen_p_init, sgen_p_init)
        
        # create the grid2op env
        env_lightsim = make_grid2op_env(case, case_name)

        # compute the right things
        time_serie = TimeSerie(env_lightsim)
        computer = time_serie.computer
        v_init = env_lightsim.backend.V
        status = computer.compute_Vs(gen_p,
                                    sgen_p,
                                    load_p,
                                    load_q,
                                    v_init,
                                    env_lightsim.backend.max_it,
                                    env_lightsim.backend.tol)
        time_serie._TimeSerie__computed = True
        a_or = time_serie.compute_A()
        assert status, f"some powerflow diverge for {case_name}"

        if VERBOSE:
            # print detailed results if needed
            print(f"For environment: {case_name} ({env_lightsim.n_sub} substations) [{computer.nb_solved()} powerflows]")
            print(f"Total time spent in \"computer\" to solve everything: {res_time*computer.total_time():.2f}{res_unit} "
                f"({computer.nb_solved() / computer.total_time():.0f} pf / s), "
                f"{1000.*computer.total_time() / computer.nb_solved():.2f} ms / pf)")
            print(f"\t - time to pre process the injections: {res_time * computer.preprocessing_time():.2f}{res_unit}")
            print(f"\t - time to perform powerflows: {res_time * computer.solver_time():.2f} {res_unit} "
                f"({computer.nb_solved() / computer.solver_time():.0f} pf / s, "
                f"{1000.*computer.solver_time() / computer.nb_solved():.2f} ms / pf)")
            print(f"In addition, it took {res_time * computer.amps_computation_time():.2f} {res_unit} to retrieve the current "
                f"from the complex voltages (in total "
                f"{computer.nb_solved() / ( computer.total_time() + computer.amps_computation_time()):.1f} "
                "pf /s, "
                f"{1000.*( computer.total_time() + computer.amps_computation_time()) / computer.nb_solved():.2f} ms / pf)")
        
        times.append(computer.total_time() + computer.amps_computation_time())
        speeds.append(computer.nb_solved() / (computer.total_time() + computer.amps_computation_time()) )
        sizes.append(env_lightsim.n_sub)
        env_lightsim.close()

    # make the plot summarizing all results
    plt.plot(sizes, times, linestyle='solid', marker='+', markersize=8)
    plt.xlabel("Size (number of substation)")
    plt.ylabel("Time taken (s)")
    plt.title(f"Time to compute {computer.nb_solved()} powerflows")
    plt.show()

    plt.plot(sizes, speeds, linestyle='solid', marker='+', markersize=8)
    plt.xlabel("Size (number of substation)")
    plt.ylabel("Speed (pf / s)")
    plt.title(f"Computation speed")
    plt.yscale("log")
    plt.show()
