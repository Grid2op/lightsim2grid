# Copyright (c) 2023-2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import subprocess
import os
import shutil
from tqdm import tqdm
from tabulate import tabulate
import pandas as pd


numpy_ver = '1.24.4'  # latest on Apr. 22nd 2024 for all tested python
pds_ver = '2.0.3'  # latest on Apr. 22nd 2024 for all tested python
pyarrow_ver = '16.0.0 '  # latest on Apr. 22nd 2024 for all tested python
pp_ver = '2.14.6'  # latest on Apr. 22nd 2024 for all tested python

py_ver = "3.9"
g2op_ver = "1.10.1"

for py_ver in tqdm([# "3.8" , "3.9", 
                    "3.10", "3.11", "3.12"]):
    
    # create the venv (one for each python / lightsim2grid version, but reaused for grid2op)
    venv_nm = f"venv_py{py_ver}"
    # print(f"Creation of the virtual env for python {py_ver}")
    subprocess.run([f"python{py_ver}", "-m", "venv", venv_nm])
    py_exec = f"{venv_nm}/bin/python"

    my_env = {}
    my_env["VIRTUAL_ENV"] = f"{os.path.abspath('.')}/{venv_nm}"
    my_env["PATH"] = f"{my_env['VIRTUAL_ENV']}/bin:{os.environ['PATH']}"
    my_env["PATH_NICSLU"] = "/home/benjamin/Documents/powerflow_klu/nicslu/nicslu202110"
    my_env["PATH_CKTSO"] = "/home/benjamin/Documents/powerflow_klu/cktso"
    my_env["__COMPILE_MARCHNATIVE"] = "1"
    my_env["__O3_OPTIM"] = "1"

    # install everything in the venv
    # print(f"... Setting up the dependencies before installing lightsim2grid")
    res = subprocess.run([py_exec, "-m", "pip", "install", "--upgrade", 
                          f"grid2op=={g2op_ver}",
                          f"numpy=={numpy_ver}",
                          f"pandas=={pds_ver}",
                          f"pyarrow=={pyarrow_ver}",
                          f"pandapower=={pp_ver}",
                          "pybind11",
                          "tabulate",
                          "py-cpuinfo",
                          "distro"
                          ],
                          capture_output=True,
                          env=my_env)
    
    # prepare compilation of lightsim2grid
    # print(f"... Compile lightsim2grid")
    res = subprocess.run([py_exec, "-m", "pip", "install", ".."],
                          env=my_env,
                          capture_output=True)
    res = subprocess.run([py_exec, "-c", "import lightsim2grid; print(lightsim2grid.__version__)"],
                          env=my_env,
                          capture_output=True)
    ls_ver = "0.8.2"
    # ls_ver = res.stdout.decode("utf-8").lstrip().rstrip()
    
    for g2op_ver in tqdm([# "1.7.0",
                          # "1.7.1",
                          # "1.7.2", 
                          # "1.8.0",
                          # "1.8.1",
                          "1.9.0", 
                          "1.9.1", 
                          "1.9.2",
                          "1.9.3",
                          "1.9.4",
                          "1.9.5",
                          "1.9.6",
                          "1.9.7",
                          "1.9.8",
                          "1.10.0",
                          "1.10.1"],
                          leave=False):

        # install grid2op
        res = subprocess.run([py_exec, "-m", "pip", "install",
                              f"grid2op=={g2op_ver}",
                              "--no-deps"
                              ],
                              capture_output=True,
                              env=my_env)
        
        # fix the "dtype" in grid2op issue
        shutil.copyfile("../../grid2op/grid2op/dtypes.py",
                        f"{venv_nm}/lib/python{py_ver}/site-packages/grid2op/dtypes.py")
            
        # run the benchmarks
        res = subprocess.run([py_exec, "benchmark_solvers.py",
                              "--env_name", "l2rpn_case14_sandbox", 
                              "--no_test",
                              "--number", "1000",
                              "--save_results", f"py{py_ver}_gop{g2op_ver}_ls{ls_ver}_case14"],
                              capture_output=True,
                              env=my_env)
        res = subprocess.run([py_exec, "benchmark_solvers.py",
                              "--env_name", "l2rpn_neurips_2020_track2_small",                             
                              "--no_test",
                              "--number", "1000",
                              "--save_results", f"py{py_ver}_gop{g2op_ver}_ls{ls_ver}_case118"],
                              capture_output=True,
                              env=my_env)
        # import back the result
        df_res_14 = pd.read_csv(f"py{py_ver}_gop{g2op_ver}_ls{ls_ver}_case14"+"speed.csv", sep=";")
        df_diff_14 = pd.read_csv(f"py{py_ver}_gop{g2op_ver}_ls{ls_ver}_case14"+"diff.csv", sep=";")
        df_res_118 = pd.read_csv(f"py{py_ver}_gop{g2op_ver}_ls{ls_ver}_case118"+"speed.csv", sep=";")
        df_diff_118 = pd.read_csv(f"py{py_ver}_gop{g2op_ver}_ls{ls_ver}_case118"+"diff.csv", sep=";")
        with open(f"py{py_ver}_gop{g2op_ver}_ls{ls_ver}_case14"+"config_info.txt", "r", encoding="utf-8") as f:
            header_14 = f.readlines()
        with open(f"py{py_ver}_gop{g2op_ver}_ls{ls_ver}_case118"+"config_info.txt", "r", encoding="utf-8") as f:
            header_118 = f.readlines()
        # shape them into a proper md file
        bench_file = [f"Lightsim2grid {ls_ver} and grid2op {g2op_ver} (python {py_ver})",
                      "=================================================================",
                      "",
                      "l2rpn_case14_sandbox",
                      "---------------------",
                      "",
                      "Configuration:",
                      "",
                      ]
        bench_file += [el.rstrip() for el in header_14]
        bench_file += [""]
        res_14 = tabulate(df_res_14,
                          headers=df_res_14.columns, 
                          tablefmt="rst",
                          showindex="never")
        bench_file += [res_14] + [""]
        diff_14 = tabulate(df_diff_14,
                           headers=df_diff_14.columns, 
                           tablefmt="rst",
                           showindex="never")
        bench_file += [diff_14] + [""]
        
        bench_file += ["l2rpn_neurips_2020_track2_small",
                       "---------------------------------",
                       "",
                       "Configuration:",
                       "",
                       ]
        bench_file += [el.rstrip() for el in header_118]
        bench_file += [""]
        res_118 = tabulate(df_res_118,
                           headers=df_res_118.columns, 
                           tablefmt="rst",
                           showindex="never")
        bench_file += [res_118] + [""]
        diff_118 = tabulate(df_diff_118,
                            headers=df_diff_118.columns, 
                            tablefmt="rst",
                            showindex="never")
        bench_file += [diff_118]
        # write this md file
        with open(f"ls{ls_ver}_glop{g2op_ver}_py{py_ver}.rst", "w", encoding="utf-8") as f:
            f.write("\n".join(bench_file))
            
    # remove the venv
    shutil.rmtree(venv_nm)
  