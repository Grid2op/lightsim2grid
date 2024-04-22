# Copyright (c) 2023, RTE (https://www.rte-france.com)
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

numpy_ver = '1.24.4'  # latest on sept. 8th 2023 for all tested python
pds_ver = '2.0.3'  # latest on sept. 8th 2023 for all tested python
pyarrow_ver = '13.0.0'  # latest on sept. 8th 2023 for all tested python
pp_ver = '2.13.1'  # latest on sept. 8th 2023 for all tested python

py_ver = "3.9"
g2op_ver = "1.7.0"
for py_ver in tqdm([# "3.8", 
                    "3.10", "3.11"]):
  for g2op_ver in tqdm([
                        "1.7.0", "1.7.1", "1.7.2", 
                        "1.8.0", "1.8.1",
                        "1.9.0", 
                        "1.9.1", 
                        "1.9.2", "1.9.3", "1.9.4"],
                       leave=False):
      # create the venv
      venv_nm = f"venv_py{py_ver}_g2op{g2op_ver}"
      subprocess.run([f"python{py_ver}", "-m", "venv", venv_nm])
      py_exec = f"{venv_nm}/bin/python"

      my_env = {}
      my_env["VIRTUAL_ENV"] = f"{os.path.abspath('.')}/venv_nm"
      my_env["PATH"] = f"{my_env['VIRTUAL_ENV']}/bin:{os.environ['PATH']}"
      my_env["PATH_NICSLU"] = "/home/benjamin/Documents/powerflow_klu/nicslu/nicslu202110"
      my_env["PATH_CKTSO"] = "/home/benjamin/Documents/powerflow_klu/cktso"
      my_env["__COMPILE_MARCHNATIVE"] = "1"
      my_env["__O3_OPTIM"] = "1"

      # install everything in the venv
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

      # fix the "dtype" in grid2op issue
      shutil.copyfile("../../grid2op/grid2op/dtypes.py",
                      f"{venv_nm}/lib/python{py_ver}/site-packages/grid2op/dtypes.py")

      # prepare compilation of lightsim2grid
      res = subprocess.run([py_exec, "-m", "pip", "install", ".."],
                            env=my_env,
                            capture_output=True)
      res = subprocess.run([py_exec, "-c", "import lightsim2grid; print(lightsim2grid.__version__)"],
                            env=my_env,
                            capture_output=True)
      ls_ver = res.stdout.decode("utf-8")

      # run the benchmarks
      res = subprocess.run([py_exec, "benchmark_solvers.py",
                            "--name", "l2rpn_case14_sandbox", 
                            "--no_test",
                            "--number", "1000",
                            "--save_results", f"py{py_ver}_gop{g2op_ver}_ls{ls_ver}_case14.md"],
                            capture_output=True,
                            env=my_env)
      res = subprocess.run([py_exec, "benchmark_solvers.py",
                            "--name", "l2rpn_neurips_2020_track2_small",                             
                            "--no_test",
                            "--number", "1000",
                            "--save_results", f"py{py_ver}_gop{g2op_ver}_ls{ls_ver}_case118.md"],
                            capture_output=True,
                            env=my_env)

      # remove the venv
      shutil.rmtree(venv_nm)
