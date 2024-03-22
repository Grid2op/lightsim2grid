# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import setuptools
from setuptools import setup
import sys
import os
import warnings
from pybind11.setup_helpers import Pybind11Extension, build_ext


__version__ = "0.8.1.dev0"
KLU_SOLVER_AVAILABLE = False

# Try to link against SuiteSparse (if available)
# check that they exist (if SuiteSparse has been built with "make")
suitesparse_path_make = os.path.abspath("./SuiteSparse")
LIBS_MAKE = ["{}/KLU/Lib/libklu.a",
             "{}/BTF/Lib/libbtf.a",
             "{}/AMD/Lib/libamd.a",
             "{}/COLAMD/Lib/libcolamd.a",
             "{}/CXSparse/Lib/libcxsparse.a",
             "{}/SuiteSparse_config/libsuitesparseconfig.a"
             ]
LIBS_MAKE = [el.format(suitesparse_path_make) for el in LIBS_MAKE]
exists_libs_make = True
for el in LIBS_MAKE:
    if not os.path.exists(el):
        exists_libs_make = False

# check that they exist (if SuiteSparse has been built with "cmake" on macos / linux or windows)
suitesparse_path_cmake = os.path.abspath("./build_cmake/built/")
for ext in ["a", "lib"]:
    LIBS_CMAKE = [f"libklu.{ext}",
                  f"libbtf.{ext}",
                  f"libamd.{ext}",
                  f"libcolamd.{ext}",
                  f"libcxsparse.{ext}"]
    if ext == "a":
        # unix like system
        LIBS_CMAKE.append(f"libsuitesparseconfig.{ext}")
    else:
        # windows like system
        LIBS_CMAKE.append(f"suitesparseconfig.{ext}")

    LIBS_CMAKE = [os.path.join(suitesparse_path_cmake, "lib", el) for el in LIBS_CMAKE]

    exists_libs_cmake = True
    for el in LIBS_CMAKE:
        if not os.path.exists(el):
            exists_libs_cmake = False
            break

    if exists_libs_cmake:
        break


if exists_libs_make:
    # you will be able to use "SuiteSparse" and the faster "KLU" linear solver
    KLU_SOLVER_AVAILABLE = True

    # include directory
    INCLUDE_suitesparse = ["{}/SuiteSparse_config",
                           "{}/CXSparse/Include",
                           "{}/AMD/Include",
                           "{}/BTF/Include",
                           "{}/COLAMD/Include",
                           "{}/KLU/Include"
                           ]
    INCLUDE_suitesparse = [el.format(suitesparse_path_make) for el in INCLUDE_suitesparse]

    # compiled libraries location
    LIBS = LIBS_MAKE
elif exists_libs_cmake:
    # you will be able to use "SuiteSparse" and the faster "KLU" linear solver
    KLU_SOLVER_AVAILABLE = True

    # include directory
    INCLUDE_suitesparse = [os.path.join(suitesparse_path_cmake, "include", "suitesparse")]

    # compiled libraries location
    LIBS = LIBS_CMAKE
else:
    # SuiteSparse, and in particular the KLU linear solver is not available.
    # we'll use a default solver (a bit slower)
    LIBS = []
    INCLUDE_suitesparse = []
    warnings.warn("SuiteSparse is not available on your system, or has not been compiled. The faster "
                  "\"KLU\" linear algebra solver will not be available. The \"SparseLU\" solver will however "
                  "be available, which is maybe ~30% slower than \"KLU\". If you are using grid2op there "
                  "will still be a huge benefit.")

INCLUDE = ["src"] + INCLUDE_suitesparse

# now add the Eigen library (header only)
eigen_path = os.path.abspath(".")
INCLUDE.append("{}/eigen".format(eigen_path))

include_dirs = []
include_dirs += INCLUDE

# extract the version information
VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR = __version__.split(".")[:3]

# compiler options
extra_compile_args_tmp = ["-DNDEBUG"]
IS_LINUX = False
IS_MACOS = False
IS_WINDOWS = False
if sys.platform.startswith('linux'):
    # extra_compile_args_tmp = ["-fext-numeric-literals"]
    # -fext-numeric-literals is used for definition of complex number by some version of gcc
    extra_compile_args_tmp += []
    IS_LINUX = True
elif sys.platform.startswith("darwin"):
    # extra_compile_args_tmp = ["-fsized-deallocation"]
    extra_compile_args_tmp += []
    # fix a bug in pybind11
    # https://github.com/pybind/pybind11/issues/1604
    IS_MACOS = True
elif sys.platform.startswith("win32"):
    extra_compile_args_tmp += [# otherwise windows compiler does not import "M_PI" from the math header
                               "-D_USE_MATH_DEFINES"]
    IS_WINDOWS = True

# if you have installed some BLAS or LAPACKE libraries (on ubuntu sudo apt-get install libblas-dev liblapacke-dev)
# you can also trigger their use when using eigen.
# extra_compile_args_tmp += ["-DEIGEN_USE_BLAS", "-DEIGEN_USE_LAPACKE"]
extra_compile_args = extra_compile_args_tmp
# add the version information
extra_compile_args += [f"-DVERSION_MAJOR={VERSION_MAJOR}",
                       f"-DVERSION_MEDIUM={VERSION_MEDIUM}",
                       f"-DVERSION_MINOR={VERSION_MINOR}"]
extra_compile_args += [f"-DVERSION=\"{__version__}\""]
src_files = ['src/main.cpp',
             "src/powerflow_algorithm/GaussSeidelAlgo.cpp",
             "src/powerflow_algorithm/GaussSeidelSynchAlgo.cpp",
             "src/powerflow_algorithm/BaseAlgo.cpp",
             "src/linear_solvers/SparseLUSolver.cpp",
             "src/help_fun_msg.cpp",
             "src/BaseConstants.cpp",
             "src/GridModel.cpp",
             "src/ChooseSolver.cpp",
             "src/Solvers.cpp",
             "src/Utils.cpp",
             "src/DataConverter.cpp",
             "src/batch_algorithm/BaseBatchSolverSynch.cpp",
             "src/batch_algorithm/TimeSeries.cpp",
             "src/batch_algorithm/ContingencyAnalysis.cpp",
             "src/element_container/LineContainer.cpp",
             "src/element_container/GenericContainer.cpp",
             "src/element_container/ShuntContainer.cpp",
             "src/element_container/TrafoContainer.cpp",
             "src/element_container/LoadContainer.cpp",
             "src/element_container/GeneratorContainer.cpp",
             "src/element_container/SGenContainer.cpp",
             "src/element_container/DCLineContainer.cpp"]

if KLU_SOLVER_AVAILABLE:
    src_files.append("src/linear_solvers/KLUSolver.cpp")
    extra_compile_args_tmp.append("-DKLU_SOLVER_AVAILABLE")
    print("INFO: Using KLU package")

# Try to locate the NICSLU sparse linaer solver
if "PATH_NICSLU" in os.environ:
    # user indicate the path for the NICSLU library (see https://github.com/chenxm1986/nicslu)
    # eg "/home/user/Documents/nicslu/nicslu202103/"

    path_nicslu = os.path.abspath(os.environ["PATH_NICSLU"])
    include_nicslu = True
    # check for appropriate license
    if not os.path.exists(path_nicslu):
        print(f"WARNING: nothing for NICSLU at at: {path_nicslu}")
        include_nicslu = False
    license_path = os.path.join(path_nicslu, "license")
    if include_nicslu and not os.path.exists(license_path):
        # license not located at the right directory
        print(f"WARNING: no license path found for NICSLU at: {license_path}")
        include_nicslu = False
    license_file = os.path.join(license_path, "nicslu.lic")
    if include_nicslu and not os.path.exists(license_file):
        # no license found
        print(f"WARNING: no license file is found for NICSLU at: {license_file}")
        include_nicslu = False
    libnicslu_path = None
    if include_nicslu:
        if sys.platform.startswith("linux") or sys.platform.startswith("darwin"):
            libnicslu_path = os.path.join(path_nicslu, "linux/lib_centos6_x64_gcc482_fma/int32/libnicslu.so")
            if not os.path.exists(libnicslu_path):
                print(f"WARNING: cannot locate the NICSLU shared object that should be at: {libnicslu_path}")
                include_nicslu = False
                libnicslu_path = None
        elif sys.platform.startswith("win"):
            libnicslu_path = os.path.join(path_nicslu, "windows/lib_win7_x64_fma/int32/nicslu.lib")
            if not os.path.exists(libnicslu_path):
                print(f"WARNING: cannot locate the NICSLU shared object that should be at: {libnicslu_path}")
                include_nicslu = False
                libnicslu_path = None
        else:
            print(f"WARNING: NICSLU can only be added when using linux, darwin (MacOS) or win (Windows) python version, you are using {sys.platform}")
            include_nicslu = False
            libnicslu_path = None

    if include_nicslu and libnicslu_path is not None:
        LIBS.append(os.path.join(path_nicslu, libnicslu_path))
        include_dirs.append(os.path.join(path_nicslu, "include"))
        src_files.append("src/linear_solvers/NICSLUSolver.cpp")
        extra_compile_args.append("-DNICSLU_SOLVER_AVAILABLE")
        extra_compile_args += [f"-DNICSLU_PATH=\"{os.path.join(path_nicslu, libnicslu_path)}\""]
        print("INFO: Using NICSLU package")

# Try to locate the CKTSO sparse linear solver
if "PATH_CKTSO" in os.environ:
    # user indicate the path for the CKTSO library (see https://github.com/chenxm1986/cktso)
    # eg "/home/user/Documents/cktso/"
    
    path_cktso = os.path.abspath(os.environ["PATH_CKTSO"])
    include_cktso = True
    # check for appropriate license
    if not os.path.exists(path_cktso):
        print(f"WARNING: nothing for CKTSO at at: {path_cktso}")
        include_cktso = False
    license_path = os.path.join(path_cktso, "license")
    if include_cktso and not os.path.exists(license_path):
        # license not located at the right directory
        print(f"WARNING: no license path found for CKTSO at: {license_path}")
        include_cktso = False
    license_file = os.path.join(license_path, "cktso.lic")
    if include_cktso and not os.path.exists(license_file):
        # no license found
        print(f"WARNING: no license file is found for CKTSO at: {license_file}")
        include_cktso = False
    libcktso_path = None
    if include_cktso:
        if sys.platform.startswith("linux") or sys.platform.startswith("darwin"):
            libcktso_path = os.path.join(path_cktso, "ubuntu1804_x64_gcc750/libcktso.so")  # TODO CHANGE THAT !
            if not os.path.exists(libcktso_path):
                print(f"WARNING: cannot locate the CKTSO shared object that should be at: {libcktso_path}")
                include_cktso = False
                libcktso_path = None
        elif sys.platform.startswith("win"):
            libcktso_path = os.path.join(path_cktso, "win7_x64/cktso.lib")  # TODO CHANGE THAT !
            if not os.path.exists(libcktso_path):
                print(f"WARNING: cannot locate the CKTSO shared object that should be at: {libcktso_path}")
                include_cktso = False
                libcktso_path = None
        else:
            print(f"WARNING: CKTSO can only be added when using linux, darwin (MacOS) or win (Windows) python version, you are using {sys.platform}")
            include_cktso = False
            libcktso_path = None

    if include_cktso and libcktso_path is not None:
        LIBS.append(os.path.join(path_cktso, libcktso_path))
        include_dirs.append(os.path.join(path_cktso, "include"))
        src_files.append("src/linear_solvers/CKTSOSolver.cpp")
        extra_compile_args.append("-DCKTSO_SOLVER_AVAILABLE")
        extra_compile_args += [f"-DCKTSO_PATH=\"{os.path.join(path_cktso, libcktso_path)}\""]
        print("INFO: Using CKTSO package")
        
        
if "__COUT_TIMES" in os.environ:
    # to add extra info in cout for the computation times, we do not recommend to use it !
    if os.environ["__COUT_TIMES"] == "1":
        extra_compile_args.append("-D__COUT_TIMES")
        print("WARNING: Using the \"cout times\" compatibility mode, do not use the generated package outside of testing !")

if "__COMPILE_MARCHNATIVE" in os.environ:
    if os.environ["__COMPILE_MARCHNATIVE"] == "1":
        extra_compile_args.append("-march=native")
        extra_compile_args.append("-D__COMPILE_MARCHNATIVE")
        print("INFO: Using \"-march=native\" compiler flag")

# $Env:_READ_THE_DOCS = "1" in powershell
if "_READ_THE_DOCS" in os.environ:
    # generation is made on readthedocs.org for documentation, everything must be added, even though some packages will
    # not be available (eg KLU, NICSLU, etc.)

    if os.environ["_READ_THE_DOCS"] == "1":
        extra_compile_args.append("-D_READ_THE_DOCS")
        print("WARNING: Using the \"read the docs\" compatibility mode, do not use the package for something else than generating documentation.")
        
if False:
    path_iidm = ""
    lib_iidm = [os.path.join(path_iidm, "lib", "libiidm.a")]

    src_files.append("src/IIDMConverter.cpp")
    extra_compile_args_tmp.append("-DIIDM_CONVERTER_AVAILABLE")

    include_dirs.append(os.path.join(path_iidm, "include"))
    include_dirs.append("/usr/include/libxml2/")  # for libxml2
    LIBS += lib_iidm

if "__O3_OPTIM" in os.environ:
    # to add extra info in cout for the computation times, we do not recommend to use it !
    if os.environ["__O3_OPTIM"] == "1":
        if IS_LINUX or IS_MACOS:
            extra_compile_args.append("-O3")
            print("INFO: Using \"-O3\" compiler flag")
        elif IS_WINDOWS:
            extra_compile_args.append("/O2")
        extra_compile_args.append("-D__O3_OPTIM")

ext_modules = [
    Pybind11Extension(
        'lightsim2grid_cpp',
        src_files,
        include_dirs=include_dirs,
        extra_objects=LIBS,
        extra_compile_args=extra_compile_args
    )
]

# before pandapower 2.8 scipy was forced to be <= 1.6 which does not work with
# python 3.10+
req_pkgs = [
        "pandapower" if sys.version_info < (3, 10) else "pandapower>=2.8",
        "pytest",  # for pandapower see https://github.com/e2nIEE/pandapower/issues/1988
    ]

if sys.version_info.major == 3 and sys.version_info.minor <= 7:
    # typing "Literal" not available on python 3.7
    req_pkgs.append("typing_extensions")
    # do not use pandapower 2.12 (broken on python 3.7 
    # see https://github.com/e2nIEE/pandapower/issues/1985
    req_pkgs[0] = "pandapower>=2.2.2,<2.12"

pkgs = {
    "required": req_pkgs,
    "extras": {
        "docs": [
            "numpydoc>=0.9.2",
            "sphinx>=2.4.4,<7",
            "sphinx-rtd-theme>=0.4.3",
            "sphinxcontrib-trio>=1.1.0",
            "autodocsumm>=0.1.13",
            "grid2op>=1.6.4",
            "recommonmark",
            "pypowsybl"
        ],
        "benchmark": [
            "tabulate",
            "grid2op>=1.6.4",
            "numpy>=1.20",
            "distro",
            "py-cpuinfo"
        ],
        "recommended": [
            "grid2op>=1.6.4",
            "numba"
        ],
        "test": [
            "grid2op>=1.6.4",
            "numba",
            "pandapower>=2.8.0",
            "packaging", 
            "pypowsybl"
        ]
    }
}

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='LightSim2Grid',
      version=__version__,
      author='Benjamin Donnot',
      author_email='benjamin.donnot@rte-france.com',
      url='https://github.com/BDonnot/lightsim2grid/',
      description='LightSim2Grid implements a c++ backend targeting the Grid2Op platform.',
      long_description=long_description,
      long_description_content_type='text/markdown',
      ext_modules=ext_modules,
      install_requires=pkgs["required"],
      extras_require=pkgs["extras"],
      # setup_requires=['pybind11>=2.4'],  # in the pyproject.toml directly now
      cmdclass={'build_ext': build_ext},
      zip_safe=False,
      license='MPL 2.0',
      platforms=["Windows", "Linux", "Mac OS-X", "Unix"],
      packages=setuptools.find_packages(),
      keywords='pandapower powergrid simulator KLU Eigen c++',
      classifiers=[
            'Development Status :: 4 - Beta',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
            'Programming Language :: Python :: 3.12',
            "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "Natural Language :: English"
      ]
      )
