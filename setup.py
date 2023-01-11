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


__version__ = "0.7.1"
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

INCLUDE = INCLUDE_suitesparse

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
src_files = ['src/main.cpp',
             "src/help_fun_msg.cpp",
             "src/SparseLUSolver.cpp",
             "src/BaseConstants.cpp",
             "src/GridModel.cpp",
             "src/DataConverter.cpp",
             "src/DataLine.cpp",
             "src/DataGeneric.cpp",
             "src/DataShunt.cpp",
             "src/DataTrafo.cpp",
             "src/DataLoad.cpp",
             "src/DataGen.cpp",
             "src/DataSGen.cpp",
            #  "src/BaseNRSolver.cpp",  # moved as a template class
            #  "src/BaseNRSolverSingleSlack.cpp",  # moved as a template class
            #  "src/DCSolver.cpp",  # moved as a template class
             "src/ChooseSolver.cpp",
             "src/GaussSeidelSolver.cpp",
             "src/GaussSeidelSynchSolver.cpp",
             "src/BaseSolver.cpp",
             "src/BaseMultiplePowerflow.cpp",
             "src/Computers.cpp",
             "src/SecurityAnalysis.cpp"]

if KLU_SOLVER_AVAILABLE:
    src_files.append("src/KLUSolver.cpp")
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
        src_files.append("src/NICSLUSolver.cpp")
        extra_compile_args.append("-DNICSLU_SOLVER_AVAILABLE")
        print("INFO: Using NICSLU package")

if "__COUT_TIMES" in os.environ:
    # to add extra info in cout for the computation times, we do not recommend to use it !
    if os.environ["__COUT_TIMES"] == "1":
        extra_compile_args.append("-D__COUT_TIMES")
        print("WARNING: Using the \"cout times\" compatibility mode, do not use the generated package outside of testing !")

if "__COMPILE_MARCHNATIVE" in os.environ:
    if os.environ["__COMPILE_MARCHNATIVE"] == "1":
        extra_compile_args.append("-march=native")
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
    ]

pkgs = {
    "required": req_pkgs,
    "extras": {
        "docs": [
            "numpydoc>=0.9.2",
            "sphinx>=2.4.4",
            "sphinx-rtd-theme>=0.4.3",
            "sphinxcontrib-trio>=1.1.0",
            "autodocsumm>=0.1.13",
            "grid2op>=1.6.4",
            "recommonmark",
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
            "pandapower>=2.8.0"
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
            "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "Natural Language :: English"
      ]
      )
