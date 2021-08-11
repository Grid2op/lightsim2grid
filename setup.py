from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os
import warnings
from pybind11.setup_helpers import Pybind11Extension, build_ext

__version__ = "0.5.3"
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
if sys.platform.startswith('linux'):
    # extra_compile_args_tmp = ["-fext-numeric-literals"]
    # -fext-numeric-literals is used for definition of complex number by some version of gcc
    extra_compile_args_tmp += []
elif sys.platform.startswith("darwin"):
    # extra_compile_args_tmp = ["-fsized-deallocation"]
    extra_compile_args_tmp += []
    # fix a bug in pybind11
    # https://github.com/pybind/pybind11/issues/1604
elif sys.platform.startswith("win32"):
    extra_compile_args_tmp += [# otherwise windows compiler does not import "M_PI" from the math header
                               "-D_USE_MATH_DEFINES"]


# for even greater speed, you can add the "-march=native" flag. It does not work on all platform, that is
# why we deactivated it by default
# extra_compile_args_tmp += ["-march=native"]

# if you have installed some BLAS or LAPACKE libraries (on ubuntu sudo apt-get install libblas-dev liblapacke-dev)
# you can also trigger their use when using eigen.
# extra_compile_args_tmp += ["-DEIGEN_USE_BLAS", "-DEIGEN_USE_LAPACKE"]
extra_compile_args = extra_compile_args_tmp
# add the version information
extra_compile_args += [f"-DVERSION_MAJOR={VERSION_MAJOR}",
                       f"-DVERSION_MEDIUM={VERSION_MEDIUM}",
                       f"-DVERSION_MINOR={VERSION_MINOR}"]
src_files = ['src/main.cpp',
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
             "src/BaseNRSolver.cpp",
             "src/ChooseSolver.cpp",
             "src/GaussSeidelSolver.cpp",
             "src/GaussSeidelSynchSolver.cpp",
             "src/BaseSolver.cpp",
             "src/DCSolver.cpp"]

if KLU_SOLVER_AVAILABLE:
    src_files.append("src/KLUSolver.cpp")
    extra_compile_args_tmp.append("-DKLU_SOLVER_AVAILABLE")

ext_modules = [
    Pybind11Extension(
        'lightsim2grid_cpp',
        src_files,
        include_dirs=include_dirs,
        extra_objects=LIBS,
        extra_compile_args=extra_compile_args
    )
]

pkgs = {
    "required": [
        "pandapower"
    ],
    "extras": {
        "docs": [
            "numpydoc>=0.9.2",
            "sphinx>=2.4.4",
            "sphinx-rtd-theme>=0.4.3",
            "sphinxcontrib-trio>=1.1.0",
            "autodocsumm>=0.1.13",
            # "m2r"
            "recommonmark",
        ],
        "benchmark": [
            "tabulate",
            "grid2op>=1.5.0",
            "numpy"
        ],
        "recommended": [
            "grid2op>=1.5.0"
        ],
        "test": [
            "grid2op>=1.5.0"
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
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "Natural Language :: English"
      ]
      )
