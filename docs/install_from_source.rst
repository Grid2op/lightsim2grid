.. _page_install_from_source:

Installation from source
==========================

.. toctree::
  :maxdepth: 2
  :numbered:
  :caption: Main Steps

At a glance, to install from source you will have to:

- :ref:`clone_repo` clone this repository and get the code of Eigen (mandatory for compilation) and SuiteSparse (optional, but recommended)

  - :ref:`include_nicslu`: retrieve and get a proper license for the NICSLU linear solver (see https://github.com/chenxm1986/nicslu)
  - :ref:`include_cktso`: retrieve and get a proper license for the CKTSO linear solver (see https://github.com/chenxm1986/cktso)
  - :ref:`other_customization` specify some compilation flags to make the package run faster on your machine see 

- :ref:`install_python`: install the python package see 

.. warning::
    Everything here is better done inside a python "virtual environment" (also known as `venv`) or even better
    using the excellent uv package. 

    We only provide the commands related to lightsim2grid in this documentation. 

    It's probably safer and better to execute them in a venv, *eg* (on macos / linux) with :

    .. code-block::

        uv venv
        source .venv/bin/activate

    And then proceed as described below.


Important note
----------------
This package relies on the excellent `pybind11` package to integrate c++ code into python easily. 

So to install lightsim2grid you need `pybind11` and its requirement, which include a working compiler: for example 
(as of writing) 
gcc (default on ubuntu, version >= 4.8), clang (default on MacOS, version >= 5.0.0) or 
Microsoft visual studio (Microsoft Visual Studio 2015 Update 3 or newer). 

This readme does not cover the install of such compilers. Please refer to the documentation of 
`pybind11 <https://pybind11.readthedocs.io/en/latest/>`_ for more information. Do not hesitate to write github issues
if you encounter a problem in installing such compiler (**nb** on windows you have to install
visual studio, on linux of MacOs you might already have a working compiler installed).

.. _clone_repo:

1. Retrieve the sources
--------------------------
First, you can download it with git with:

.. code-block::

    git clone https://github.com/Grid2Op/lightsim2grid.git
    cd lightsim2grid

    # retrieve the code of SuiteSparse and Eigen (dependencies, mandatory)
    git submodule init
    git submodule update

The build itself is driven by `scikit-build-core <https://scikit-build-core.readthedocs.io/>`_ from
``pyproject.toml``, which in turn runs CMake to compile the ``lightsim2grid_core`` C++ library and its
pybind11 bindings (see :ref:`cpp_library` if you want to build and use that C++ library directly, without
python, or run its own unit test suite). SuiteSparse (for the faster ``KLU`` linear solver) is compiled
automatically as part of this process, now that its sources are available -- you don't have anything more
to do than run the install command in :ref:`install_python` below.

.. _include_nicslu:

(optional) Include NICSLU linear solver (experimental)
++++++++++++++++++++++++++++++++++++++++++++++++++++++

Another linear solver that can be used with lighsim2grid is the "NICSLU" linear solver that might, in some cases, be
even faster than the KLU linear solver. This can lead to more speed up if using lighsim2grid.

To use it, you need to:

1) retrieve the sources (only available as a freeware) from https://github.com/chenxm1986/nicslu and save
   it on your machine. Say you clone this github repository in `NICSLU_GIT` 
   (*eg* NICSLU_GIT="/home/user/Documents/nicslu/"). Also note that you need to check that your usage
   is compliant with their license !
2) define the "PATH_NICSLU" environment variable **before** compiling lightsim2grid, on linux you can do
   `export PATH_NICSLU=NICSLU_GIT/nicsluDATE` 
   (for example `export PATH_NICSLU=/home/user/Documents/nicslu/nicslu202103` if you cloned the repository 
   as the example of `step 1)` and use the version of nicslu compiled by the author on March 2021 [version 
   distributed at time of writing the readme] )

And this is it. Lightsim will be able to use this linear solver.

Be carefull though, you require a license file in order to use it. As of now, the best way is to copy paste the license
file at the same location that the one you execute python from (*ie* you need to copy paste it each time).

.. _include_cktso:

(optional) Include CKTSO linear solver (experimental)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Another linear solver that can be used with lighsim2grid is the "CKTSO" linear solver (a newer version of NICSLU) that might, 
in some cases, be even faster than the NICSLU linear solver. This can lead to more speed up if using lighsim2grid.

To use it, you need to:

1) retrieve the sources (only available as a freeware) from https://github.com/chenxm1986/cktso and save
   it on your machine. Say you clone this github repository in `$CKTSO_GIT` 
   (*eg* CKTSO_GIT="/home/user/Documents/cktso/"). Also note that you need to check that your usage
   is compliant with their license !
2) define the "PATH_CKTSO" environment variable **before** compiling lightsim2grid, on linux you can do
   `export PATH_CKTSO=$CKTSO_GIT` 
   (for example `export PATH_CKTSO=/home/user/Documents/cktso` if you cloned the repository 
   as the example of `step 1`)

And this is it. Lightsim will be able to use this linear solver.

Be carefull though, you require a license file in order to use it. As of now, the best way is to copy paste the license
file at the same location as the library used when ligthsim2grid is compiled. It should be on `$CKTSO_GIT/ubuntu1804_x64_gcc750`
on linux and `$CKTSO_GIT/win7_x64` on windows.

.. _other_customization:

(optional) customization of the installation
++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you bother to compile from source the package, you might also want to benefit from some
extra speed ups.

This can be achieve by specifying the `__O3_OPTIM` and `__COMPILE_MARCHNATIVE` environment variables.

The first one controls whether the package is compiled using the `-O3` compiler flag (`/O2` on windows) which tells the compiler to optimize the code for speed even more. It is **enabled by default**.

The second one will compile the package using the `-march=native` flag (on macos and linux). It is **disabled by default**.

And example to enable `-march=native` on a linux based machine is:

.. code-block::

    export __COMPILE_MARCHNATIVE=1

If you want to disable the O3 optimization (for example to reduce compilation time), you need to explicitly set it to 0:

.. code-block::

    export __O3_OPTIM=0

.. note::
    By default the package (including the one on pypi) is compiled with O3 optimization activated.

.. _install_python:

2. Installation of the python package
----------------------------------------

Now you simply need to install the lightsim2grid package this way, like any python package:

.. code-block::

    # compile and install the python package
    pip install -U .

And you are done :-)

