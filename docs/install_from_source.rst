Installation from source
==========================

.. toctree::
  :maxdepth: 2
  :numbered:
  :caption: Main Steps

At a glance, to install from source you will have to:

- :ref:`clone_repo` clone this repository and get the code of Eigen (mandatory for compilation) and SparseSuite (optional, but recommended)

  - :ref:`compile_suitesparse` compile a piece of SparseSuite
  - :ref:`include_nicslu`: retrieve and get a proper license for the NICSLU linear solver (see https://github.com/chenxm1986/nicslu)
  - :ref:`include_cktso`: retrieve and get a proper license for the CKTSO linear solver (see https://github.com/chenxm1986/cktso)
  - :ref:`other_customization` specify some compilation flags to make the package run faster on your machine see 

- :ref:`install_python`: install the python package see 


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
    # it is recommended to do a python virtual environment
    python -m virtualenv venv  # optional
    source venv/bin/activate  # optional

    # retrieve the code of SparseSuite and Eigen (dependencies, mandatory)
    git submodule init
    git submodule update

.. _compile_suitesparse:

(optional, recommended) Compilation of SuiteSparse
++++++++++++++++++++++++++++++++++++++++++++++++++++++

SuiteSparse comes with the faster KLU linear solver. 

Since version 0.3.0 this requirement has been removed. This entails
that on linux / macos you can still benefit from the faster KLU solver. You can still benefit from the
speed up of lightsim (versus the default PandaPowerBackend) but this speed up will be less than if you manage
to compile SuiteSparse (see the subsection [Benchmark](#benchmark) for more information).

**NB** in both cases the algorithm to compute the powerflow is exactly the same. It is a 
Newton-Raphson based method. But to carry out this algorithm, one need to solver some linear equations. The only
difference in the two version (with KLU and without) is that the linear equation solver is different. Up to the
double float precision, both results (with and without KLU) should match.

There are 2 ways to install this package. Either you use "make" (preferred method on linux / unix -- including MacOS) or you use "cmake", which works on all platforms but takes more time and is less automatic (mainly because SuiteSparse
cannot be directly built with "cmake" so we need extra steps to make it possible.)

(optional) option A. Compilation of SuiteSparse using "make"
**************************************************************

This is the easiest method to compile SuiteSparse on your system but unfortunately it only works on OS where "make" is
available (*eg* Linux or MacOS) but this will not work on Windows... The compilation on windows is covered in the next
paragraph :ref:`optionB`

Anyway, in this case, it's super easy. Just do:

.. code-block::

    make

And you are good to go. Nothing more.

.. _optionB:

(optional) option B. Compilation of SuiteSparse using "cmake"
**************************************************************

This works on most platform including MacOS, Linux and Windows.

It requires to install the free `cmake` program and to do a bit more works than for other system. This is why we
only recommend to use it on Windows.

The main steps (for windows, somme commands needs to be adapted on linux / macos) are:

1) `cd build_cmake`
2) `py generate_c_files.py`
3) `mkdir build` and cd there: `cd build`
4) `cmake -DCMAKE_INSTALL_PREFIX=..\built -DCMAKE_BUILD_TYPE=Release ..`
5) `cmake --build . --config Release`
6) `cmake --build . --config Release --target install`

For more information, feel free to read the dedicated [README](build_cmake/README.md).

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

The first one will compile the package using the `-O3` compiler flag (`/O2` on windows) which will tell the compiler to optimize the code for speed even more.

The second one will compile the package using the `-march=native` flag (on macos and linux)

And example to do such things on a linux based machine is:

.. code-block::

    export __O3_OPTIM=1
    export __COMPILE_MARCHNATIVE=1

If you want to disable them, you simply need to set their respective value to 0.

.. _install_python:

2. Installation of the python package
----------------------------------------

Now you simply need to install the lightsim2grid package this way, like any python package:

.. code-block::

    # install the dependency
    pip install -U pybind11
    # compile and install the python package
    pip install -U .

And you are done :-)
