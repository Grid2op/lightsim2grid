# LightSim2Grid
Provide a fast backend for grid2op using c++ KLU and Eigen librairies. Its primary goal is to serve as a fast
backend for the grid2op platform, used primarily as a testbed platform for sequential decision making in
the world of power system.

See the [Disclaimer](DISCLAIMER.md) to have a more detailed view on what is and what is not this package. For example
this package should not be used for detailed power system computations or simulations.


*   [1 Usage](#Usage)
    *   [1.1. As a grid2op backend (preferred method)](#1-as-a-grid2op-backend-preferred-method)
    *   [1.2. replacement of pandapower "newtonpf" method (advanced method)](#2-replacement-of-pandapower-newtonpf-method-advanced-method)
*   [2 Installation (from pypi official repository, recommended)](#Installation-from-pypi-official-repository-recommended)
*   [3 Installation (from source, for more advanced user)](#Installation-from-source-for-more-advanced-user)
    *   [3.0 Important note](#Important-note)
    *   [3.1. Retrieve the sources](#1-Retrieve-the-sources)
    *   [(optional, recommended) Compilation of SuiteSparse](#optional-recommended-Compilation-of-SuiteSparse)
        *   [option A. Compilation of SuiteSparse using "make"](#optional-option-a-Compilation-of-SuiteSparse-using-make)
        *   [option B. Compilation of SuiteSparse using "cmake"](#optional-option-b-Compilation-of-SuiteSparse-using-cmake)
    *   [(optional) Include NICSLU linear solver (experimental)](#optional-Include-NICSLU-linear-solver-experimental)
    *   [(optional) customization of the installation](#optional-customization-of-the-installation)
    *   [3.2 Installation of the python package](#2-Installation-of-the-python-package)
*   [4. Benchmarks](#Benchmarks)
*   [5. Philosophy](#Philosophy)
*   [6. Miscellaneous](#Miscellaneous)
    * [6.1 Customization of the compilation](#Customization-of-the-compilation)
    * [6.2 Profile the code](#Profile-the-code)
    * [6.3 Local testing](#Local-testing)
    * [6.4 Tests performed automatically](#Tests-performed-automatically)
    * [5.5 Known issues](#Known-issues)

## Usage
Once installed (don't forget, if you used the optional virtual env
above you need to load it with `source venv/bin/activate`) you can
use it as any python package.

### 1. As a grid2op backend (preferred method)
This functionality requires you to have grid2op installed, with at least version 0.7.0. You can install it with

```commandline
pip install grid2op>=1.6.4
```

Then you can use a LightSimBackend instead of the default PandapowerBackend this way:

```python3
import grid2op
from lightsim2grid import LightSimBackend
env_name = "l2rpn_case14_sandbox"  # or any other name.
env = grid2op.make(env_name, backend=LightSimBackend())

# do regular computation as you would with grid2op
```
And you are good to go.

### 2. replacement of pandapower "newtonpf" method (advanced method)
It is also possible to use directly the "solver" part of lightsim2grid.

Suppose you somehow get:
- `Ybus` the admittance matrix of your powersystem, for example given by pandapower
  (will be converted to a scipy `sparse.csc_matrix` )
- `V0` the (complex) voltage vector at each bus, for example given by pandapower
- `Sbus` the (complex) power absorb at each bus, for example as given by pandapower
- `ref` Ids of the slack buses (added in version 0.5.6 to match recent pandapower changes)
- `pv` list of PV buses
- `pq` list of PQ buses
- `ppci` a ppc internal pandapower test case (or dictionary, is used to retrieve the coefficients associated to each slack bus)
- `options` list of pandapower "options" (or dictionary with keys `max_iteration` and `tolerance_mva`)

You can define replace the `newtonpf` function of `pandapower.pandapower.newtonpf` function with the following
piece of code:
```python
from lightsim2grid.newtonpf import newtonpf
V, converged, iterations, J = newtonpf(Ybus, V, Sbus, ref, weights, pv, pq, ppci, options)
```

This function uses the KLU algorithm and a c++ implementation of a Newton solver for speed.

## Installation (from pypi official repository, recommended)

Since version 0.5.3, lightsim2grid is can be installed like most python packages, with a call to:
`python -m pip install lightsim2grid`

It includes faster grid2op backend and the `SuiteSparse` faster `KLU` solver, even on windows. This is definitely the 
easiest method to install lightsim2grid on your system and have it running without too much issues.

Note though that these packages have been compiled on a different platform that the one you are using. You might still
get some benefit (in terms of performances) to install it from your on your machine.

Pypi packages are available for linux, windows and macos with python versions: 

- 3.7
- 3.8
- 3.9
- 3.10 (lightsim2grid >= 0.6.1)

## Installation (from source, for more advanced user)
You need to:
- clone this repository and get the code of Eigen (mandatory for compilation) and SparseSuite (optional, but recommended)
- (optional, but recommended) compile a piece of SparseSuite
- (optional) [experimental] retrieve and get a proper license for the NICSLU linear solver (see https://github.com/chenxm1986/nicslu)
- (optional) specify some compilation flags to make the package run faster on your machine
- install the package

### Important note
This package relies on the excellent `pybind11` package to integrate c++ code into python easily. 

So to install lightsim2grid you need `pybind11` and its requirement, which include a working compiler: for example 
(as of writing) 
gcc (default on ubuntu, version >= 4.8), clang (default on MacOS, version >= 5.0.0) or 
Microsoft visual studio (Microsoft Visual Studio 2015 Update 3 or newer). 

This readme does not cover the install of such compilers. Please refer to the documentation of 
[pybind11](https://pybind11.readthedocs.io/en/latest/) for more information. Do not hesitate to write github issues
if you encounter a problem in installing such compiler (**nb** on windows you have to install
visual studio, on linux of MacOs you might already have a working compiler installed).

### 1. Retrieve the sources
First, you can download it with git with:
```commandline
git clone https://github.com/BDonnot/lightsim2grid.git
cd lightsim2grid
# it is recommended to do a python virtual environment
python -m virtualenv venv  # optional
source venv/bin/activate  # optional

# retrieve the code of SparseSuite and Eigen (dependencies, mandatory)
git submodule init
git submodule update
```

### (optional, recommended) Compilation of SuiteSparse
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

#### (optional) option A. Compilation of SuiteSparse using "make"
This is the easiest method to compile SuiteSparse on your system but unfortunately it only works on OS where "make" is
available (*eg* Linux or MacOS) but this will not work on Windows... The compilation on windows is covered in the next
paragraph 
[(optional) option B. Compilation of SuiteSparse using "cmake"](#\(optional\)-option-B.-Compilation-of-SuiteSparse-using-"cmake")

Anyway, in this case, it's super easy. Just do:

```commandline
# compile static libraries of SparseSuite
make
```
And you are good to go. Nothing more.

#### (optional) option B. Compilation of SuiteSparse using "cmake"
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

### (optional) Include NICSLU linear solver (experimental)
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

### (optional) customization of the installation

If you bother to compile from source the package, you might also want to benefit from some
extra speed ups.

This can be achieve by specifying the `__O3_OPTIM` and `__COMPILE_MARCHNATIVE` environment variables. 

The first one will compile the package using the `-O3` compiler flag (`/O2` on windows) which will tell the compiler to optimize the code for speed even more.

The second one will compile the package using the `-march=native` flag (on macos and linux)

And example to do such things on a linux based machine is:

```commandline
export __O3_OPTIM=1
export __COMPILE_MARCHNATIVE=1
```

If you want to disable them, you simply need to set their respective value to 0.

### 2. Installation of the python package
Now you simply need to install the lightsim2grid package this way, like any python package:

```commandline
# install the dependency
pip install -U pybind11
# compile and install the python package
pip install -U .
```

And you are done :-)


## Benchmarks

Lightsim2grid is significantly faster than pandapower when used with grid2op for all kind of environment size.

First on an environment based on the IEEE case14 grid:

| case14_sandbox     |   grid2op speed (it/s) |   grid2op 'backend.runpf' time (ms) |   solver powerflow time (ms) |
|--------------------|------------------------|-------------------------------------|------------------------------|
| PP                 |                   70.5 |                              11     |                       4.27   |
| LS+GS              |                  881   |                               0.447 |                       0.327  |
| LS+GS S            |                  877   |                               0.446 |                       0.327  |
| LS+SLU (single)    |                 1110   |                               0.191 |                       0.0655 |
| LS+SLU             |                 1120   |                               0.195 |                       0.0683 |
| LS+KLU (single)    |                 1200   |                               0.138 |                       0.0176 |
| LS+KLU             |                 1180   |                               0.141 |                       0.0188 |
| LS+NICSLU (single) |                 1200   |                               0.139 |                       0.0179 |
| LS+NICSLU          |                 1200   |                               0.139 |                       0.0184 |

Then on an environment based on the IEEE case 118:

| neurips_2020_track2   |   grid2op speed (it/s) |   grid2op 'backend.runpf' time (ms) |   solver powerflow time (ms) |
|-----------------------|------------------------|-------------------------------------|------------------------------|
| PP                    |                   39.6 |                              13.3   |                        5.58  |
| LS+GS                 |                    5.3 |                             188     |                      188     |
| LS+GS S               |                   36.5 |                              26.6   |                       26.4   |
| LS+SLU (single)       |                  642   |                               0.775 |                        0.607 |
| LS+SLU                |                  588   |                               0.932 |                        0.769 |
| LS+KLU (single)       |                  945   |                               0.277 |                        0.116 |
| LS+KLU                |                  918   |                               0.306 |                        0.144 |
| LS+NICSLU (single)    |                  947   |                               0.274 |                        0.11  |
| LS+NICSLU             |                  929   |                               0.298 |                        0.134 |

For more information (including the exact way to reproduce these results, as well as the computer used), you can consult the dedicated [Benchmarks](https://lightsim2grid.readthedocs.io/en/latest/benchmarks.html) page on the documentation.

## Philosophy
Lightsim2grid aims at providing a somewhat efficient (in terms of computation speed) backend targeting the 
grid2op platform.

It provides a c++ api, compatible with grid2op that is able to compute flows (and voltages and reactive power) from
a given grid. This grid can be modified according to grid2op mechanism (see more information in the [official
grid2Op documentation](https://grid2op.readthedocs.io/en/latest/index.html) ).

This code do not aim at providing state of the art solver in term of performances nor in terms of realism in the
modeling of power system elements (*eg* loads, generators, powerlines, transformers, etc.).

Lightsim2grid codebase is "organized" in 4 different parts:

1) modify the elements (*eg* disconnecting a powerline or changing the voltage magnitude setpoint of a 
   generator, or any other action made possible by grid2op)
2) generate the `Ybus` (sparse) complex admitance matrix and `Sbus` complex injection vector from the state of the
   powergrid (*eg* physical properties of each elements, which elements are in service, which power is produce at 
   each generators and consumed at each loads, what is the grid topology etc.)
3) solving for the complex voltage `V` (and part of the `Sbus` vector) the equation `V.(Ybus.V)* = Sbus` with the 
   "standard" "powerflow constraints"
   (*eg* the voltage magnitude of `V` is set at given components, and on other it's the imaginary part of `Sbus`)
4) computes the active power, reactive power, flow on powerllines etc. from the `V` and `Sbus` complex vectors computed
   at step 3).

Step 1, 2 and 4 are done in the [GridModel](https://lightsim2grid.readthedocs.io/en/latest/gridmodel.html#lightsim2grid.gridmodel.GridModel) class.

Step 3 is performed thanks to a "powerflow solver".

### Using a custom powerflow solver
For now some basic "solver" (*eg* the program that performs points `3.` above) are available, based on the
Gauss Seidel or the Newton-Raphson methods to perform "powerflows". 

Nothing prevents any other "solver" to be used with lightsim2grid and thus with grid2op. For this, you simply need to
implement, in c++ a "lightsim2grid solver" which mainly consists in defining a function:
```c
bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,  // the admittance matrix
                CplxVect & V,  // store the results of the powerflow and the Vinit !
                const CplxVect & Sbus,  // the injection vector
                const Eigen::VectorXi & ref,  // bus id participating to the distributed slack
                const RealVect & slack_weights,  // slack weights for each bus
                const Eigen::VectorXi & pv,  // (might be ignored) index of the components of Sbus should be computed
                const Eigen::VectorXi & pq,  // (might be ignored) index of the components of |V| should be computed
                int max_iter,  // maximum number of iteration (might be ignored)
                real_type tol  // solver tolerance 
                );
```

The types used are:

- `real_type`: double => type representing the real number
- `cplx_type` :  std::complex<real_type> => type representing the complex number
- `CplxVect` : Eigen::Matrix<cplx_type, Eigen::Dynamic, 1> => type representing a vector of complex numbers
- `RealVect` : Eigen::Matrix<real_type, Eigen::Dynamic, 1> => type representing a vector of real numbers
- `Eigen::VectorXi` => represents a vector of integer
- `Eigen::SparseMatrix<cplx_type>` => represents a sparse matrix

See for example [BaseNRSolver](./src/BaseNRSolver.h) for the implementation of a Newton Raphson solver (it requires some "linear solvers", more details about that are given in the section bellow)

Any contribution in this area is more than welcome.

**NB** For now the "solver" only uses these above information to perform the powerflow. If a more
"in depth" solution needs to be implemented, let us know with a github issue. For example, it could be totally fine that a proposed "solver" uses direct information about the elements (powerline, topology etc.) of the grid in order to perform some powerflow.

**NB** It is not mandatory to "embed" all the code of the solver in lightsim2grid. Thanks to different customization, 
it is perfectly possible to install a given "lightsim solver" only if certain conditions are met. For example, on
windows based machine, the SuiteSparse library cannot be easily compiled, and the KLUSolver is then not available.

**NB** It would be totally fine if some "lightsim2grid" solvers are available only if some packages are installed on the
machine for example.

### Using custom linear solvers to solve powerflows

In lightsim2grid (c++ part) it is also possible, thanks to the use of "template meta programming" to
not recode the Newton Raphson algorithm (or the DC powerflow algorithm) and to leverage the 
use of a linear solver.

A "linear solver" is anything that can implement 3 basic functions:

- `initialize(const Eigen::SparseMatrix<real_type> & J)` : initialize the solver and prepare it to solve for linear systems `J.x = b` (usually called once per powerflow)
- `ErrorType solve(const Eigen::SparseMatrix<real_type> & J, RealVect & b, bool has_just_been_inialized)`: effectively solves `J.x = b` (usually called multiple times per powerflow)
- `ErrorType reset()`: clear the state of the solver (usually performed at the end of a powerflow
  to reset the state to a "blank" / "as if it was just initialized" state)

Some example are given in the c++ code "KLUSolver.h", "SparLUSolver.h" and "NICSLU.h"

This usage usually takes approximately around 20 / 30 lines of c++ code (not counting the comments, and boiler code for exception handling for example).

## Citing

If you use this package in one of your work, please cite:
```
@misc{lightsim2grid,
    author = {B. Donnot},
    title = {{Lightsim2grid - A c++ backend targeting the Grid2Op platform. }},
    year = {2020},
    publisher = {GitHub},
    journal = {GitHub repository},
    howpublished = {\url{https://GitHub.com/bdonnot/lightsim2grid}},
}
```

## Miscellaneous

### Customization of the compilation
#### Enable NICSLU
For that, you need to declare the environment variables `PATH_NICSLU` that points to a valid installation of
the NICSLU package (see https://github.com/chenxm1986/nicslu). 
For example: `export PATH_NICSLU=/home/user/Documents/nicslu/nicslu202103`

#### Enable 03 optimization
By default, at least on ubuntu, only the "-O2" compiler flags is used. To use the O3 optimization flag, you need
to specify the `__COMPLILE_O3` environment variable: `set __COMPLILE_O3=1` before the compilation (so before
`python3 setup.py build` or `python -m pip install -e .`)

This compilation argument will increase the compilation time, but will make the package faster.

#### Enable "-march=native" optimization
By default, for portability, we do not compile with `-march=native` flags. This lead to some error on some platform.
If you want to further improve the performances.

You can `set __COMPILE_MARCHNATIVE=1` to enable it before the compilation (so before
`python3 setup.py build` or `python -m pip install -e .`)

### Profile the code
This is a work in progress for now. And it is far from perfect, and probably only work on linux.

See https://github.com/xflash96/pybind11_package_example/blob/main/tutorial.md#perf for more details.

```commandline
cd benchmarks
perf record ./test_profile.py
perf report
```

### Local testing
And some official tests, to make sure the solver returns the same results as pandapower
are performed in "lightsim2grid/tests"
```bash
cd lightsim2grid/tests
python -m unittest discover
```

This tests ensure that the results given by this simulator are consistent with the one given by pandapower when
using the Newton-Raphson algorithm, with a single slack bus, without enforcing q limits on the generators etc.

**NB** to run these tests you need to install grid2op from source otherwise all the test of the LightSim2gridBackend 
will fail. In order to do so you can do:
```
git clone https://github.com/rte-france/Grid2Op.git
cd Grid2Op
pip3 install -U -e .
cd ..
```
### Tests performed automatically

Some tests are performed automatically on standard platform each time modifications are made in the lightsim2grid code.

These tests include, for now, compilation on gcc (version 8, 10, 11 and 12) and clang (version 11, 13 and 14).

**NB** Intermediate versions of clang and gcc (*eg* gcc 9 or clang 12) are not tested regularly, but lightsim2grid used to work on these. We suppose that if it works on *eg* clang 10 and clang 14 then it compiles also on all intermediate versions.

**NB** Package might work (we never tested it) on earlier version of these compilers. 
The only "real" requirement for lightsim2grid is to have a compiler supporting c++11
(at least).

### Known issues

#### Storage units
There are discrepency in the handling of storage units, when the are not asked to produce / consume anything (setpoint
is 0.) between pandapower and lightsim2grid only in the case where the storage unit is alone on its bus.

Pandapower does not detect it and the episode can continue. On the other side, lightsim2grid detects it and raise an
error because in that case the grid is not connex anymore (which is the desired behaviour).

#### Compilation issue

On the clang compiler (default one on MacOS computer) it is sometime require to downgrade the pybind11 version
to 2.6.2 to install the package.

You can downgrade pybind11 with: `python -m pip install -U pybind11==2.6.2`
