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
- 3.11 (lightsim2grid >= 0.7.1)
- 3.12 (lightsim2grid >= 0.7.5)

**NB** on some version of MacOs (thanks Apple !), especially the one using M1 or M2 chip, lightsim2grid is only available
on pypi starting from version 0.7.3

## Installation (from source, for more advanced user)

See the official documentation at [Install from source](https://lightsim2grid.readthedocs.io/en/latest/install_from_source.html) for more information

## Benchmarks

Lightsim2grid is significantly faster than pandapower when used with grid2op for all kind of environment size (sometimes 
more than 30x faster - making 30 steps while pandapower makes one).


If you prefer to use the dedicated lightsim2grid `SecurityAnalysis` or `TimeSerie` classes you can even expect another 10-20x
speed ups compared to grid2op with lightsim2grid (sometimes more than 300x faster than grid2op with pandapower). 

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

#### Enable CKTSO
For that, you need to declare the environment variables `PATH_CKTSO` that points to a valid installation of
the NICSLU package (see https://github.com/chenxm1986/cktso). 
For example: `export PATH_NICSLU=/home/user/Documents/cktso`

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

These tests include, for now, compilation on gcc (version 8, 12 and 13) and clang (version 11, 16 and 17).

**NB** Intermediate versions of clang and gcc (*eg* gcc 9 or clang 12) are not tested regularly, but lightsim2grid used to work on these. 
We suppose that if it works on *eg* clang 10 and clang 14 then it compiles also on all intermediate versions.

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
