# LightSim2Grid
Provide a fast backend for grid2op using c++ KLU and Eigen librairies. Its primary goal is to serve as a fast
backend for the grid2op platform, used primarily as a testbed platform for sequential decision making in
the world of power system.

See the [Disclaimer](DISCLAIMER.md) to have a more detailed view on what is and what is not this package. For example
this package should not be used for detailed power system computations or simulations.

*   [1 Key features](#Key-features)
*   [2 Usage](#Usage)
    *   [2.1. As a grid2op backend (preferred method)](#1-as-a-grid2op-backend-preferred-method)
    *   [2.2. replacement of pandapower "newtonpf" method (advanced method)](#2-replacement-of-pandapower-newtonpf-method-advanced-method)
*   [3 Installation (from pypi official repository, recommended)](#Installation-from-pypi-official-repository-recommended)
*   [4 Installation (from source, for more advanced user)](#Installation-from-source-for-more-advanced-user)
*   [5. Benchmarks](#Benchmarks)
*   [6. Philosophy](#Philosophy)
*   [7. Miscellaneous](#Miscellaneous)
    * [7.1 Customization of the compilation](#Customization-of-the-compilation)
    * [7.2 Profile the code](#Profile-the-code)
    * [7.3 Local testing](#Local-testing)
    * [7.4 Tests performed automatically](#Tests-performed-automatically)
    * [7.5 Known issues](#Known-issues)

## Key features

- **Multiple powerflow algorithms**: Newton-Raphson (single or distributed slack), Fast-Decoupled and
  Gauss-Seidel, each selectable with a choice of linear solver (Eigen SparseLU, SuiteSparse KLU, NICSLU
  or CKTSO) via the `lightsim2grid.algorithm` module -- see the
  [Available powerflow algorithms](https://lightsim2grid.readthedocs.io/en/latest/solvers.html) and
  [Naming conventions](https://lightsim2grid.readthedocs.io/en/latest/algorithm_names.html) pages
  (`lightsim2grid.solver` is the older, deprecated name for this same module).
- **Load your own solver as a plugin** at runtime, without forking or recompiling lightsim2grid -- see
  [External Algorithm Plugins](https://lightsim2grid.readthedocs.io/en/latest/solver_plugin.html).
- **Multiple grid source formats**: pandapower, pypowsybl (iidm), MATPOWER and PowerModels.jl, in
  addition to grid2op's own environments -- see
  [LSGrid module](https://lightsim2grid.readthedocs.io/en/latest/network.html) (`lightsim2grid.gridmodel`
  / `GridModel` are the older, deprecated names for `lightsim2grid.network` / `LSGrid`).
- **Multi-threaded contingency (n-1) analysis**, time-series, injection-sweep and scenario-sweep
  computation, all much faster than stepping through grid2op one scenario / contingency at a time --
  see [Contingency Analysis](https://lightsim2grid.readthedocs.io/en/latest/security_analysis.html) and
  [Time series](https://lightsim2grid.readthedocs.io/en/latest/time_series.html). Every one of these
  batch classes exposes `nb_solved()` / `nb_converged()` and a per-row `converged_mask()`, so a failing
  row is diagnosable without re-running the batch one row at a time.
- **Fast binary serialization** (`save_binary` / `load_binary`) of a grid, in addition to the standard
  pickle support -- see
  [Fast binary serialization](https://lightsim2grid.readthedocs.io/en/latest/binary_serialization.html).
- **PTDF / LODF matrices** for near-instant linear (DC) sensitivity analysis -- see the
  [LSGrid module](https://lightsim2grid.readthedocs.io/en/latest/network.html#ptdf-lodf-section) page.
- A standalone `lightsim2grid_core` **C++ library** (no python required) for embedding the solver in a
  pure C++ project -- see the
  [C++ library](https://lightsim2grid.readthedocs.io/en/latest/cpp_library.html) page.

## Usage
Once installed (don't forget, if you used the optional virtual env
above you need to load it with `source venv/bin/activate`) you can
use it as any python package.

### 1. As a grid2op backend (preferred method)
This functionality requires you to have grid2op installed, with at least version 1.6.4. You can install it with

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

# when pandapower version <= 2.7.0
# V, converged, iterations, J, Vm_it, Va_it = newtonpf(Ybus, Sbus, V0, pv, pq, ppci, options)

# when pandapower version > 2.7.0
V, converged, iterations, J, Vm_it, Va_it = newtonpf(Ybus, Sbus, V0, ref, pv, pq, ppci, options)
```

This function uses the KLU algorithm (when available) and a c++ implementation of a Newton solver for speed.

`lightsim2grid.newtonpf` is a re-export of `lightsim2grid.pandapower_compat.newtonpf`, which also provides a
DC counterpart, `dcpf` (`from lightsim2grid.pandapower_compat import dcpf`). See the ["Use as Pandapower
Solver"](https://lightsim2grid.readthedocs.io/en/latest/use_solver.html) page of the documentation for more.

## Installation (from pypi official repository, recommended)

Since version 0.5.3, lightsim2grid is can be installed like most python packages, with a call to:
`python -m pip install lightsim2grid`

It includes faster grid2op backend and the `SuiteSparse` faster `KLU` solver, even on windows. This is definitely the 
easiest method to install lightsim2grid on your system and have it running.

Note though that these packages have been compiled on a different platform that the one you are using. You might still
get some benefit (in terms of performances) to install it from your on your machine with the proper compilations flags (
see section [6.1 Customization of the compilation](#Customization-of-the-compilation) for more information)

Pypi packages are available for linux (`x86_64` cpu architecture), windows (`x86_64` cpu architecture) and macos (`x86_64` cpu architecture) with python versions: 

- 3.7  (lightsim2grid < 0.10.4)
- 3.8
- 3.9
- 3.10 (lightsim2grid >= 0.6.1)
- 3.11 (lightsim2grid >= 0.7.1)
- 3.12 (lightsim2grid >= 0.7.5)
- 3.13 (lightsim2grid >= 0.9.2.post2)
- 3.14 (lightsim2grid >= 0.10.4)


As from version 0.8.2, we also distribute windows `arm64` and macos `arm64` binaries of lightsim2grid that can be installed
directly with pip too (requires python >= 3.8 for macos and python >= 3.9 for windows). 

We do not currently produce `arm64` (`aarch64`) linux binaries because it takes too long to build. If you really want them, let us know and we'll see what we can do.


**NB** on some version of MacOs (thanks Apple !), especially the one using M1 or M2 chip, lightsim2grid is only available
on pypi starting from version 0.7.3 We attempted to deliver arm64 lightsim2grid version but we could not test them. So if you 
want a reliable working and tested version of lightsim2grid on newest version of macos (with M1 or M2 chips for example) please
use lightsim2grid >= 0.8.2

**NB** we do not currently build any 32 bits lightsim2grid libraries.

## Installation (from source, for more advanced user)

See the official documentation at [Install from source](https://lightsim2grid.readthedocs.io/en/latest/install_from_source.html) for more information

## Benchmarks

Lightsim2grid is significantly faster than pandapower when used with grid2op for all kind of environment size (sometimes 
more than 30x faster - making 30 steps while pandapower makes one).


If you prefer to use the dedicated lightsim2grid `ContingencyAnalysis` or `TimeSerie` classes you can even expect another 10-20x
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

Step 1, 2 and 4 are done in the [LSGrid](https://lightsim2grid.readthedocs.io/en/latest/network.html#lightsim2grid.network.LSGrid) class.

Step 3 is performed thanks to a "powerflow solver".

### Using a custom powerflow solver
Several "solvers" (*eg* the program that performs point `3.` above) are available out of the box, based on
Gauss-Seidel, Newton-Raphson or Fast-Decoupled Power Flow methods, each combined with a choice of linear
solver (Eigen's SparseLU, SuiteSparse KLU, NICSLU or CKTSO) -- see the
[Available powerflow algorithms](https://lightsim2grid.readthedocs.io/en/latest/solvers.html) page.

Beyond that, lightsim2grid also supports loading your **own** solver at runtime as a plugin, without
having to fork or recompile lightsim2grid itself: implement a C++ class deriving from `BaseAlgo`, build it
as a separate shared library against an installed lightsim2grid, then load it with
`lightsim2grid.load_algorithm_plugin(...)` and select it with `grid.change_algorithm("YourSolverName")`.
See the [External Algorithm Plugins](https://lightsim2grid.readthedocs.io/en/latest/solver_plugin.html)
page of the documentation for the full mechanism, worked examples (`examples/external_algorithm/`,
`examples/dist_slack_algorithm/`, `examples/lm_algorithm/`), and the build-flag pitfalls to avoid.

### Using custom linear solvers to solve powerflows

In lightsim2grid (c++ part) it is also possible, thanks to the use of "template meta programming" to
not recode the Newton Raphson algorithm (or the DC powerflow algorithm) and to leverage the 
use of a linear solver.

A "linear solver" is anything that can implement these basic functions:

- `ErrorType analyze(const EigenRefConstRealSpMat & J)`: reordering + symbolic factorization of `J` (structure only, usually called once per powerflow)
- `ErrorType factorize(const EigenRefConstRealSpMat & J)`: numeric factorization of `J` (requires values, usually called once per powerflow)
- `ErrorType refactorize(const EigenRefConstRealSpMat & J)`: re-numeric factorization, reuses the symbolic factorization from `analyze` (usually called multiple times per powerflow, when only the values of `J` changed)
- `ErrorType solve(Eigen::Ref<RealVect> b) const`: effectively solves `J.x = b` in place, given the previous `analyze` / `factorize` (or `refactorize`) call
- `ErrorType reset()`: clear the state of the solver (usually performed at the end of a powerflow
  to reset the state to a "blank" / "as if it was just initialized" state)

Some example are given in the c++ code "KLUSolver.hpp", "SparseLUSolver.hpp" and "NICSLUSolver.hpp" (in `src/core/linear_solvers`)

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
    howpublished = {\url{https://GitHub.com/Grid2Op/lightsim2grid}},
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
the CKTSO package (see https://github.com/chenxm1986/cktso). 
For example: `export PATH_CKTSO=/home/user/Documents/cktso`

#### Enable O3 optimization
The "-O3" compiler flag is now used by default. To disable it (and fall back to "-O2"), you need
to specify the `__O3_OPTIM` environment variable before compiling (`export __O3_OPTIM=0` on linux / macos,
`set __O3_OPTIM=0` on windows cmd, or `$Env:__O3_OPTIM=0` in powershell), then `python -m pip install -e .`

This compilation argument will increase the compilation time, but will make the package faster.

#### Enable "-march=native" optimization
By default, for portability, we do not compile with `-march=native` flags. This lead to some error on some platform.
If you want to further improve the performances.

You can set the `__COMPILE_MARCHNATIVE` environment variable to `1` to enable it before compiling
(`export __COMPILE_MARCHNATIVE=1` on linux / macos, `set __COMPILE_MARCHNATIVE=1` on windows cmd, or
`$Env:__COMPILE_MARCHNATIVE=1` in powershell), then `python -m pip install -e .`

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
git clone https://github.com/Grid2Op/grid2op.git
cd grid2op
pip3 install -U -e .
cd ..
```
### Tests performed automatically

Some tests are performed automatically on standard platform each time modifications are made in the lightsim2grid code.

These tests include, for now, compilation on gcc (version 8, 14 and 15) and clang (version 11, 20 and 21).

**NB** Intermediate versions of clang and gcc (*eg* gcc 9 or clang 12) are not tested regularly, but lightsim2grid used to work on these. 
We suppose that if it works on *eg* clang 11 and clang 21 then it compiles also on all intermediate versions.

**NB** Package might work (we never tested it) on earlier version of these compilers. 
The only "real" requirement for lightsim2grid is to have a compiler supporting c++14
(at least).

CI also explicitly compiles lightsim2grid (both the C++ unit test suite and the python
bindings) pinned to C++14 and to C++26: the oldest and the latest standard claimed to be
supported (see the `LS2G_CXX_STANDARD` CMake option and
`.github/workflows/cpp-standards.yml`). Absent that option, the build auto-detects and
uses the newest standard the compiler supports, from C++26 down to C++14.

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
