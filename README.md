# LightSim2Grid
Provide a fast backend for grid2op using c++ KLU and Eigen librairies. Its primary goal is to serve as a fast
backend for the grid2op platform, used primarily as a testbed platform for sequential decision making in
the world of power system.

See the [Disclaimer](DISCLAIMER.md) to have a more detailed view on what is and what is not this package. For example
this package should not be used for detailed power system computations or simulations.

## Usage
Once installed (don't forget, if you used the optional virtual env
above you need to load it with `source venv/bin/activate`) you can
use it as any python package.

### 1. As a grid2op backend (preferred method)
This functionality requires you to have grid2op installed, with at least version 0.7.0. You can install it with
```commandline
pip install grid2op>=0.7.0
```

Then you can use a LightSimBackend instead of the default PandapowerBackend this way:

```python3
import grid2op
from lightsim2grid import LightSimBackend
backend = LightSimBackend()
env = grid2op.make(backend=backend)
# do regular computation as you would with grid2op
```
And you are good to go.

### 2. replacement of pandapower "newtonpf" method (advanced method)
It is also possible to use directly the "solver" part of lightsim2grid.

Suppose you somehow get:
- `Ybus` the admittance matrix of your powersystem given by pandapower
- `V0` the (complex) voltage vector at each bus given by pandapower
- `Sbus` the (complex) power absorb at each bus as given by pandapower
- `ppci` a ppc internal pandapower test case
- `pv` list of PV buses
- `pq` list of PQ buses
- `options` list of pandapower "options"

You can define replace the `newtonpf` function of `pandapower.pandapower.newtonpf` function with the following
piece of code:
```python
from lighsim2grid.newtonpf import newtonpf
V, converged, iterations, J = newtonpf(Ybus, V, Sbus, pv, pq, ppci, options)
```

This function uses the KLU algorithm and a c++ implementation of a Newton solver for speed.

## Installation (from source, recommended)
You need to:
- clone this repository and get the code of Eigen and SparseSuite (mandatory for compilation)
- (optional) compile a piece of SparseSuite
- install the package

### Important note
This package relies on the excellent pybind11 package to integrate c++ code into python easily. 

So to install lightsim2grid you need pybind and its requirement, which include a working compiler: for example 
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

### (optional) Compilation of SuiteSparse
SuiteSparse comes with the faster KLU linear solver. Since version 0.3.0 this requirement has been removed. This entails
that on linux / macos you can still benefit from the faster KLU solver. On windows you will still benefit from the
speed up of lightsim (versus the default PandaPowerBackend) but this speed up will be less than if you manage
to compile SuiteSparse (see the subsection [Benchmark](#benchmark) for more information).

**NB** in both cases the algorithm to compute the powerflow is exactly the same. It is a 
Newton Raphson based method. But to carry out this algorithm, one need to solver some linear equations. The only
difference in the two version (with KLU and without) is that the linear equation solver is different. Up to the
double float precision, both results (with and without KLU) should match.

We only detail the compilation on a system using "make" (so most likely GNU-Linux and MacOS). If you manage to
do this step on Windows, you can continue (and let us know!). If you don't feel comfortable with this, we
provided a docker version. See the next section [Installation Using Docker](#installation-using-docker) 
for more information.

```commandline
# compile static libraries of SparseSuite
make
```

### 2. Installation of the python package
Now you simply need to install the lightsim2grid package this way, like any python package:

```commandline
# install the dependency
pip install -U pybind11
# compile and install the python package
pip install -U .
```

And you are done :-)

### Benchmark
In this section we will expose some brief benchmarks about the use of lightsim2grid in the grid2op settings.
The code to run these benchmarks are given with this package int the [benchmark](./benchmarks) folder.

All of them has been run on a computer with a `Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz` processor. The command
to run the benchmark is (once `cd` in the [benchmark](./benchmarks) folder):
```commandline
cd benchmarks  # cd in the lightsim2grid benchmarks directory if not already
python3 benchmark_solvers.py --name l2rpn_case14_sandbox --no_test --number 1000
python3 benchmark_solvers.py --name l2rpn_neurips_2020_track2_small --no_test --number 1000
```
(results may vary depending on the hard drive, the ram etc. and are presented here for illustration only)

(we remind that these simulations correspond to simulation on one core of the CPU. Of course it is possible to
make use of all the available cores, which would increase the number of steps that can be performed)

Results were obtained with grid2op version 1.5.0.


We compare 5 different solvers:

- **PP**: PandaPowerBackend (default grid2op backend) which is the reference in our benchmarks (uses the numba
  acceleration)
- **LS+GS** (LightSimBackend+Gauss Seidel): the grid2op backend based on lightsim2grid that uses the "Gauss Seidel"
  solver to compute the powerflows It is implemented in
  [GaussSeidelSolver](./src/GaussSeidelSolver.h).
- **LS+GS S** (LightSimBackend+Gauss Seidel Synchronous): the grid2op backend based on lightsim2grid that uses a
  variant of the "Gauss Seidel" method to compute the powerflows. It is implemented in
  [GaussSeidelSynchSolver](./src/GaussSeidelSynchSolver.h).
- **LS+SLU** (LightSimBackend+SparseLU): the grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver "SparseLU" from the
  Eigen c++ library (available on all platform) and is implemented in
  [SparseLUSolver](./src/SparseLUSolver.h).
- **LS+KLU** (LightSimBackend+KLU): he grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "KLU" from the SuiteSparse c package, available only (for now) on Linux and Mac OS and is implemented in
  [KLUSolver](./src/KLUSolver.h).

First on an environment based on the IEEE case14 grid:

| case14_sandbox   |   grid2op speed (it/s) |   grid2op powerflow time (ms) |   solver powerflow time (ms) |
|------------------|------------------------|-------------------------------|------------------------------|
| PP               |                     67 |                        11.5   |                       4.38   |
| LS+GS            |                    779 |                         0.481 |                       0.356  |
| LS+GS S          |                    759 |                         0.51  |                       0.383  |
| LS+SLU           |                    988 |                         0.205 |                       0.0771 |
| LS+KLU           |                   1018 |                         0.177 |                       0.0473 |

From a grid2op perspective, lightsim2grid allows to compute up to ~1000 steps each second on the case 14 and
"only" 67 for the default PandaPower Backend, leading to a speed up of **~17** in this case
(lightsim2grid is ~15 times faster than Pandapower). For such a small environment, there is no difference in using 
KLU linear solver (not available on
Windows based machine) compared to using the SparseLU solver of Eigen (1018 vs 988 iterations on the reported
runs, might slightly vary across runs)

Then on an environment based on the IEEE case 118:

| neurips_2020_track2   |   grid2op speed (it/s) |   grid2op powerflow time (ms) |   solver powerflow time (ms) |
|-----------------------|------------------------|-------------------------------|------------------------------|
| PP                    |                     38 |                        13.9   |                        5.71  |
| LS+GS                 |                      4 |                       201     |                      201     |
| LS+GS S               |                     29 |                        32.3   |                       32     |
| LS+SLU                |                    546 |                         0.932 |                        0.645 |
| LS+KLU                |                    673 |                         0.586 |                        0.298 |

For an environment based on the IEEE 118, the speed up in using lightsim + KLU (LS+KLU)
[for now only available on linux and MacOS] is **~18** time faster than
using the default PandaPower backend. The speed up of lightsim + SparseLU is a bit lower, but it is still **~14**
times faster than using the default backend [the `LS+KLU` solver is ~2-3 times faster than the `LS+SLU` solver 
(`0.298` ms per powerflow for `L2+KLU`  compared to `0.645` ms for `LS+SLU`), but it only translates to `LS+KLU` 
providing ~20-25% more
iterations per second in the total program (`673` vs `546`) mainly because grid2op itself takes some times to modify the
grid and performs all the check it does.]

If we look now only at the time to compute one powerflow (and don't take into account the time to load the data, to
initialize the solver, to modify the grid, read back the results, to perform the other update in the
grid2op environment etc.) we can notice that it takes on average (over 1000 different states) approximately **0.30ms**
to compute a powerflow with the LightSimBackend (if using the KLU linear solver) compared to the **5.7 ms** when using
the PandaPowerBackend (speed up of **~20** times)

**NB** pandapower performances heavily depends on the pandas version used, we used here pandas version 1.1.4 which
we found gave the best performances on our machine.

## Philosophy
Lightsim2grid aims at providing a somewhat efficient (in terms of computation speed) backend targeting the 
grid2op platform.

It provides a c++ api, compatible with grid2op that is able to compute flows (and voltages and reactive power) from
a given grid. This grid can be modified according to grid2op mechanism (see more information in the [official
grid2Op documentation](https://grid2op.readthedocs.io/en/latest/index.html`) ).

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

For now some basic "solver" (*eg* the program that performs points `3.` above) are available, based on the
Gauss Seidel or the Newton-Raphson methods to perform "powerflows". 

Nothing prevents any other "solver" to be used with lightsim2grid and thus with grid2op. For this, you simply need to
implement, in c++ a "lightsim2grid solver" which mainly consists in defining a function:
```c
bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,  // the admittance matrix
                CplxVect & V,  // store the results of the powerflow and the Vinit !
                const CplxVect & Sbus,  // the injection vector
                const Eigen::VectorXi & pv,  // (might be ignored) index of the components of Sbus should be computed
                const Eigen::VectorXi & pq,  // (might be ignored) index of the components of |V| should be computed
                int max_iter,  // maximum number of iteration (might be ignored)
                real_type tol  // solver tolerance 
                );
```

The types used are:

- `real_type`: double => type representing the real number
- `cplx_type` :  std::complex<real_type> => type representing the complex number
- `CplxVect` : Eigen::Matrix<cplx_type, Eigen::Dynamic, 1> => type representing a vector of complex elements
- `Eigen::VectorXi` => represents a vector of integer
- `Eigen::SparseMatrix<cplx_type>` => represents a sparse matrix

See for example [BaseNRSolver](./src/BaseNRSolver.h) for the implementation of a Newton Raphson solver (and
its derived classes [KLUSolver](./src/KLUSolver.h) and [SparseLUSolver](./src/SparseLUSolver.h) that uses
different routine to implement this algorithm) for examples on how to implement a solver.

Any contribution in this area is more than welcome.

**NB** It is not mandatory to "embed" all the code of the solver in lightsim2grid. Thanks to different customization, 
it is perfectly possible to install a given "lightsim solver" only if certain conditions are met. For example, on
windows based machine, the SuiteSparse library cannot be easily compiled, and the KLUSolver is then not available.

It would be totally fine if some "lightsim2grid" solvers are available only if some packages are installed on the
machine for example.

## Installation (using docker)
In this section we cover the use of docker with grid2op.

### 1. Install docker
First, you need to install docker. You can consult the 
[docker on windows](https://hub.docker.com/editions/community/docker-ce-desktop-windows) if you use a windows like
operating system, if you are using MacOs you can consult 
[docker on Mac](https://hub.docker.com/editions/community/docker-ce-desktop-mac/). The installation of docker on linux
depends on your linux distribution, we will not list them all here.

### 2. Get the lightsim2grid image

Once done, you can simply "install" the lightsim2grid image with:
```commandline
docker pull bdonnot/lightsim2grid:latest
```

This step should be done only once (unless you delete the image) it will download approximately 4 or 5GB from the
internet. The lightsim2grid image contains lightsim and grid2op python packages (as well as their
dependencies), equivalent of what would be installed if you typed:
```commandline
pip install -U grid2op[optional] pybind11
# and do steps detailed in section "Installation (from source)"
# that we will not repeat
```

### 3. Run a code on this container
You can skip this section if you know how to use docker. We will present here "the simplest way" to use. This is NOT
a tutorial on docker, and you can find better use of this technology on 
[the docker website](https://www.docker.com/get-started).

For this tutorial, we suppose you have a script named `my_script.py` located in the directory (complete path) 
`DIR_PATH` (*e.g.* on windows you can have `DIR_PATH` looking like "c:\User\MyName\L2RPNCompeitionCode" or 
on Linux `DIR_PATH` will look like "/home/MyName/L2RPNCompeitionCode", this path is your choice, you can name it
the way you like)


#### 3.1) Start a docker container
You first need to start a docker container and tell docker that the container can access your local files with the 
following command:

```commandline
docker run -t -d -p 8888:8888 --name lightsim_container -v DIR_PATH:/L2RPNCompeitionCode -w /L2RPNCompeitionCode bdonnot/lightsim2grid
```
More information on this command 
[in the official docker documentation](https://docs.docker.com/engine/reference/commandline/run/)

After this call you can check everything went smoothly with by invoking:
```commandline
docker ps
```
And the results should look like:
```
CONTAINER ID        IMAGE                   COMMAND             CREATED             STATUS              PORTS               NAMES
89750964ca55        bdonnot/lightsim2grid   "python3"           5 seconds ago       Up 4 seconds        80/tcp              lightsim_container
```

**NB** I insist, `DIR_PATH` should be replaced by the path on which you are working, see again the introduction of this
section for more information, in the example above this can look like:
```commandline
docker run -t -d -p 8888:8888 --name lightsim_container -v /home/MyName/L2RPNCompeitionCode:/L2RPNCompeitionCode -w /L2RPNCompeitionCode bdonnot/lightsim2grid
```

#### 3.2) Execute your code on this container
Once everything is set-up you can execute anything you want on this container. Note that doing so, the execution
of the code will be totally independant of your system. Only the things located in `DIR_PATH` will be visible 
by your script, only the python package installed in the container will be usable, only the python interpreter
of the containter (python 3.6 at time of writing) will be usable etc.

```commandline
docker exec lightsim_container python my_script.py
```

Of course, the "my_script.py" should save its output somewhere on the hard drive.

If you rather want to execute a python REPL (read-eval-print loop), corresponding to the "interactive python 
interpreter", you can run this command:
```commandline
docker exec -it lightsim_container python
```

We also added the possibility to run jupyter notebook from this container. To do so, you can run the command:
```commandline
docker exec -it lightsim_container jupyter notebook --port=8888 --no-browser --ip='*' --allow-root
```

More information is provided in the official documentation of 
[docker exec](https://docs.docker.com/engine/reference/commandline/exec/).

#### 3.3) Disclaimer
Usually, docker run as root on your machine, be careful, you can do irreversible things with it. "A great power 
comes with a great responsibility".

Also, we recall that we presented a really short introduction to docker and its possibility. We have not implied
that this was enough, nor explain (on purpose, to make this short) any of the commands. 
We strongly encourage you to have a look for yourself. 

We want to recall the paragraph `7. Limitation of Liability` under which lightsim2grid, and this "tutorial" 
is distributed:

*Under no circumstances and under no legal 
theory, whether tort (including negligence), 
contract, or otherwise, shall any Contributor, or 
anyone who distributes Covered Software as 
permitted above, be liable to You for any direct, 
indirect, special, incidental, or consequential 
damages of any character including, without 
limitation, damages for lost profits, loss of 
goodwill, work stoppage, __**computer failure or**__
__**malfunction**__, or any and all other commercial 
damages or losses, even if such party shall have 
been informed of the possibility of such damages.*

### 4) Clean-up
Once you are done with your experiments, you can stop the docker container:
```commandline
docker container stop lightsim_container
```
This will free all the CPU / GPU resources that this container will use. If you want to start it again, for another 
experiment for example, just use the command:
```commandline
docker container start lightsim_container
```
This will allow you to run another batch of `dcoker exec` (see `3.2) Execute your code on this container`) 
without having to re run the container.


If you want to go a step further, you can also delete the container with the command:
```commandline
docker container rm lightsim_container
```
This will remove the container, and all your code executed there, the history of commands etc. If you want to use
lightsim2grid with docker again you will have to go through section `3. Run a code on this container` all over
again.

And if you also want to remove the image, you can do:
```commandline
docker rmi bdonnot/lightsim2grid 
```
**NB** this last command will completely erase the lightsim2grid image from your machine. This means that 
if you want to use it again, you will have to download it again (see section `2. Get the lightsim2grid image`)

Finally, you can see the official documentation in case you need to uninstall docker completely from your system.

## Miscellaneous

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

These tests include, for now, compilation on gcc (version 7, 8, 9 and 10) and clang (version 7, 8, 9 and 10).