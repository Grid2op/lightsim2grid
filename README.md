# LightSim2Grid
Provide a fast backend for grid2op using c++ KLU and Eigen librairies. Its primary goal is to serve as a fast
backend for the grid2op platform, used primarily as a testbed platform for sequential decision making in
the world of power system.

See the [Disclaimer](DISCLAIMER.md) to have a more detailed view on what is and what is not this package. For example
this package should not be used for detailed power system computations or simulations.

## installation
You need to:
- clone this repository
- get the code of Eigen and SparseSuite (mandatory for compilation)
- compile a piece of SparseSuite
- install the package

All of that can be done, on a computer with make (so most likely MacOS or GNU/Linux based):
```bash
git clone https://github.com/BDonnot/lightsim2grid.git
cd lightsim2grid
# it is recommended to do a python virtual environment
python -m virtualenv venv  # optional
source venv/bin/activate  # optional

# retrieve the code of SparseSuite and Eigen
git submodule init
git submodule update

# compile static libraries of SparseSuite
make

# compile and install the python package
pip install -U pybind11
pip install -U .
```

## Usage
Once installed (don't forget, if you used the optional virtual env
above you need to load it with `source venv/bin/activate`) you can
use it as any python package.

### 1. As a grid2op backend (preferred method)
This functionality requires you to have grid2op installed, with at least version 0.7.0. You can install it with
```bash
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

## Miscelanous
And some official test, to make sure the solver returns the same results as pandapower
are performed in "pyklu/tests"
```bash
cd pyklu/tests
python -m unittest discover
```

This tests ensure that the results given by this simulator are consistent with the one given by pandapower when
using the Newton-Raphson algorithm, with a single slack bus, without enforcing q limits on the generators etc.
