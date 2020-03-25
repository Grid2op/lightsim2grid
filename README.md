# powerflow_klu
Provide a faster powerflow for pandapower using c++ klu (I hope)

## installation
You need to:
- clone this repository
- get the code of Eigen and SparseSuite (mandatory for compilation)
- compile a piece of SparseSuite
- install the package

All of that can be done, on a computer with make (so most likely MacOS or GNU/Linux based):
```bash
git clone https://github.com/BDonnot/powerflow_klu.git
cd powerflow_klu
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

### 1. replacement of pandapower "newtonpf" method
Suppose you somehow get:
- `Ybus` the admittance matrix of your powersystem
- `V0` the (complex) voltage vector at each bus
- `Sbus` the (complex) power absorb at each bus
- `ppci` a ppc internal pandapower test case
- `pv` list of PV buses
- `pq` list of PQ buses
- `options` list of pandapower powerflow "options"

You can define replace the `newtonpf` function of `pandapower.pandapower.newtonpf` function with the following
piece of code:
```python
from pyklu.newtonpf import newtonpf
V, converged, iterations, J = newtonpf(Ybus, V, Sbus, pv, pq, ppci, options)
```

This function uses the KLU algorithm and a c++ implementation of a Newton solver for speed.

### 2. attempt to optimize the code
You can also use a different class, that should be able to avoid
un necessary conversion from pandapower to the format specified above.

WORK IN PROGRESS USE AT YOUR OWN RISK

```python
from pyklu.compute_powerflow import KLU4Pandapower
import pandapower.networks as pn
grid2 = pn.case9241pegase()
cpp_solver = KLU4Pandapower()
nb_max_newton_it = 10  #maximum number of iteration for newton raphson
# first call, you need to "reset the solver"
cpp_solver.runpp(grid2, 
                 max_iteration=nb_max_newton_it,
                 need_reset=True   # reset the KLU solver to an original state, need to be done each time the Ymatrix is changed (might be slow)
                 )

# second and further calls, as long as the admittance matrix does not change
cpp_solver.runpp(grid2, 
                 max_iteration=nb_max_newton_it,
                 need_reset=False  # reuse previous voltage angles, jacobian matrix, and admittance matrix
                 )
```

## Miscelanous
WORK IN PROGRESS USE AT YOUR OWN RISK

You can run some tests with:
```bash
python old/test_klu.py
```

And some official test, to make sure the solver returns the same results as pandapower
are performed in "pyklu/tests"
```bash
cd pyklu/tests
python -m unittest discover
```


