__all__ = ["SparseLUSolver", "GaussSeidelSolver", "GaussSeidelSynchSolver", "DCSolver"]

try:
    from lightsim2grid_cpp import KLUSolver
    __all__.append("KLUSolver")
except Exception as exc_:
    # KLU is not available
    pass

from lightsim2grid_cpp import SparseLUSolver
from lightsim2grid_cpp import GaussSeidelSolver
from lightsim2grid_cpp import GaussSeidelSynchSolver
from lightsim2grid_cpp import DCSolver
