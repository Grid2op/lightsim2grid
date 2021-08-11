__version__ = "0.5.3"

__all__ = ["newtonpf", "SolverType", "solver"]

# import directly from c++ module
from lightsim2grid_cpp import SolverType
import lightsim2grid.solver

try:
    from lightsim2grid.LightSimBackend import LightSimBackend
    __all__.append("LightSimBackend")
except ImportError:
    # grid2op is not installed, the Backend will not be available
    pass

try:
    from lightsim2grid.physical_law_checker import PhysicalLawChecker
    __all__.append("PhysicalLawChecker")
except ImportError as exc_:
    # grid2op is not installed, the Backend will not be available
    pass
    print(exc_)
