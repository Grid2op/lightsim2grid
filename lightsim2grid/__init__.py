__version__ = "0.4.0"

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
