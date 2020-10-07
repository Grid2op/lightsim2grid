__version__ = "0.3.0"

__all__ = ["newtonpf", "SolverType"]

from lightsim2grid.initGridModel import SolverType

try:
    from lightsim2grid.LightSimBackend import LightSimBackend
    __all__.append("LightSimBackend")
except ImportError:
    # grid2op is not installed, the Backend will not be available
    pass
