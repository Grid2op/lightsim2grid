__version__ = "0.5.4"

__all__ = ["newtonpf", "SolverType", "solver"]

# import directly from c++ module
import lightsim2grid.solver
from lightsim2grid.solver import SolverType

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
    # grid2op is not installed, the PhysicalLawChecker will not be available
    pass

try:
    from lightsim2grid.timeSerie import TimeSerie
    __all__.append("TimeSerie")
except ImportError as exc_:
    # grid2op is not installed, the TimeSeries module will not be available
    pass

try:
    from lightsim2grid.securtiyAnalysis import SecurityAnalysis
    __all__.append("SecurityAnalysis")
except ImportError as exc_:
    # grid2op is not installed, the SecurityAnalysis module will not be available
    pass
