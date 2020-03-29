from grid2op import make
from pyklu.PyKLUBackend import PyKLUBackend
backend = PyKLUBackend()
env = make(backend=backend)