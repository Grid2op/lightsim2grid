import grid2op
from lightsim2grid.LightSimBackend import LightSimBackend

env = grid2op.make(backend=LightSimBackend())
# obs = env.reset()
