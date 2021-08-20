import grid2op
from lightsim2grid import LightSimBackend

env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend())
