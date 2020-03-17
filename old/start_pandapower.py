import pandapower as pp
import pandapower.networks as pn

grid = pn.case30()
grid = pn.case118()
pp.runpp(grid, numba=False)