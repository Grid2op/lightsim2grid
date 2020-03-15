import pandapower as pp
import pandapower.networks as pn

grid = pn.case30()
pp.runpp(grid)