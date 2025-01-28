Disclaimer
==================

Disclaimer
----------

This disclaimer only serves as a complement to the
`LICENSE <https://github.com/Grid2Op/lightsim2grid/blob/master/LICENSE>`_ file provided with it. It can not serve as a
replacement of this file.

The simulator implemented in this package is made for speed mainly to serve as a grid2op backend (see
https://github.com/rte-france/grid2op for more information). We recall here that grid2op is a research testbed platform
aiming at emulating sequential decisions making in powergrids targeting mainly the "reinforcement learning"
community (though open to anyone).

This simulator is free and uses Eigen (optionally KLU) for performances optimization. It is not meant to be used as an
independant for powersystem focus analysis. Indeed, the code provided in this package was made mainly for
developping AI focused controlers. Some of its limitations include, but are not limited to:

- it does not enforce reactive power limits on generators
- it does not model AC/DC converters
- transformers have fixed tap ratio (though it can be changed at initialization of the solver)
- shunts have fixed tap during the Newton-Raphson algorithm (though it can be changed at the initialization of the solver)
- only powerflow ("steady state") can be performed

Open source options for powerflow analysis
--------------------------------------------------
To get free of these limitations and be able to perform state of the art powerflow analysis,
while still using open source software, we kindly recommend you to have a look at:

- `Matpower <https://matpower.org/>`_ which is a "free, open-source tools for electric power system simulation and
  optimization "
- `Pandapower <https://www.pandapower.org/>`_ that is "An easy to use open source tool for power system modeling,
  analysis and optimization with a high degree of automation. "
- `PowerModels <https://lanl-ansi.github.io/PowerModels.jl/stable/>`_ described by its authors as "Julia/JuMP package
  for Steady-State Power Network Optimization".
- `GridPack <https://www.pnnl.gov/projects/gridpacktm-open-source-framework-developing-high-performance-computing-simulations-power>`_ which is described as "An open source toolkit for
  developing power grid simulation applications for high performance computing architectures"
- `Dynaωo <https://github.com/dynawo/dynawo>`_ an hybrid C++/Modelica suite of simulation tools for
  power systems.

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`