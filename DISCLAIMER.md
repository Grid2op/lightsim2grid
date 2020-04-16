Copyright (c) 2020, RTE (https://www.rte-france.com)

See [AUTHORS.txt](AUTHORS.txt)

This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.

If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,

you can obtain one at http://mozilla.org/MPL/2.0/.

SPDX-License-Identifier: MPL-2.0

This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

Disclaimer
----------
This disclaimer only serves as a complement to the [`LICENSE`](LICENSE.md) file provided with it. It can not serve as a 
replacement of this file.

The simulator implemented in this package is made for speed mainly to serve as a grid2op backend (see
https://github.com/rte-france/grid2op for more information). We recall here that grid2op is a research testbed platform
aiming at emulating sequential decisions making in powergrids targeting mainly the "reinforcement learning"
community (though open to anyone).

This simulator is free and uses Eigen and KLU for performances optimization. It is not meant to be used as an 
independant for powersystem focus analysis. Indeed, the code provided in this package was made mainly for
developping AI focused controlers. Some of its limitations include, but are not limited to:

- it does not support multiple slack buses
- it can only use Newton-Raphson algorithm
- it does not enforce reactive power limits on generators
- it does not model AC/DC converters
- it does not model phase shifters
- transformers have fixed tap ratio (though it can be changed at initialization of the solver)
- shunts have fixed tap during the Newton-Raphson algorithm (though it can be changed at the initialization of the solver)
- only powerflow ("steady state") can be performed

Open source options for powerflow analysis
--------------------------------------------------
To get free of these limitations and be able to perform state of the art powerflow analysis, 
while still using open source softwares, we kindly recommend you to have a look at:
- [PowerModels](https://lanl-ansi.github.io/PowerModels.jl/stable/) described by its authors as "Julia/JuMP package 
  for Steady-State Power Network Optimization".
- [GridPack](https://www.gridpack.org/wiki/index.php/Main_Page) which is described as "An open source toolkit for 
  developing power grid simulation applications for high performance computing architectures"



