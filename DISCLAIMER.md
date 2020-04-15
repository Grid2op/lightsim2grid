Copyright (c) 2020, RTE (https://www.rte-france.com)
See AUTHORS.txt
This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
you can obtain one at http://mozilla.org/MPL/2.0/.
SPDX-License-Identifier: MPL-2.0
This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

This disclaimer only serves as a complement to the `LICENSE.md` file provided with it. It can not serve as a 
replacement of this file.


The simulator implemented in this package is made for speed mainly to serve as a grid2op backend (see
https://github.com/rte-france/grid2op for more information). We recall here that grid2op is a research testbed platform
aiming at emulating sequential decisions making in powergrids.

This powerflow is free and uses Eigen and KLU for performances optimization. However, it is a really unrealistic
powerflow with many limitations. Some of these limitations include, but are not limited to:

- it does not support multiple slack bus
- it can only use Newton-Raphson algorithm
- it does not enforce reactive power limits on generators
- it does not model AC/DC converters
- it does not model phase shifters
- transformers have fixed tap ratio (though it can be changed at initialization of the solver)
- shunts have fixed tap during the Newton-Raphson algorithm (though it can be changed at the initialization of the solver)
- only powerflow ("steady state") can be performed




