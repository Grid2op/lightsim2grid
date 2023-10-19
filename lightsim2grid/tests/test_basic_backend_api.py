# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import unittest

try:
    from grid2op._create_test_suite import create_test_suite
    CAN_PERFORM_THESE = True
except ImportError as exc_:
    CAN_PERFORM_THESE = False

from lightsim2grid import LightSimBackend

if CAN_PERFORM_THESE:
    def this_make_backend(self, detailed_infos_for_cascading_failures=False):
        return LightSimBackend(
                detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures
            )
    add_name_cls = "test_LightSimBackend"

    res = create_test_suite(make_backend_fun=this_make_backend,
                            add_name_cls=add_name_cls,
                            add_to_module=__name__,
                            extended_test=False,  # for now keep `extended_test=False` until all problems are solved
                            )
else:
    print("Have you installed grid2op in dev / editable mode ? We cannot make the `create_test_suite` :-(")
    
# and run it with `python -m unittest gridcal_backend_tests.py`
if __name__ == "__main__":
    unittest.main()
