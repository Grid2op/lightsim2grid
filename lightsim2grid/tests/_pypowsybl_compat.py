# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Build the same IIDM data on both sides of pypowsybl 1.17.

pypowsybl 1.17 moved two pieces of voltage-control data out of extensions and into
columns of the element tables: the ``voltagePerReactivePowerControl`` extension of a static
var compensator became its ``slope`` column, and the ``voltageRegulation`` extension of a
battery became its ``regulating`` / ``regulation_mode`` / ``regulation_value`` columns.
The extensions no longer exist there. The test fixtures go through these helpers so that
they describe one network whichever API the installed pypowsybl offers; which one it is
is detected from the extensions pypowsybl lists, not from its version number."""

import pypowsybl as pp


def _has_extension(name):
    return name in pp.network.get_extensions_names()


def set_svc_slope(net, svc_id, slope):
    """The voltage / reactive-power slope of a static var compensator, in kV/MVar."""
    if _has_extension("voltagePerReactivePowerControl"):
        net.create_extensions("voltagePerReactivePowerControl", id=svc_id, slope=slope)
    else:
        net.update_static_var_compensators(id=svc_id, slope=slope)


def set_battery_voltage_regulation(net, battery_id, regulating, target_v):
    """A battery regulating the voltage of its own bus at ``target_v`` (kV)."""
    if _has_extension("voltageRegulation"):
        net.create_extensions("voltageRegulation", id=battery_id, voltage_regulator_on=bool(regulating),
                              target_v=target_v)
    else:
        net.update_batteries(id=battery_id, regulation_mode="VOLTAGE", regulation_value=target_v,
                             regulating=bool(regulating))
