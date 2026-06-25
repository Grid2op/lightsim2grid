# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Canonical PowSyBl Open Load Flow (OLF) parameters that disable *every* outer
loop, so that OLF solves the same single-shot power-flow problem lightsim2grid
solves.

This is the single source of truth used by the lightsim2grid benchmarks, the
documentation, the test suite and the OLF-vs-lightsim2grid comparison helper.
Disabling the outer loops is what makes a *raw* (un-baked) network give
consistent results between the two engines on grids where no outer loop would
trigger anyway; on grids where the loops do act, combine it with
:func:`lightsim2grid.network.bake_outer_loops`.
"""

import inspect
from typing import Iterable, Optional, Union

import pypowsybl.loadflow as _lf

# pypowsybl >= 1.15 renamed ``connected_component_mode`` -> ``component_mode``
# (with the ``ComponentMode`` enum). Pick whichever the installed version
# exposes so we neither trigger the deprecation warning on recent pypowsybl nor
# break on older releases. Both spellings select the *main connected* component.
_HAS_COMPONENT_MODE = "component_mode" in inspect.signature(_lf.Parameters.__init__).parameters


def get_pypowsybl_loopfree_parameters(
    slack_bus_ids: Optional[Union[str, Iterable[str]]] = None,
    **overrides,
) -> _lf.Parameters:
    """Build a fresh :class:`pypowsybl.loadflow.Parameters` with every OLF outer
    loop disabled.

    The returned parameters disable distributed slack, reactive limits,
    transformer / shunt / phase-shifter voltage control, area-interchange and
    secondary-voltage control, automation systems, etc. The mechanism that
    makes this future-proof is the empty ``outerLoopNames`` allow-list: it
    registers *zero* outer loops regardless of which loops a future OLF release
    adds. Every individual trigger flag is also forced off (belt and
    suspenders), and the control modes are pushed to the inner Newton-Raphson
    loop so nothing escalates to an outer loop.

    Parameters
    ----------
    slack_bus_ids : str or iterable of str, optional
        If given, the slack is pinned by NAME on these bus(es)
        (``slackBusSelectionMode = "NAME"``) and ``read_slack_bus`` is turned
        off. This is what the benchmarks need to use the same slack as
        lightsim2grid. If ``None`` (default), the slack is read from the network
        (``read_slack_bus = True``).
    **overrides
        Any top-level :class:`pypowsybl.loadflow.Parameters` keyword to override
        on the returned object (e.g. ``distributed_slack=True``).

    Returns
    -------
    pypowsybl.loadflow.Parameters
        A new object on each call (no shared mutable state).

    Notes
    -----
    Some ``provider_parameters`` keys are OLF-specific; other load-flow
    providers may ignore or warn on them. ``maxOuterLoopIterations`` is set to
    ``"1"`` because OLF rejects ``"0"``; with an empty ``outerLoopNames`` the
    iteration cap is moot (there is no outer loop to iterate). Verified against
    pypowsybl 1.15.0.
    """
    provider_parameters = {
        # ---- empty allow-list: ZERO outer loops, future-proof ----
        "outerLoopNames": "",
        # ---- every outer-loop creation trigger forced off (belt and suspenders) ----
        "areaInterchangeControl": "false",
        "svcVoltageMonitoring": "false",
        "secondaryVoltageControl": "false",
        "transformerReactivePowerControl": "false",
        "simulateAutomationSystems": "false",
        "voltageRemoteControl": "false",
        "generatorVoltageRemoteControl": "false",
        # ---- control modes pushed to inner-loop / harmless ----
        "transformerVoltageControlMode": "WITH_GENERATOR_VOLTAGE_CONTROL",
        "shuntVoltageControlMode": "WITH_GENERATOR_VOLTAGE_CONTROL",
        "phaseShifterControlMode": "CONTINUOUS_WITH_DISCRETISATION",
        # ---- NR inner-solver settings forced (the only thing left active) ----
        "maxNewtonRaphsonIterations": "15",
        "newtonRaphsonConvEpsPerEq": "1.0E-4",
        "stateVectorScalingMode": "NONE",
        "lineSearchStateVectorScalingMaxIteration": "10",  # inert with NONE; forced anyway
        # ---- hard backstop (0 is rejected by OLF; 1 is moot with empty allow-list) ----
        "maxOuterLoopIterations": "1",
        "slackDistributionFailureBehavior": "FAIL",
    }

    read_slack_bus = slack_bus_ids is None
    if slack_bus_ids is not None:
        if not isinstance(slack_bus_ids, str):
            slack_bus_ids = ",".join(str(s) for s in slack_bus_ids)
        provider_parameters["slackBusSelectionMode"] = "NAME"
        provider_parameters["slackBusesIds"] = f"{slack_bus_ids}"

    kwargs = dict(
        distributed_slack=False,
        balance_type=_lf.BalanceType.PROPORTIONAL_TO_GENERATION_P_MAX,  # forced; inert with slack off
        use_reactive_limits=False,
        phase_shifter_regulation_on=False,
        transformer_voltage_control_on=False,
        shunt_compensator_voltage_control_on=False,
        voltage_init_mode=_lf.VoltageInitMode.UNIFORM_VALUES,
        read_slack_bus=read_slack_bus,
        write_slack_bus=True,
        dc_use_transformer_ratio=True,
        twt_split_shunt_admittance=True,
        provider_parameters=provider_parameters,
    )
    if _HAS_COMPONENT_MODE:
        kwargs["component_mode"] = _lf.ComponentMode.MAIN_CONNECTED
    else:
        kwargs["connected_component_mode"] = _lf.ConnectedComponentMode.MAIN
    kwargs.update(overrides)
    return _lf.Parameters(**kwargs)
