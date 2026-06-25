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

The factory is written to be consistent across **every** pypowsybl version from
1.0.0 to 1.15.0. The pypowsybl ``Parameters`` constructor and the OLF provider
parameters changed over time (renamed kwargs, new enums, parameters that did not
exist yet, and values rejected by older releases). Rather than hard-code one
version, the factory introspects the installed pypowsybl and only passes what it
accepts:

* top-level kwargs are filtered against ``Parameters.__init__`` and renamed
  parameters fall back to their old spelling (``use_reactive_limits`` ->
  ``no_generator_reactive_limits``, ``shunt_compensator_voltage_control_on`` ->
  ``simul_shunt``, ``component_mode`` -> ``connected_component_mode``);
* enums are resolved from wherever the installed version exposes them;
* provider parameters are filtered against the OLF metadata (valid name *and*,
  for enum-valued ones, valid value);
* the empty ``outerLoopNames`` allow-list is only used from pypowsybl 1.4.0 on,
  because older OLF rejects the empty value ("Unknown outer loop ''"). Below
  1.4.0 the explicit trigger flags already disable every outer loop that exists.
"""

import inspect
from typing import Iterable, Optional, Union

import pypowsybl as _pp
import pypowsybl.loadflow as _lf


def _ver_tuple(v: str):
    out = []
    for part in str(v).split("."):
        num = ""
        for ch in part:
            if ch.isdigit():
                num += ch
            else:
                break
        out.append(int(num) if num else 0)
    return tuple(out)


_PARAM_SIG = set(inspect.signature(_lf.Parameters.__init__).parameters)
_PYPOWSYBL_VER = _ver_tuple(_pp.__version__)
# OLF rejects an empty ``outerLoopNames`` value before 1.4.0 (and the parameter
# does not exist at all in 1.0.0); only use the empty-allow-list trick where it
# is accepted.
_OUTERLOOP_EMPTY_OK = _PYPOWSYBL_VER >= (1, 4)

# Lazily-computed cache of the OLF provider parameter metadata
# ({name: set(possible_values) or None}); ``False`` means "not yet computed".
_OLF_PROVIDER_PARAMS = False


def _enum(name: str):
    """Resolve an enum class by name from wherever pypowsybl exposes it."""
    for mod in (_lf, getattr(_pp, "_pypowsybl", None)):
        if mod is not None and hasattr(mod, name):
            return getattr(mod, name)
    return None


def _olf_provider_params():
    """Return ``{name: set(possible_values) or None}`` for OpenLoadFlow, or
    ``None`` if the installed pypowsybl does not expose the metadata."""
    global _OLF_PROVIDER_PARAMS
    if _OLF_PROVIDER_PARAMS is not False:
        return _OLF_PROVIDER_PARAMS

    result = None
    fn = getattr(_lf, "get_provider_parameters", None)
    if fn is not None:
        try:
            df = fn("OpenLoadFlow")
            has_pv = "possible_values" in df.columns
            result = {}
            for name in df.index:
                vals = None
                if has_pv:
                    raw = df.loc[name, "possible_values"]
                    if isinstance(raw, str) and raw.strip() not in ("", "[]"):
                        vals = {x.strip() for x in raw.strip("[]").split(",") if x.strip()}
                result[name] = vals
        except Exception:
            result = None
    if result is None:
        fn = getattr(_lf, "get_provider_parameters_names", None)
        if fn is not None:
            try:
                result = {n: None for n in fn("OpenLoadFlow")}
            except Exception:
                result = None

    _OLF_PROVIDER_PARAMS = result
    return result


# The full "ZERO outer loops" OLF provider-parameter wish-list. Entries that the
# installed OLF does not know (or whose value it rejects) are filtered out.
_DESIRED_PROVIDER = {
    # ---- empty allow-list: ZERO outer loops, future-proof (>= 1.4.0) ----
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


def get_pypowsybl_loopfree_parameters(
    slack_bus_ids: Optional[Union[str, Iterable[str]]] = None,
    **overrides,
) -> _lf.Parameters:
    """Build a fresh :class:`pypowsybl.loadflow.Parameters` with every OLF outer
    loop disabled.

    The returned parameters disable distributed slack, reactive limits,
    transformer / shunt / phase-shifter voltage control, area-interchange and
    secondary-voltage control, automation systems, etc. From pypowsybl 1.4.0 on,
    the empty ``outerLoopNames`` allow-list also registers *zero* outer loops
    regardless of which loops a future OLF release adds; on older pypowsybl the
    explicit trigger flags below disable every outer loop that exists.

    The factory is consistent across pypowsybl 1.0.0 .. 1.15.0: it only passes
    constructor keywords and provider parameters the installed version accepts
    (see the module docstring).

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
        on the returned object (e.g. ``distributed_slack=True``). Only applied if
        the installed pypowsybl accepts that keyword.

    Returns
    -------
    pypowsybl.loadflow.Parameters
        A new object on each call (no shared mutable state).

    Notes
    -----
    ``maxOuterLoopIterations`` is set to ``"1"`` because OLF rejects ``"0"``;
    with an empty ``outerLoopNames`` the iteration cap is moot. Verified to
    converge on IEEE-14 against every pypowsybl from 1.0.0 to 1.15.0.
    """
    kwargs = {}

    def put(name, value):
        if name in _PARAM_SIG:
            kwargs[name] = value
            return True
        return False

    put("distributed_slack", False)
    put("transformer_voltage_control_on", False)
    put("phase_shifter_regulation_on", False)
    put("write_slack_bus", True)
    put("dc_use_transformer_ratio", True)
    put("twt_split_shunt_admittance", True)

    vim = _enum("VoltageInitMode")
    if vim is not None:
        put("voltage_init_mode", vim.UNIFORM_VALUES)
    bt = _enum("BalanceType")
    if bt is not None:
        # forced; inert with distributed slack off
        put("balance_type", bt.PROPORTIONAL_TO_GENERATION_P_MAX)

    # reactive limits: new ``use_reactive_limits`` vs old ``no_generator_reactive_limits``
    if not put("use_reactive_limits", False):
        put("no_generator_reactive_limits", True)
    # shunt voltage control: new ``shunt_compensator_voltage_control_on`` vs old ``simul_shunt``
    if not put("shunt_compensator_voltage_control_on", False):
        put("simul_shunt", False)

    read_slack_bus = slack_bus_ids is None
    put("read_slack_bus", read_slack_bus)

    # connected component: ``component_mode`` (>= 1.15) vs ``connected_component_mode``
    if "component_mode" in _PARAM_SIG:
        cm = _enum("ComponentMode")
        if cm is not None and hasattr(cm, "MAIN_CONNECTED"):
            kwargs["component_mode"] = cm.MAIN_CONNECTED
    elif "connected_component_mode" in _PARAM_SIG:
        ccm = _enum("ConnectedComponentMode")
        if ccm is not None and hasattr(ccm, "MAIN"):
            kwargs["connected_component_mode"] = ccm.MAIN

    desired = dict(_DESIRED_PROVIDER)
    if not _OUTERLOOP_EMPTY_OK:
        desired.pop("outerLoopNames", None)
    if slack_bus_ids is not None:
        if not isinstance(slack_bus_ids, str):
            slack_bus_ids = ",".join(str(s) for s in slack_bus_ids)
        desired["slackBusSelectionMode"] = "NAME"
        desired["slackBusesIds"] = f"{slack_bus_ids}"

    valid = _olf_provider_params()
    if valid is not None:
        provider_parameters = {}
        for name, value in desired.items():
            if name not in valid:
                continue  # parameter unknown in this pypowsybl version
            possible = valid[name]
            if possible is not None and value not in possible:
                continue  # value not accepted in this version
            provider_parameters[name] = value
    else:
        provider_parameters = desired
    put("provider_parameters", provider_parameters)

    # only forward overrides the installed Parameters constructor accepts
    for name, value in overrides.items():
        put(name, value)
    return _lf.Parameters(**kwargs)
