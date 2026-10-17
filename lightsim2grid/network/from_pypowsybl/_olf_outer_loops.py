# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The outer loops OpenLoadFlow would run with some pypowsybl loadflow parameters, as
lightsim2grid outer loops, with the parameters of its outer-loop driver.

OpenLoadFlow builds its list in ``AbstractAcOuterLoopConfig`` / ``DefaultAcOuterLoopConfig``
(or ``ExplicitAcOuterLoopConfig`` when ``outerLoopNames`` is given): which loops exist follows
its switches, each loop takes its parameters from OpenLoadFlow's. The same rules are applied
here, on the parameters given -- by default those of the installed pypowsybl, whose builds
ship different defaults.

The order of the default list has been the same since OpenLoadFlow 1.4: a version only
inserts the loops it adds at their place. :data:`_OLF_DEFAULT_ORDER` is that order (2.3.0's).
"""

import warnings
from typing import List, Optional, Tuple

import pypowsybl.loadflow as _lf

from lightsim2grid.lightsim2grid_cpp import (
    BaseOuterLoop,
    DistributedSlack,
    HvdcAcEmulationLimits,
    OuterLoopDriverParams,
    PhaseControl,
    ReactiveLimits,
    ShuntVoltageControl,
    TransformerVoltageControl,
    VoltageMonitoring,
)

#: OpenLoadFlow's default order (DefaultAcOuterLoopConfig), by its loop names
_OLF_DEFAULT_ORDER = (
    "DistributedSlack",
    "FreezingHvdcACEmulation",
    "AcHvdcAcEmulationLimits",
    "AreaInterchangeControl",
    "SecondaryVoltageControl",
    "VoltageMonitoring",
    "ReactiveLimits",
    "PhaseControl",
    "TransformerVoltageControl",
    "TransformerReactivePowerControl",
    "ShuntVoltageControl",
    "AutomationSystem",
)

#: names ``outerLoopNames`` accepts for the same loop (OpenLoadFlow renamed its HVDC one)
_OLF_NAME_ALIASES = {"HvdcAcEmulationLimits": "AcHvdcAcEmulationLimits"}

# PerUnit.SB: the base OpenLoadFlow turns a mismatch in MW / MVar into pu with
_OLF_SB_MVA = 100.

_PROVIDER_DEFAULTS = None


def _provider_defaults() -> dict:
    """``{name: default}`` (as strings) of the installed OpenLoadFlow's provider parameters"""
    global _PROVIDER_DEFAULTS
    if _PROVIDER_DEFAULTS is None:
        res = {}
        fn = getattr(_lf, "get_provider_parameters", None)
        if fn is not None:
            try:
                df = fn("OpenLoadFlow")
                if "default" in df.columns:
                    res = {name: str(df.loc[name, "default"]) for name in df.index}
            except Exception:
                res = {}
        _PROVIDER_DEFAULTS = res
    return _PROVIDER_DEFAULTS


class _Reader:
    """reads one parameter, the given value first, then the installed OpenLoadFlow's default,
    then lightsim2grid's (OpenLoadFlow's own) when the installed pypowsybl does not know it"""

    def __init__(self, parameters):
        self.parameters = parameters
        self.provider = dict(parameters.provider_parameters or {})
        self.defaults = _provider_defaults()

    def _raw(self, name):
        if name in self.provider:
            return str(self.provider[name])
        return self.defaults.get(name)

    def float(self, name, fallback):
        raw = self._raw(name)
        return fallback if raw is None else float(raw)

    def int(self, name, fallback):
        raw = self._raw(name)
        return fallback if raw is None else int(float(raw))

    def bool(self, name, fallback):
        raw = self._raw(name)
        return fallback if raw is None else raw.strip().lower() == "true"

    def str(self, name, fallback):
        raw = self._raw(name)
        return fallback if raw is None else raw.strip()

    def names(self, name):
        """the list given for `name`, None when it is not given (an empty value is an empty
        list: OpenLoadFlow then runs no loop at all)"""
        if name not in self.provider:
            return None
        raw = str(self.provider[name]).strip().strip("[]")
        return [el.strip() for el in raw.split(",") if el.strip()]

    def top(self, name, fallback):
        try:
            value = getattr(self.parameters, name)
        except AttributeError:
            return fallback
        return fallback if value is None else value


def _use_reactive_limits(read: _Reader) -> bool:
    value = read.top("use_reactive_limits", None)
    if value is not None:
        return bool(value)
    # its name before pypowsybl 1.x
    return not bool(read.top("no_generator_reactive_limits", False))


def _creators(read: _Reader, explicit: bool):
    """``{OpenLoadFlow loop name: a function returning the lightsim2grid loop, None when the
    parameters do not create it, or raising _Missing when lightsim2grid has no such loop}``.

    In the default list a control's mode is its ``*ControlMode`` parameter; in an explicit
    one (``outerLoopNames``) the name says the mode, as in ExplicitAcOuterLoopConfig."""

    def distributed_slack():
        if not bool(read.top("distributed_slack", True)):
            return None
        balance = read.top("balance_type", None)
        if balance is not None and getattr(balance, "name", str(balance)) != "PROPORTIONAL_TO_GENERATION_P_MAX":
            warnings.warn(f"olf_outer_loops: OpenLoadFlow's balance type {getattr(balance, 'name', balance)} "
                          "is not lightsim2grid's: its DistributedSlack loop shares the slack in proportion "
                          "to the units' maximum active power (PROPORTIONAL_TO_GENERATION_P_MAX).")
        behavior = read.str("slackDistributionFailureBehavior", "FAIL")
        if behavior not in ("FAIL", "LEAVE_ON_SLACK_BUS"):
            warnings.warn(f"olf_outer_loops: slackDistributionFailureBehavior {behavior} is not "
                          "lightsim2grid's: a slack its DistributedSlack loop cannot share fails the solve (FAIL).")
        default = DistributedSlack()
        return DistributedSlack(
            slack_bus_p_max_mismatch_mw=read.float("slackBusPMaxMismatch", default.slack_bus_p_max_mismatch_mw),
            fail_on_residue=behavior != "LEAVE_ON_SLACK_BUS")

    def hvdc_limits():
        return HvdcAcEmulationLimits() if bool(read.top("hvdc_ac_emulation", True)) else None

    def voltage_monitoring():
        return VoltageMonitoring() if read.bool("svcVoltageMonitoring", True) else None

    def reactive_limits():
        if not _use_reactive_limits(read):
            return None
        default = ReactiveLimits()
        criteria = read.str("newtonRaphsonStoppingCriteriaType", "UNIFORM_CRITERIA")
        if criteria == "PER_EQUATION_TYPE_CRITERIA":
            mismatch = read.float("maxReactivePowerMismatch", 1e-2) / _OLF_SB_MVA
        else:
            mismatch = read.float("newtonRaphsonConvEpsPerEq", default.max_reactive_power_mismatch)
        return ReactiveLimits(
            max_pq_pv_switch=read.int("reactiveLimitsMaxPqPvSwitch", default.max_pq_pv_switch),
            max_reactive_power_mismatch=mismatch,
            robust_mode=read.bool("voltageRemoteControlRobustMode", default.robust_mode),
            min_realistic_voltage=read.float("minRealisticVoltage", default.min_realistic_voltage),
            max_realistic_voltage=read.float("maxRealisticVoltage", default.max_realistic_voltage))

    def phase_control(mode=None):
        if not bool(read.top("phase_shifter_regulation_on", False)):
            return None
        mode = mode or read.str("phaseShifterControlMode", "CONTINUOUS_WITH_DISCRETISATION")
        if mode != "CONTINUOUS_WITH_DISCRETISATION":
            raise _Missing(f"PhaseControl in mode {mode}")
        return PhaseControl()

    def transformer_voltage_control(mode=None):
        if not bool(read.top("transformer_voltage_control_on", False)):
            return None
        mode = mode or read.str("transformerVoltageControlMode", "INCREMENTAL_VOLTAGE_CONTROL")
        if mode != "AFTER_GENERATOR_VOLTAGE_CONTROL":
            raise _Missing(f"TransformerVoltageControl in mode {mode}")
        default = TransformerVoltageControl()
        return TransformerVoltageControl(
            use_initial_tap_position=read.bool("transformerVoltageControlUseInitialTapPosition",
                                               default.use_initial_tap_position),
            max_controlled_nominal_voltage=read.float("generatorVoltageControlMinNominalVoltage",
                                                      default.max_controlled_nominal_voltage))

    def shunt_voltage_control(mode=None):
        if not bool(read.top("shunt_compensator_voltage_control_on", False)):
            return None
        mode = mode or read.str("shuntVoltageControlMode", "WITH_GENERATOR_VOLTAGE_CONTROL")
        if mode != "WITH_GENERATOR_VOLTAGE_CONTROL":
            raise _Missing(f"ShuntVoltageControl in mode {mode}")
        return ShuntVoltageControl()

    def missing_if(flag, name):
        def create():
            if flag():
                raise _Missing(name)
            return None
        return create

    return {
        "DistributedSlack": distributed_slack,
        "FreezingHvdcACEmulation": missing_if(lambda: read.bool("startWithFrozenACEmulation", False),
                                              "FreezingHvdcACEmulation"),
        "AcHvdcAcEmulationLimits": hvdc_limits,
        "AreaInterchangeControl": missing_if(lambda: read.bool("areaInterchangeControl", False),
                                             "AreaInterchangeControl"),
        "SecondaryVoltageControl": missing_if(lambda: read.bool("secondaryVoltageControl", False),
                                              "SecondaryVoltageControl"),
        "VoltageMonitoring": voltage_monitoring,
        "ReactiveLimits": reactive_limits,
        "PhaseControl": lambda: phase_control("CONTINUOUS_WITH_DISCRETISATION" if explicit else None),
        "IncrementalPhaseControl": lambda: phase_control("INCREMENTAL"),
        "TransformerVoltageControl": lambda: transformer_voltage_control(
            "AFTER_GENERATOR_VOLTAGE_CONTROL" if explicit else None),
        "SimpleTransformerVoltageControl": lambda: transformer_voltage_control("WITH_GENERATOR_VOLTAGE_CONTROL"),
        "IncrementalTransformerVoltageControl": lambda: transformer_voltage_control("INCREMENTAL_VOLTAGE_CONTROL"),
        "TransformerReactivePowerControl": missing_if(lambda: read.bool("transformerReactivePowerControl", False),
                                                      "TransformerReactivePowerControl"),
        "IncrementalTransformerReactivePowerControl": missing_if(
            lambda: read.bool("transformerReactivePowerControl", False), "IncrementalTransformerReactivePowerControl"),
        "ShuntVoltageControl": lambda: shunt_voltage_control("WITH_GENERATOR_VOLTAGE_CONTROL" if explicit else None),
        "IncrementalShuntVoltageControl": lambda: shunt_voltage_control("INCREMENTAL_VOLTAGE_CONTROL"),
        "AutomationSystem": missing_if(lambda: read.bool("simulateAutomationSystems", False), "AutomationSystem"),
    }


class _Missing(Exception):
    """an outer loop OpenLoadFlow would run that lightsim2grid does not have"""


def olf_outer_loops(parameters: Optional[_lf.Parameters] = None) -> Tuple[List[BaseOuterLoop], OuterLoopDriverParams]:
    """The outer loops OpenLoadFlow would run with ``parameters``, in its order, and the
    parameters of its outer-loop driver.

    Which loops are created, and with what parameters, follows OpenLoadFlow's own rules:
    ``distributed_slack``, ``hvdc_ac_emulation``, ``svcVoltageMonitoring``,
    ``use_reactive_limits``, ``phase_shifter_regulation_on``, ``transformer_voltage_control_on``
    and ``shunt_compensator_voltage_control_on`` create theirs; ``outerLoopNames``, when given,
    gives the list and its order instead of the default one (empty: no loop at all). A loop these parameters create
    that lightsim2grid does not have (an incremental mode, the area interchange control, the
    automation systems...) is left out, with a warning.

    Parameters
    ----------
    parameters: ``pypowsybl.loadflow.Parameters``, optional
        The OpenLoadFlow parameters. ``None`` (default): the defaults of the installed
        pypowsybl, which differ between its builds. A provider parameter left unset takes the
        installed OpenLoadFlow's default too.

    Returns
    -------
    loops: ``list``
        The lightsim2grid outer loops, in OpenLoadFlow's order (see
        :func:`lightsim2grid.network.LSGrid.add_outer_loop`)
    driver: :class:`lightsim2grid.algorithm.OuterLoopDriverParams`
        The parameters of the outer-loop driver (see
        :func:`lightsim2grid.network.LSGrid.set_outer_loop_driver_params`)
    """
    if parameters is None:
        parameters = _lf.Parameters()
    read = _Reader(parameters)
    explicit = read.names("outerLoopNames")
    creators = _creators(read, explicit is not None)
    names = list(_OLF_DEFAULT_ORDER) if explicit is None else [_OLF_NAME_ALIASES.get(name, name) for name in explicit]
    unknown = [name for name in names if name not in creators]
    if unknown:
        raise ValueError(f"olf_outer_loops: unknown OpenLoadFlow outer loop(s) {unknown} in outerLoopNames")
    if len(set(names)) != len(names):
        raise ValueError(f"olf_outer_loops: outerLoopNames holds the same loop twice: {explicit}")

    loops, missing = [], []
    for name in names:
        try:
            loop = creators[name]()
        except _Missing as exc:
            missing.append(str(exc))
            continue
        if loop is not None:
            loops.append(loop)
    if missing:
        warnings.warn(f"olf_outer_loops: OpenLoadFlow would also run {missing}, which lightsim2grid does "
                      "not have: they are left out.")

    default = OuterLoopDriverParams()
    driver = OuterLoopDriverParams(
        max_outer_iterations=read.int("maxOuterLoopIterations", default.max_outer_iterations),
        voltage_remote_control_robust_mode=read.bool("voltageRemoteControlRobustMode",
                                                     default.voltage_remote_control_robust_mode),
        min_realistic_voltage=read.float("minRealisticVoltage", default.min_realistic_voltage),
        max_realistic_voltage=read.float("maxRealisticVoltage", default.max_realistic_voltage),
        min_nominal_voltage_realistic_check=read.float("minNominalVoltageRealisticVoltageCheck",
                                                       default.min_nominal_voltage_realistic_check))
    return loops, driver


def set_olf_outer_loops(grid, parameters: Optional[_lf.Parameters] = None) -> None:
    """Make ``grid`` run the outer loops OpenLoadFlow would run with ``parameters``, with the
    same driver parameters (see :func:`olf_outer_loops`): its list of outer loops is replaced.
    ``None`` (default): the defaults of the installed pypowsybl."""
    loops, driver = olf_outer_loops(parameters)
    grid.clear_outer_loops()
    for loop in loops:
        grid.add_outer_loop(loop)
    grid.set_outer_loop_driver_params(driver)
