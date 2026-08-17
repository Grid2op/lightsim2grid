# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Remove PowSyBl Open Load Flow (OLF) *outer loops* from a set of pypowsybl
loadflow parameters, without touching anything else.

:func:`remove_outer_loops` is the single source of truth for "make OLF solve
the same single-shot power-flow problem lightsim2grid solves". It takes a
caller-supplied :class:`pypowsybl.loadflow.Parameters` and strips out the
outer-loop *mechanism* only: it does not change reactive limits, remote
voltage control, voltage initialization, balance type, or any other modeling
choice, because none of those are outer loops -- forcing them changes the
solved answer for reasons unrelated to outer loops (and, for
``voltage_init_mode`` in particular, can push OLF into a different
Newton-Raphson root on a multi-root grid, contaminating a comparison meant to
isolate the outer-loop effect).

Background: OLF's outer loops fall into two families.

* Some have an inline alternative: ``transformerVoltageControlMode``,
  ``shuntVoltageControlMode`` and ``phaseShifterControlMode`` can each be
  switched from an outer-loop-based mode to a mode that folds the same
  control into the inner Newton-Raphson equations. Switching the mode keeps
  the modeled control (e.g. transformer voltage regulation) *active*, it just
  removes the separate outer-loop pass -- this is real outer-loop removal
  that changes nothing else, and it works on every pypowsybl version.
* Others (``DistributedSlack``, ``ReactiveLimits``, ``VoltageMonitoring``,
  ``SecondaryVoltageControl``, ``AreaInterchangeControl``,
  ``AutomationSystem``) have no inline equivalent. From pypowsybl 1.4.0 on,
  the ``outerLoopNames`` provider parameter is an *allow-list*: only loops
  named there are ever created, regardless of any other trigger flag
  (verified empirically: ``distributed_slack=True`` with an empty
  ``outerLoopNames`` does not distribute). So on 1.4.0+, this allow-list is
  the only lever this module needs to touch for that family for an AC solve.
  Below 1.4.0, where OLF rejects an empty ``outerLoopNames`` (and the
  parameter does not exist at all in 1.0.0), the individual creation-trigger
  flags are the only way to stop a loop, so they are forced off there instead.

  ``DistributedSlack`` and ``ReactiveLimits`` are the two exceptions to
  "the allow-list is the only lever": both also have a top-level
  ``Parameters`` flag (``distributed_slack``, ``use_reactive_limits``) that
  is not reliably gated by ``outerLoopNames``. ``run_dc`` always honors
  ``distributed_slack`` directly, with no outer-loop pass and no
  ``outerLoopNames`` gating at all -- it still distributes the slack across
  generators with ``distributed_slack`` left at its default (``True``) even
  with an empty allow-list. ``use_reactive_limits`` is not DC-specific but is
  just as unreliable: on pypowsybl 1.9.0, ``run_ac`` still clips a generator
  to its reactive limit with an empty ``outerLoopNames`` and
  ``use_reactive_limits`` left at its default (``True``), producing a
  generator ``Q`` a few 1e-3 MVAr off from lightsim2grid's (which never
  clips) -- confirmed fixed by forcing ``use_reactive_limits=False``, and
  confirmed absent on pypowsybl 1.15.0 for the same grid/scenario, so this
  looks like an OLF-side version difference in how much the allow-list
  actually covers, not a stable "AC is always allow-list-gated" rule. So
  these two top-level flags are always forced off here (unless kept), on
  every pypowsybl version, in addition to whatever the allow-list/trigger-flag
  branch below does for the other, allow-list-only loops.

The factory (both this module's ``remove_outer_loops`` and the
``get_pypowsybl_loopfree_*`` convenience wrappers below) is written to be
consistent across every pypowsybl version from 1.0.0 to 1.15.0: it only
passes constructor keywords and provider parameters the installed version
accepts, resolving renamed kwargs (``use_reactive_limits`` ->
``no_generator_reactive_limits``, ``shunt_compensator_voltage_control_on`` ->
``simul_shunt``) and enums from wherever the installed version exposes them.
"""

import inspect
from typing import Iterable, Optional, Union

import pypowsybl as _pp
import pypowsybl.loadflow as _lf
from packaging import version as _version


_PARAM_SIG = set(inspect.signature(_lf.Parameters.__init__).parameters)
_PYPOWSYBL_VER = _version.parse(_pp.__version__)
# OLF rejects an empty ``outerLoopNames`` value before 1.4.0 (and the parameter
# does not exist at all in 1.0.0); only use the allow-list where it is accepted.
_OUTERLOOP_EMPTY_OK = _PYPOWSYBL_VER >= _version.parse("1.4")

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


# Outer loops with an inline (inner-Newton-Raphson) alternative: switching the
# mode keeps the modeled control active while removing the separate outer-loop
# pass. Works on every pypowsybl version, so these are always applied unless
# the caller explicitly asks to ``keep`` the outer-loop-based mode. This *is*
# a deliberate override of whatever mode ``parameters`` already had -- not a
# bug: this whole function's job is to make OLF solve the same single-shot,
# outer-loop-free problem lightsim2grid solves (see module docstring), so an
# input that requests an outer-loop-based mode is exactly what gets changed
# by default. Yes, that can in principle move OLF to a different NR root than
# the caller's original (with-loops) parameters would have (see
# ``_olf_bake.py``'s "should agree" caveat); pass the loop's name in ``keep``
# if that is not acceptable for a given comparison.
_INLINE_MODE = {
    "TransformerVoltageControl": ("transformerVoltageControlMode", "WITH_GENERATOR_VOLTAGE_CONTROL"),
    "ShuntVoltageControl": ("shuntVoltageControlMode", "WITH_GENERATOR_VOLTAGE_CONTROL"),
    "PhaseControl": ("phaseShifterControlMode", "CONTINUOUS_WITH_DISCRETISATION"),
}

# Outer loops with no inline alternative. ``None`` marks a top-level
# Parameters kwarg (handled specially in ``remove_outer_loops`` because of the
# renamed-kwarg fallback); a tuple marks a provider parameter trigger.
_LEGACY_TRIGGER = {
    "DistributedSlack": None,  # top-level `distributed_slack`
    "ReactiveLimits": None,  # top-level `use_reactive_limits` / `no_generator_reactive_limits`
    "VoltageMonitoring": ("svcVoltageMonitoring", "false"),
    "SecondaryVoltageControl": ("secondaryVoltageControl", "false"),
    "AreaInterchangeControl": ("areaInterchangeControl", "false"),
    "AutomationSystem": ("simulateAutomationSystems", "false"),
}

_ALL_LOOP_NAMES = set(_INLINE_MODE) | set(_LEGACY_TRIGGER)


def _copy_parameters(parameters: _lf.Parameters) -> _lf.Parameters:
    """Return an independent copy of ``parameters``.

    Built by copying each constructor-exposed attribute one at a time
    (``getattr``/``setattr``) rather than via ``copy.deepcopy``: ``Parameters``
    is a pybind11-wrapped object backed by a native (Java) instance, and
    whether the generic copy protocol produces a truly independent object --
    as opposed to a new Python wrapper aliasing the same underlying native
    state -- is not guaranteed to hold the same way on every pypowsybl version
    this module supports (1.0.0 to 1.15.0). Explicit getattr/setattr only
    relies on the property interface already used everywhere else here.
    """
    out = _lf.Parameters()
    for name in _PARAM_SIG:
        if name in ("self", "provider_parameters"):
            continue
        try:
            value = getattr(parameters, name)
        except AttributeError:
            continue
        setattr(out, name, value)
    out.provider_parameters = dict(parameters.provider_parameters)
    return out


def remove_outer_loops(
    parameters: _lf.Parameters,
    keep: Iterable[str] = (),
    slack_bus_ids: Optional[Union[str, Iterable[str]]] = None,
    max_outer_loop_iterations: int = 20,
) -> _lf.Parameters:
    """Return a copy of ``parameters`` with OLF's outer loops removed.

    Only the outer-loop mechanism is touched: reactive limits, remote voltage
    control, voltage initialization, balance type, connected-component mode,
    and every other field of ``parameters`` are left exactly as given.

    Parameters
    ----------
    parameters : pypowsybl.loadflow.Parameters
        The parameters to strip outer loops from. Not mutated; a modified
        copy is returned.
    keep : iterable of str, optional
        Names of outer loops to leave active instead of removing. Valid
        names: ``"TransformerVoltageControl"``, ``"ShuntVoltageControl"``,
        ``"PhaseControl"`` (inline-mode family), ``"DistributedSlack"``,
        ``"ReactiveLimits"``, ``"VoltageMonitoring"``,
        ``"SecondaryVoltageControl"``, ``"AreaInterchangeControl"``,
        ``"AutomationSystem"`` (no-inline-alternative family). Keeping a loop
        only stops *this function* from removing it -- for a no-inline-
        alternative loop this function does not force its own creation
        trigger on, so it still runs only if ``parameters`` already enables
        it (e.g. ``distributed_slack=True`` for ``"DistributedSlack"``).
    slack_bus_ids : str or iterable of str, optional
        If given, the slack is pinned by NAME on these bus(es)
        (``slackBusSelectionMode = "NAME"``) and ``read_slack_bus`` is turned
        off. Unrelated to outer-loop removal; a convenience for matching
        lightsim2grid's slack choice.
    max_outer_loop_iterations : int
        Outer-loop iteration budget, only set when ``keep`` is non-empty (so
        a kept loop, e.g. distributed slack, has room to converge). Ignored
        otherwise.

    Returns
    -------
    pypowsybl.loadflow.Parameters
        A new object; ``parameters`` is not modified.
    """
    keep = set(keep)
    unknown = keep - _ALL_LOOP_NAMES
    if unknown:
        raise ValueError(
            f"Unknown outer loop(s) in `keep`: {sorted(unknown)}. "
            f"Valid names: {sorted(_ALL_LOOP_NAMES)}"
        )

    params = _copy_parameters(parameters)

    def put(name, value):
        if name in _PARAM_SIG:
            setattr(params, name, value)
            return True
        return False

    valid = _olf_provider_params()
    prov = dict(params.provider_parameters)

    def put_provider(name, value):
        if valid is not None:
            if name not in valid:
                return  # parameter unknown in this pypowsybl/OLF version
            possible = valid[name]
            if possible is not None and value not in possible:
                return  # value not accepted in this version
        prov[name] = value

    # 1. Inline-mode loops: fold into the inner NR loop unless kept as an
    #    outer loop. Version-independent.
    for loop, (pname, value) in _INLINE_MODE.items():
        if loop not in keep:
            put_provider(pname, value)

    # 2. Loops with no inline alternative.

    # DistributedSlack / ReactiveLimits: top-level flags OLF's DC solve honors
    # directly (see the module docstring) -- forced regardless of pypowsybl
    # version and regardless of the outerLoopNames allow-list below, which
    # only gates the AC outer-loop mechanism.
    if "DistributedSlack" not in keep:
        put("distributed_slack", False)
    if "ReactiveLimits" not in keep:
        if not put("use_reactive_limits", False):
            put("no_generator_reactive_limits", True)

    if _OUTERLOOP_EMPTY_OK:
        # The allow-list alone gates AC outer-loop creation, regardless of any
        # trigger flag, so nothing else needs touching here.
        put_provider("outerLoopNames", ",".join(sorted(keep & set(_LEGACY_TRIGGER))))
        if keep:
            put_provider("maxOuterLoopIterations", str(int(max_outer_loop_iterations)))
    else:
        # Pre-1.4 OLF: no allow-list: the only way to stop a loop is its own
        # creation trigger, forced off unless the caller asked to keep it.
        # DistributedSlack/ReactiveLimits (spec is None) were already handled
        # above, uniformly across versions.
        for loop, spec in _LEGACY_TRIGGER.items():
            if loop in keep or spec is None:
                continue
            pname, value = spec
            put_provider(pname, value)

    if slack_bus_ids is not None:
        if not isinstance(slack_bus_ids, str):
            slack_bus_ids = ",".join(str(s) for s in slack_bus_ids)
        put_provider("slackBusSelectionMode", "NAME")
        put_provider("slackBusesIds", slack_bus_ids)
        put("read_slack_bus", False)

    params.provider_parameters = prov
    return params


def get_pypowsybl_loopfree_parameters(
    slack_bus_ids: Optional[Union[str, Iterable[str]]] = None,
    **overrides,
) -> _lf.Parameters:
    """Build a fresh :class:`pypowsybl.loadflow.Parameters` with every OLF
    outer loop removed (see :func:`remove_outer_loops`).

    Everything other than the outer-loop mechanism is left at the installed
    pypowsybl's own defaults (or at ``**overrides``, if given) -- in
    particular ``voltage_init_mode``, reactive limits and remote voltage
    control are **not** forced, since none of those are outer loops. If a
    test needs a specific value for one of those to be reproducible across
    pypowsybl builds, pass it explicitly via ``**overrides`` (different
    pypowsybl builds are known to ship different defaults for these).

    Parameters
    ----------
    slack_bus_ids : str or iterable of str, optional
        Forwarded to :func:`remove_outer_loops`.
    **overrides
        Any top-level :class:`pypowsybl.loadflow.Parameters` keyword to set
        before removing the outer loops (e.g. ``voltage_init_mode=...``).
        Only applied if the installed pypowsybl accepts that keyword.

    Returns
    -------
    pypowsybl.loadflow.Parameters
        A new object on each call (no shared mutable state).
    """
    kwargs = {}
    for name, value in overrides.items():
        if name in _PARAM_SIG:
            kwargs[name] = value
    return remove_outer_loops(_lf.Parameters(**kwargs), slack_bus_ids=slack_bus_ids)


def get_pypowsybl_loopfree_distributed_slack_parameters(
    slack_bus_ids: Optional[Union[str, Iterable[str]]] = None,
    max_outer_loop_iterations: int = 20,
    **overrides,
) -> _lf.Parameters:
    """Loop-free OLF parameters EXCEPT the active-power slack distribution.

    Identical to :func:`get_pypowsybl_loopfree_parameters` but keeps OLF's
    ``DistributedSlack`` outer loop active, with ``balance_type`` set to
    ``PROPORTIONAL_TO_GENERATION_P_MAX`` -- what lightsim2grid's default
    distributed slack reproduces. Every *other* outer loop is removed.

    Parameters
    ----------
    slack_bus_ids : str or iterable of str, optional
        Forwarded to :func:`remove_outer_loops`.
    max_outer_loop_iterations : int
        Outer-loop iteration cap for the distribution (default 20, the
        upstream OLF default).
    **overrides
        Any top-level :class:`pypowsybl.loadflow.Parameters` keyword to set
        before removing the outer loops. Only applied if the installed
        pypowsybl accepts that keyword.
    """
    overrides.setdefault("distributed_slack", True)
    bt = _enum("BalanceType")
    if bt is not None:
        overrides.setdefault("balance_type", bt.PROPORTIONAL_TO_GENERATION_P_MAX)

    kwargs = {}
    for name, value in overrides.items():
        if name in _PARAM_SIG:
            kwargs[name] = value
    params = _lf.Parameters(**kwargs)
    return remove_outer_loops(
        params,
        keep=("DistributedSlack",),
        slack_bus_ids=slack_bus_ids,
        max_outer_loop_iterations=max_outer_loop_iterations,
    )
