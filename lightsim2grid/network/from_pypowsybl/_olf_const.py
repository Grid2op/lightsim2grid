# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The numerical constants PowSyBl Open Load Flow's own rules are written with, and the
tolerances this package uses to read a reference solve back.

They live here, in one place, rather than next to their first use: several modules
(``_olf_bake``, ``_aux_battery_apc``, ``_aux_add_slack``, ...) apply the same OLF rule to
a different element type, and a constant declared twice is a constant that will one day be
updated once.
"""

# ---------------------------------------------------------------------------------------
# OpenLoadFlow's own constants, mirrored
# ---------------------------------------------------------------------------------------

# OLF's per-unit base (PerUnit.SB, MVA).
_OLF_SB_MVA = 100.0

# OLF's POWER_EPSILON_SI = 1e-4 (AbstractLfGenerator), compared with active powers in MW in
# checkActivePowerControl (slack participation: a zero target, a degenerate range): the
# tolerance (MW) for deciding an active power target is "zero" there.
_ZERO_P_TOL = 1e-4

# The same POWER_EPSILON_SI, but compared with PER-UNIT values in
# checkIfGeneratorStartedForVoltageControl (``getTargetP()`` / ``getMinP()`` are per unit
# there): a generator is "not started" for voltage control below 1e-4 pu, ie this many MW.
_NOT_STARTED_P_TOL_MW = 1e-4 * _OLF_SB_MVA

# OLF's TARGET_V_EPSILON (LfNetworkLoaderImpl), compared with targets in per unit: two units
# of one bus regulating one bus with targets further apart are inconsistent.
_TARGET_V_EPSILON_PU = 1e-2

# OLF's own PlausibleValues.MIN_REACTIVE_RANGE (1 MVar): with the default
# reactiveRangeCheckMode ("MAX"), a generator whose widest reactive range across its
# whole active-power range (or its fixed min_q/max_q box, for a MIN_MAX-kind generator)
# falls below this is discarded from voltage control
# (AbstractLfGenerator.checkIfReactiveRangesAreLargeEnoughForVoltageControl).
_MIN_REACTIVE_RANGE_MVAR = 1.0

# OLF's own defaults (LfNetworkParameters): a generator's targetV, expressed in per
# unit of its (regulated bus) nominal voltage, is discarded from voltage control if
# outside this range -- but only on buses above the nominal-voltage floor below (this
# is minNominalVoltageTargetVoltageCheck, distinct from the *realistic-voltage-check*
# floor that PARAMS_STANDARD sets to 180 kV; this one is left at its own OLF default).
_MIN_PLAUSIBLE_TARGET_V_PU = 0.8
_MAX_PLAUSIBLE_TARGET_V_PU = 1.2
_MIN_NOMINAL_V_FOR_TARGET_V_CHECK_KV = 20.0

# OLF's own plausibleActivePowerLimit default (MW): a unit whose maxP exceeds this is
# discarded from active-power (slack) participation regardless of any other parameter
# (AbstractLfGenerator.checkActivePowerControl).
_MAX_PLAUSIBLE_ACTIVE_POWER_MW = 10000.0

# OLF's hardcoded fallback droop (`AbstractLfGenerator.DEFAULT_DROOP`, "why not"), used
# for every unit whose `activePowerControl` extension does not set its own, under the
# ``PROPORTIONAL_TO_GENERATION_P_MAX`` balance type.
_OLF_DEFAULT_DROOP = 4.0


# ---------------------------------------------------------------------------------------
# Tolerances for reading a reference solve back
# ---------------------------------------------------------------------------------------

# Tolerance (MVAr) for deciding a reactive injection sits "at" its Q limit: the
# larger of a small absolute floor (catches exact/near-exact hits, dominated by
# float rounding) and a fraction of the unit's own Q range (catches OLF's discrete
# reactive-limit outer loop settling a hair inside the limit -- e.g. because P kept
# shifting slightly in later outer-loop iterations after the PV->PQ switch already
# happened). A fixed absolute tolerance alone is either too tight for large-range
# units (missing genuine saturation) or too loose for small-range ones (freezing
# units that still have real headroom).
_Q_LIMIT_TOL_ABS = 1e-3
_Q_LIMIT_TOL_REL = 0.005  # 0.5% of (qmax - qmin)

# Tolerance (pu of the regulated bus nominal voltage) under which the reference
# solve is deemed to have HELD a generator's target voltage. OLF eliminates the
# voltage-magnitude unknown of a PV bus (the magnitude *is* the target), so a
# generator that actually took part in voltage control leaves |V - target_v| at
# machine precision, whatever the Newton-Raphson stopping tolerance; a generator OLF
# did NOT voltage-control -- switched to PQ at a reactive limit, "not started",
# discarded by a consistency check, on a bus below
# ``generatorVoltageControlMinNominalVoltage``, sharing its bus with a local
# controller, ... -- has a free magnitude that lands orders of magnitude further
# away. This threshold sits in the middle of that gap on a log scale, measured on
# real grid snapshots.
_TARGET_V_HELD_TOL_PU = 1e-8

# Tolerance (MVAr) under which a non-regulating generator's realized reactive output is
# deemed to be its target_q: a PQ injection is reproduced to the digit, whatever the
# Newton-Raphson tolerance, well below the smallest clamp seen on real grid snapshots.
_TARGET_Q_TOL_MVAR = 1e-6

# Tolerance (MW) for deciding an angle-droop hvdc line sits "at" its active power limit:
# OLF's ``AcHvdcAcEmulationLimits`` outer loop does not solve the saturated line to a
# tolerance, it imposes the limit as the sending converter's setpoint, so the realized
# flow reproduces the limit to float rounding, while a line still in the linear regime
# lands anywhere along its droop.
_HVDC_P_LIMIT_TOL_MW = 1e-4

# Tolerance (MVAr) under which a generator (or voltage-mode SVC) whose target the reference
# solve held is deemed saturated, for ``bake_saturated_voltage_control``: the relative
# tolerance above is made to catch a unit OLF switched and that settled inside its limit,
# whose voltage is free; on a held unit it would also freeze one OLF still regulates with
# some headroom left, which a contingency may then ask for (or relieve) while the frozen
# unit cannot follow -- on a stiff bus, by a voltage deviation too small to be reported.
_Q_SATURATED_HELD_TOL_MVAR = 0.1

# OLF's own minRealisticVoltage / maxRealisticVoltage defaults (pu) and the margin its
# ReactiveLimits loop applies to them (ReactiveLimitsOuterLoop.REALISTIC_VOLTAGE_MARGIN): in
# its default "robust" remote voltage control mode (voltageRemoteControlRobustMode), a
# controller bus holding a remote target is switched to PQ at its target reactive power as
# soon as its own voltage leaves [min * margin, max / margin].
_OLF_MIN_REALISTIC_VOLTAGE_PU = 0.8
_OLF_MAX_REALISTIC_VOLTAGE_PU = 1.2
_OLF_REALISTIC_VOLTAGE_MARGIN = 1.02
