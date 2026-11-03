# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""OpenLoadFlow parameters that do not depend on the pypowsybl build that runs them.

pypowsybl builds of one version number ship different OpenLoadFlow defaults (the realistic
voltage range, the reactive-limit loading rules, the transformer shunt split, the start...).
lightsim2grid's OpenLoadFlow rules (``OlfLoadingParameters``, the outer loops) follow one of
them, the reference build; a test comparing with OpenLoadFlow pins every parameter on which
the builds differ to that build's value, so that it compares like with like wherever it runs.

This module is intentionally NOT named ``test_*`` so that ``unittest`` discovery ignores it.
"""

import pypowsybl.loadflow as lf

#: the reference build's values of the provider parameters the builds disagree on
REFERENCE_PROVIDER_PARAMETERS = {
    "disableInconsistentVoltageControls": "true",
    "extrapolateReactiveLimits": "true",
    "forceTargetQInReactiveLimits": "true",
    "generatorVoltageControlMinNominalVoltage": "120.0",
    "maxNewtonRaphsonIterations": "30",
    "maxOuterLoopIterations": "30",
    "minRealisticVoltage": "0.8",
    "maxRealisticVoltage": "1.2",
    "minNominalVoltageRealisticVoltageCheck": "180.0",
    "stateVectorScalingMode": "MAX_VOLTAGE_CHANGE",
    "maxVoltageChangeStateVectorScalingMaxDv": "0.4",
    "maxVoltageChangeStateVectorScalingMaxDphi": "1.0",
    "transformerVoltageControlMode": "AFTER_GENERATOR_VOLTAGE_CONTROL",
    "transformerVoltageControlUseInitialTapPosition": "true",
}


def reference_parameters(provider=None, **kwargs):
    """``pypowsybl.loadflow.Parameters`` with the reference build's values wherever builds
    differ, then ``kwargs`` (top-level parameters) and ``provider`` (provider parameters, as
    strings) on top."""
    top = dict(twt_split_shunt_admittance=True, voltage_init_mode=lf.VoltageInitMode.DC_VALUES,
               component_mode=lf.ComponentMode.MAIN_SYNCHRONOUS)
    top.update(kwargs)
    params = lf.Parameters(**top)
    prov = dict(REFERENCE_PROVIDER_PARAMETERS)
    prov.update(provider or {})
    params.provider_parameters = prov
    return params
