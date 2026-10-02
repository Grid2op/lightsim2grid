import pypowsybl as pypow
import pypowsybl.loadflow as pypow_lf

from lightsim2grid.network import get_pypowsybl_loopfree_parameters


def get_same_slack(case_name):
    if case_name == "ieee9":
        return "VL1_0", 0
    if case_name == "ieee14":
        return "VL1_0", 0
    if case_name == "ieee30":
        return "VL1_0", 0
    if case_name == "ieee57":
        return "VL1_0", 0
    if case_name == "ieee118":
        return "VL69_0", 68
    if case_name == "ieee300":
        return "VL7049_0", 257
    
    raise RuntimeError(f"Unknown env {case_name}")


def get_pypowsybl_parameters(slack_voltage_level=None):
    # single source of truth: the canonical "every outer loop disabled"
    # parameters shipped with lightsim2grid (see
    # lightsim2grid.network.get_pypowsybl_loopfree_parameters). When a slack
    # voltage level is given, the slack is pinned by name so lightsim2grid and
    # pypowsybl use the same slack bus.
    #
    # twt_split_shunt_admittance=True is not an outer loop (it is a Ybus-
    # construction convention: whether a transformer's shunt admittance is
    # split half/half between its two sides), so get_pypowsybl_loopfree_parameters
    # deliberately does not force it -- OLF's own default is False, but
    # lightsim2grid's from_pypowsybl converter always splits it, so it must be
    # forced here for the two engines to solve the same problem (verified on
    # IEEE-300, which has transformers with non-negligible shunt admittance:
    # leaving this at OLF's default desyncs bus voltages by up to ~0.4 pu).
    #
    # The reference is also solved tightly: these tests compare it with a lightsim2grid
    # solve at 1e-6 and below, and OLF's own stopping criterion (newtonRaphsonConvEpsPerEq,
    # 1e-4 per equation by default) leaves angle errors above that. Its default is not
    # forced elsewhere, so a pypowsybl build stopping at it fails these tests for no
    # difference of model.
    params = get_pypowsybl_loopfree_parameters(
        slack_bus_ids=slack_voltage_level,
        twt_split_shunt_admittance=True,
    )
    provider = dict(params.provider_parameters)
    provider["newtonRaphsonConvEpsPerEq"] = "1e-12"
    provider["maxNewtonRaphsonIterations"] = "50"
    params.provider_parameters = provider
    return params

