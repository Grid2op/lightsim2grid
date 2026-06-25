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
    return get_pypowsybl_loopfree_parameters(slack_bus_ids=slack_voltage_level)

