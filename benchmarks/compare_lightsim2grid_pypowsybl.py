import argparse
import copy
import time

import numpy as np

import pypowsybl as pypow
import pypowsybl.loadflow as pypow_lf

from lightsim2grid.gridmodel import init_from_pypowsybl

from utils_benchmark import print_configuration


CASE_NAME = "ieee14"


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
        return "VL7049_0", 235
    
    raise RuntimeError(f"Unknown env {case_name}")


def get_pypowsybl_parameters(slack_voltage_level):    
    params = pypow_lf.Parameters(
        voltage_init_mode=pypow._pypowsybl.VoltageInitMode.UNIFORM_VALUES,
        transformer_voltage_control_on=False,
        use_reactive_limits=False,
        phase_shifter_regulation_on=False,
        twt_split_shunt_admittance=True,
        shunt_compensator_voltage_control_on=False,
        read_slack_bus=False,
        write_slack_bus=True,
        distributed_slack=False,
        dc_use_transformer_ratio=True,
        hvdc_ac_emulation=False,
        dc_power_factor=1.,
        provider_parameters={
            "useActiveLimits": "false",
            "useReactiveLimits": "false",
            "svcVoltageMonitoring": "false",
            "voltageRemoteControl": "false",
            "writeReferenceTerminals": "false",
            "slackBusSelectionMode" : "NAME",  # for case 118
            "slackBusesIds" : f"{slack_voltage_level}",  # for case 118
            "voltagePerReactivePowerControl": "false",
            "generatorReactivePowerRemoteControl": "false",
            "secondaryVoltageControl": "false",
            }
        )
    return params

def main(case_name):
    slack_pypowysbl, slack_ls = get_same_slack(case_name)
    
    pypow_grid = getattr(pypow.network, f"create_{case_name}")()
    res = pypow_lf.run_ac(pypow_grid)
    ls_grid = init_from_pypowsybl(
        pypow_grid,
        slack_bus_id=slack_ls,
        sort_index=False,
        buses_for_sub=True,
        n_busbar_per_sub=1
        )
    
    pypowsybl_parameters = get_pypowsybl_parameters(slack_pypowysbl)
    beg_pypow = time.perf_counter()
    pypow_lf.run_ac(pypow_grid, parameters=pypowsybl_parameters)
    end_pypow = time.perf_counter()
    
    beg_ls = time.perf_counter()
    res_ls = ls_grid.ac_pf(
        np.ones(ls_grid.get_bus_vn_kv().shape[0], dtype=complex),
        10,
        1e-6)
    end_ls = time.perf_counter()
    
    print_configuration(
        pypowbk_error=True, 
        pypowsybl_error=None)
    
    va_rad_ls = np.angle(res_ls)
    vm_pu_ls = np.abs(res_ls)
    va_rad_pypow = np.deg2rad(pypow_grid.get_buses()["v_angle"].values)
    vm_pu_pypow = (
        pypow_grid.get_buses()["v_mag"].values / 
        pypow_grid.get_voltage_levels().loc[pypow_grid.get_buses()["voltage_level_id"], "nominal_v"].values
        )
    
    print(f"For {case_name}: ")
    print("Average error (across all buses): ")
    print(f"\t- voltage angle: {np.abs(va_rad_ls - va_rad_pypow).mean():.2e} rad")
    print(f"\t- voltage magnitude: {np.abs(vm_pu_ls - vm_pu_pypow).mean():.2e} pu")
    print("Max error (across all buses): ")
    print(f"\t- voltage angle: {np.abs(va_rad_ls - va_rad_pypow).max():.2e} rad")
    print(f"\t- voltage magnitude: {np.abs(vm_pu_ls - vm_pu_pypow).max():.2e} pu")
    
    print("For the initial powerflow: ")
    print(f"Lightsim2grid computation time: {1000.*(end_ls - beg_ls):.2e} ms")
    print(f"Pypowsybl computation time: {1000.*(end_pypow - beg_pypow):.2e} ms")
    
    print("For a contingency anaylisis (results might differ): ")
    print("TODO")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark lightsim2grid with pypowsybl')
    parser.add_argument('--case_name', default=CASE_NAME, type=str,
                        help='ieee case used for the benchmark.')
    args = parser.parse_args()
    
    main(case_name=args.case_name)

    