import argparse
import copy
import time

import numpy as np

import pypowsybl as pypow
import pypowsybl.loadflow as pypow_lf

from lightsim2grid.gridmodel import init_from_pypowsybl
from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP

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
        return "VL7049_0", 257
    
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


def main(case_name, nb_extra_powerflow=100):
    slack_pypowysbl, slack_ls = get_same_slack(case_name)
    
    pypow_grid = getattr(pypow.network, f"create_{case_name}")()
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
    v_res_ls = ls_grid.ac_pf(
        np.ones(ls_grid.get_bus_vn_kv().shape[0], dtype=complex),
        10,
        1e-6)
    end_ls = time.perf_counter()
    ls_grid.unset_changes()
    
    print_configuration(
        pypowbk_error=True, 
        pypowsybl_error=None)
    
    va_rad_ls = np.angle(v_res_ls)
    vm_pu_ls = np.abs(v_res_ls)
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
    print(f"\tLightsim2grid computation time: {1000.*(end_ls - beg_ls):.2e} ms")
    print(f"\tPypowsybl computation time: {1000.*(end_pypow - beg_pypow):.2e} ms")
    
    print(f"For {nb_extra_powerflow} extra powerflows: ")
    
    load_factor = np.linspace(1, 1.02, nb_extra_powerflow, endpoint=False)
    init_loads = pypow_grid.get_loads().copy()
    time_pypow = 0.
    for i in range(nb_extra_powerflow):
        # update the grid
        pypow_grid.update_loads(init_loads[['p0', 'q0']] * load_factor[i])
        # run the powerflow
        beg_pypow = time.perf_counter()
        pypow_lf.run_ac(pypow_grid, parameters=pypowsybl_parameters)
        end_pypow = time.perf_counter()
        time_pypow += end_pypow - beg_pypow
    
    load_p_init, load_q_init, *_ = ls_grid.get_loads_res()
    load_p_init = load_p_init.copy()
    load_q_init = load_q_init.copy()
    v_init_ls = np.ones(ls_grid.get_bus_vn_kv().shape[0], dtype=complex)
    all_loads = np.ones(len(ls_grid.get_loads()), dtype=np.bool_)
    time_ls = 0.
    for i in range(nb_extra_powerflow):
        # update the grid
        new_p = (load_p_init * load_factor[i]).astype(np.float32)
        new_q = (load_q_init * load_factor[i]).astype(np.float32)
        ls_grid.update_loads_p(all_loads, new_p)
        ls_grid.update_loads_q(all_loads, new_q)
        # run the powerflow
        beg_ls = time.perf_counter()
        res_ls = ls_grid.ac_pf(
            v_init_ls,
            10,
            1e-6)
        end_ls = time.perf_counter()
        ls_grid.unset_changes()
        time_ls += end_ls - beg_ls
    print(f"\tLightsim2grid computation time: {1000.*(time_ls / nb_extra_powerflow):.2e} ms / pf")
    print(f"\tPypowsybl computation time: {1000.*(time_pypow / nb_extra_powerflow):.2e} ms / pf")
    
    print("For a contingency anaylisis (results might differ): ")
    # pypowsybl
    analysis = pypow.security.create_analysis()
    pypow_grid.update_loads(init_loads[['p0', 'q0']])
    analysis.add_single_element_contingencies(pypow_grid.get_lines().index)
    analysis.add_single_element_contingencies(pypow_grid.get_2_windings_transformers().index)
    beg_pypow = time.perf_counter()
    res_ca_pypow = analysis.run_ac(pypow_grid, parameters=pypowsybl_parameters)
    end_pypow = time.perf_counter()
    
    # lightsim2grid
    ls_grid.update_loads_p(all_loads, load_p_init.astype(np.float32))
    ls_grid.update_loads_q(all_loads, load_q_init.astype(np.float32))
    ca_ls = ContingencyAnalysisCPP(ls_grid)
    ca_ls.add_all_n1()
    beg_ls = time.perf_counter()
    res_ca_ls = ca_ls.compute(
        v_res_ls,
        10,
        1e-6)
    end_ls = time.perf_counter()
    
    nb_branch = pypow_grid.get_lines().shape[0] + pypow_grid.get_2_windings_transformers().shape[0]
        
    print(f"\tLightsim2grid computation time: {1000.*((end_ls - beg_ls) / nb_branch):.2e} ms / cont")
    print(f"\tPypowsybl computation time: {1000.*((end_pypow - beg_pypow) / nb_branch ):.2e} ms / cont")



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark lightsim2grid with pypowsybl')
    parser.add_argument('--case_name', default=CASE_NAME, type=str,
                        help='ieee case used for the benchmark.')
    args = parser.parse_args()
    
    main(case_name=args.case_name)

    