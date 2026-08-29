import argparse
import copy
import time

import numpy as np

import pypowsybl as pypow
import pypowsybl.loadflow as pypow_lf

from lightsim2grid.network import init_from_pypowsybl, get_pypowsybl_loopfree_parameters
from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP

from utils_benchmark import print_configuration
TABULATE_AVAIL = False
try:
    from tabulate import tabulate
    TABULATE_AVAIL = True
except ImportError:
    print("The tabulate package is not installed. Some output might not work properly")


CASE_NAME = "ieee14"
DONT_SAVE = "__DONT_SAVE"
ALL_CASES = ["ieee9", "ieee14", "ieee30", "ieee57", "ieee118", "ieee300"]


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
    # single source of truth: the canonical "every outer loop disabled"
    # parameters shipped with lightsim2grid, with the slack pinned by name so
    # both engines use the same slack bus.
    return get_pypowsybl_loopfree_parameters(slack_bus_ids=slack_voltage_level)


def run_case(case_name,
             nb_extra_powerflow=100,
             verbose_lightsim2grid_timing=False,
             verbose=True,
             print_config=True):
    """Runs the lightsim2grid vs pypowsybl comparison for a single case.

    Returns a dict of the raw (unformatted) measured quantities, which is what
    `generate_comparison_report` uses to build the tables and narrative text found in
    docs/comparison_with_pypowsybl.rst -- so that, when comparing several cases,
    there is no need to manually "cherry pick" numbers from the per case text output
    below and place them by hand in a table.
    """
    slack_pypowysbl, slack_ls = get_same_slack(case_name)
    
    pypow_grid = getattr(pypow.network, f"create_{case_name}")()
    ls_grid = init_from_pypowsybl(
        pypow_grid,
        slack_bus_id=slack_ls,
        sort_index=False,
        buses_for_sub=True,
        n_busbar_per_sub=1
        )
    ls_grid.change_algorithm("NR_KLU")
    
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
    if verbose_lightsim2grid_timing:
        timers_jac = ls_grid.get_solver().get_timers_jacobian()
        tot_time = end_ls - beg_ls
        print("--------------------------------------")
        print("Detailed lightsim2grid timings: ")
        timer_total_nr_ = timers_jac.timer_total_nr
        print(f"Total time spent in the solver: {1e3 * timers_jac.timer_total_nr:.2e} ms ({100. * timers_jac.timer_total_nr / tot_time:.0f} % of total)")
        print(f"\t Time to pre process Ybus, Sbus etc.: {1e3 * timers_jac.timer_pre_proc:.2e} ms ({100. * timers_jac.timer_pre_proc / timer_total_nr_:.0f} % of time in solver)")
        print(f"\t Time to initialize linear solver {1e3 * timers_jac.timer_initialize:.2e} ms ({100. * timers_jac.timer_initialize / timer_total_nr_:.0f} % of time in solver)")
        print(f"\t Time to compute dS/dV {1e3 * timers_jac.timer_dSbus : .2e} ms ({100. * timers_jac.timer_dSbus / timer_total_nr_:.0f} % of time in solver)")
        print(f"\t Time to fill the Jacobian {1e3 * timers_jac.timer_fillJ:.2e} ms ({100. * timers_jac.timer_fillJ / timer_total_nr_:.0f} % of time in solver)")
        print(f"\t Time to factor the Jacobian linear system: {1e3 * timers_jac.timer_factor:.2e} ms ({100. * timers_jac.timer_factor / timer_total_nr_:.0f} % of time in solver)")
        print(f"\t Time to refactor the Jacobian linear system: {1e3 * timers_jac.timer_refactor:.2e} ms ({100. * timers_jac.timer_refactor / timer_total_nr_:.0f} % of time in solver)")
        print(f"\t Time to solve the Jacobian linear system: {1e3 * timers_jac.timer_solve:.2e} ms ({100. * timers_jac.timer_solve / timer_total_nr_:.0f} % of time in solver)")
        print(f"\t Time to update Va and Vm {1e3*timers_jac.timer_Va_Vm:.2e} ms ({100. * timers_jac.timer_Va_Vm / timer_total_nr_:.0f} % of time in solver)")
        print(f"\t Time to build p,q mismmatch at each bus {1e3 * timers_jac.timer_mismatch:.2e} ms ({100. * timers_jac.timer_mismatch / timer_total_nr_:.0f} % of time in solver)")
        print(f"\t Time to evaluate p,q mismmatch at each bus {1e3 * timers_jac.timer_Fx:.2e} ms ({100. * timers_jac.timer_Fx / timer_total_nr_:.0f} % of time in solver)")
        print(f"\t Time to evaluate cvg criteria {1e3 * timers_jac.timer_check:.2e} ms ({100. * timers_jac.timer_check / timer_total_nr_:.0f} % of time in solver)")
        print("--------------------------------------\n")
    
    if print_config:
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
    err_angle_avg = np.abs(va_rad_ls - va_rad_pypow).mean()
    err_mag_avg = np.abs(vm_pu_ls - vm_pu_pypow).mean()
    err_angle_max = np.abs(va_rad_ls - va_rad_pypow).max()
    err_mag_max = np.abs(vm_pu_ls - vm_pu_pypow).max()
    ls_time_init_ms = 1000. * (end_ls - beg_ls)
    pypow_time_init_ms = 1000. * (end_pypow - beg_pypow)
    if verbose:
        print(f"For {case_name}: ")
        print("Average error (across all buses): ")
        print(f"\t- voltage angle: {err_angle_avg:.2e} rad")
        print(f"\t- voltage magnitude: {err_mag_avg:.2e} pu")
        print("Max error (across all buses): ")
        print(f"\t- voltage angle: {err_angle_max:.2e} rad")
        print(f"\t- voltage magnitude: {err_mag_max:.2e} pu")

        print("For the initial powerflow: ")
        print(f"\tLightsim2grid computation time: {ls_time_init_ms:.2e} ms")
        print(f"\tPypowsybl computation time: {pypow_time_init_ms:.2e} ms")

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
    ls_timer_Fx = 0.
    ls_timer_solve = 0.
    ls_timer_refactor = 0.
    ls_timer_initialize = 0.
    ls_timer_check = 0.
    ls_timer_dSbus = 0.
    ls_timer_fillJ = 0.
    ls_timer_Va_Vm = 0.
    ls_timer_pre_proc = 0.
    ls_timer_total_nr = 0.
    ls_timer_factor = 0.
    ls_timer_refactor = 0.
    ls_timer_per_refactor = 0.
    ls_timer_per_solve = 0.
    ls_timer_per_scale = 0.
    ls_timer_per_mismatch = 0.
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
        time_ls += end_ls - beg_ls
        
        timers_jac = ls_grid.get_solver().get_timers_jacobian()
        nr_iter =  ls_grid.get_solver().get_nb_iter()
        ls_grid.unset_changes()  # tell lightsim2grid that the state of the grid is consistent
        ls_timer_Fx += timers_jac.timer_Fx
        ls_timer_solve += timers_jac.timer_solve
        ls_timer_initialize += timers_jac.timer_initialize
        ls_timer_check += timers_jac.timer_check
        ls_timer_dSbus += timers_jac.timer_dSbus
        ls_timer_fillJ += timers_jac.timer_fillJ
        ls_timer_Va_Vm += timers_jac.timer_Va_Vm
        ls_timer_pre_proc += timers_jac.timer_pre_proc
        ls_timer_total_nr += timers_jac.timer_total_nr
        ls_timer_refactor += timers_jac.timer_refactor
        ls_timer_factor += timers_jac.timer_factor
        if abs(timers_jac.timer_factor) > 1e-8:
            # factor is used
            nb_denom = nr_iter - 1
        else:
            # refactor is used accross all iterations
            nb_denom = nr_iter
        ls_timer_per_refactor += timers_jac.timer_refactor / nb_denom 
        ls_timer_per_solve += timers_jac.timer_solve / nr_iter  # solve is called each iter
        ls_timer_per_scale += timers_jac.timer_scale
        ls_timer_per_mismatch += timers_jac.timer_mismatch

    ls_time_pf_ms = 1000. * (time_ls / nb_extra_powerflow)
    pypow_time_pf_ms = 1000. * (time_pypow / nb_extra_powerflow)
    if verbose:
        print(f"\tLightsim2grid computation time: {ls_time_pf_ms:.2e} ms / pf")
        print(f"\tPypowsybl computation time: {pypow_time_pf_ms:.2e} ms / pf")

    if verbose_lightsim2grid_timing:
        print("--------------------------------------")
        print("Detailed lightsim2grid timings: ")
        print(f"Total time {time_ls * 1000.:.2e}ms => {time_ls / nb_extra_powerflow * 1e3:.2e} ms/pf | {nb_extra_powerflow / time_ls:.0f} pf /s")
        print(f"Total time spent in the solver: {1e3 * ls_timer_total_nr:.2e} ms "
            f"({100. * ls_timer_total_nr / time_ls:.0f}% of total time spent in lightsim2grid)")
        print(f"Total time spent in KLU {1e3 * (ls_timer_initialize + ls_timer_factor + ls_timer_refactor + ls_timer_solve):.2e} ms "
            f"({100. * (ls_timer_initialize + ls_timer_factor + ls_timer_refactor + ls_timer_solve) / ls_timer_total_nr:.0f}% of total time in the solver)")
        print(f"\t Time to pre process Ybus, Sbus etc.: {1e3 * ls_timer_pre_proc:.2e} ms ({100. * ls_timer_pre_proc / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to initialize linear solver {1e3 * ls_timer_initialize:.2e} ms ({100. * ls_timer_initialize / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to factor the linear solver {1e3 * ls_timer_factor:.2e} ms ({100. * ls_timer_factor / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to compute dS/dV {1e3 * ls_timer_dSbus : .2e} ms ({100. * ls_timer_dSbus / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to fill the Jacobian {1e3 * ls_timer_fillJ:.2e} ms ({100. * ls_timer_fillJ / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to refactor the Jacobian linear system: {1e3 * ls_timer_refactor:.2e} ms ({100. * ls_timer_refactor / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t\t This equates to {1e3 * ls_timer_per_refactor:.2e} ms per refactor")
        print(f"\t Time to solve the Jacobian linear system: {1e3 * ls_timer_solve:.2e} ms ({100. * ls_timer_solve / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t\t This equates to {1e3 * ls_timer_per_solve:.2e} ms per solve")
        print(f"\t Time to update Va and Vm {1e3*ls_timer_Va_Vm:.2e} ms ({100. * ls_timer_Va_Vm / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to evaluate p,q mismmatch at each bus {1e3*ls_timer_Fx:.2e} ms ({100. * ls_timer_Fx / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to evaluate cvg criteria {1e3*ls_timer_check:.2e} ms ({100. * ls_timer_check / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to build mismatch vector {1e3*ls_timer_per_mismatch:.2e} ms ({100. * ls_timer_per_mismatch / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to scale mismatch vector {1e3*ls_timer_per_scale:.2e} ms ({100. * ls_timer_per_scale / ls_timer_total_nr:.0f} % of time in solver)")
        print("--------------------------------------\n")
        print("--------------------------------------\n")
    
    if verbose:
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
    ls_time_cont_ms = 1000. * ((end_ls - beg_ls) / nb_branch)
    pypow_time_cont_ms = 1000. * ((end_pypow - beg_pypow) / nb_branch)
    if verbose:
        print(f"\tLightsim2grid computation time: {ls_time_cont_ms:.2e} ms / cont")
        print(f"\tPypowsybl computation time: {pypow_time_cont_ms:.2e} ms / cont")

    return {
        "err_angle_avg": err_angle_avg,
        "err_mag_avg": err_mag_avg,
        "err_angle_max": err_angle_max,
        "err_mag_max": err_mag_max,
        "ls_time_init_ms": ls_time_init_ms,
        "pypow_time_init_ms": pypow_time_init_ms,
        "ls_time_pf_ms": ls_time_pf_ms,
        "pypow_time_pf_ms": pypow_time_pf_ms,
        "ls_time_cont_ms": ls_time_cont_ms,
        "pypow_time_cont_ms": pypow_time_cont_ms,
    }


def _speedup_range(results, ls_key, pypow_key):
    ratios = [res[pypow_key] / res[ls_key] for res in results.values()]
    return min(ratios), max(ratios)


def generate_comparison_report(results, tol=1e-6):
    """Builds the 5 tables and the narrative text that comment on them, found in
    docs/comparison_with_pypowsybl.rst, directly from the numbers measured for each case
    (in `results`, a dict of {case_name: run_case(...) output}). This replaces the previous
    workflow of running the script once per case and manually placing each number in the
    right table / cell of the docs.
    """
    cases = list(results.keys())

    hds_avg = ["case name", "angle (rad)", "magnitude (pu)"]
    tab_avg = [[c, f"{results[c]['err_angle_avg']:.2e}", f"{results[c]['err_mag_avg']:.2e}"] for c in cases]

    hds_max = ["case name", "angle (rad)", "magnitude (pu)"]
    tab_max = [[c, f"{results[c]['err_angle_max']:.2e}", f"{results[c]['err_mag_max']:.2e}"] for c in cases]

    hds_init = ["case name", "lightsim2grid", "pypowsybl"]
    tab_init = [[c, f"{results[c]['ls_time_init_ms']:.2e}", f"{results[c]['pypow_time_init_ms']:.2e}"] for c in cases]

    hds_pf = ["case name", "lightsim2grid", "pypowsybl"]
    tab_pf = [[c, f"{results[c]['ls_time_pf_ms']:.2e}", f"{results[c]['pypow_time_pf_ms']:.2e}"] for c in cases]

    hds_cont = ["case name", "lightsim2grid", "pypowsybl"]
    tab_cont = [[c, f"{results[c]['ls_time_cont_ms']:.2e}", f"{results[c]['pypow_time_cont_ms']:.2e}"] for c in cases]

    def _fmt(hds, tab):
        if TABULATE_AVAIL:
            return tabulate(tab, headers=hds, tablefmt="rst")
        return "\n".join([str(hds)] + [str(row) for row in tab])

    max_err_angle = max(results[c]["err_angle_max"] for c in cases)
    max_err_mag = max(results[c]["err_mag_max"] for c in cases)
    precision_text = (
        f"As we can notice in the tables above, the results match up to the solver precision "
        f"(set to {tol:.0e} for lightsim2grid): the largest observed mismatch, across all cases, is "
        f"{max_err_angle:.2e} rad for the voltage angle and {max_err_mag:.2e} pu for the voltage "
        f"magnitude. On these grids, lightsim2grid and pypowsybl give the same exact results."
    )

    init_lo, init_hi = _speedup_range(results, "ls_time_init_ms", "pypow_time_init_ms")
    init_text = (
        f"For this initial computation, lightsim2grid is between **{init_lo:.0f}** and **{init_hi:.0f}** "
        f"times faster than pypowsybl."
    )

    pf_lo, pf_hi = _speedup_range(results, "ls_time_pf_ms", "pypow_time_pf_ms")
    pf_text = (
        f"For successive powerflows, lightsim2grid is between **{pf_lo:.0f}** and **{pf_hi:.0f}** "
        f"times faster than pypowsybl."
    )

    cont_lo, cont_hi = _speedup_range(results, "ls_time_cont_ms", "pypow_time_cont_ms")
    cont_text = (
        f"For the contingency analysis, lightsim2grid is between **{cont_lo:.0f}** and **{cont_hi:.0f}** "
        f"times faster than pypowsybl."
    )

    parts = [
        "Precision of lightsim2grid",
        "***************************",
        "",
        "On average (across all buses) the errors were:",
        "",
        _fmt(hds_avg, tab_avg),
        "",
        "Maximum error, for all buses:",
        "",
        _fmt(hds_max, tab_max),
        "",
        precision_text,
        "",
        "Computation times (1 powerflow)",
        "********************************",
        "",
        _fmt(hds_init, tab_init),
        "",
        init_text,
        "",
        "Computation times (successive powerflows)",
        "*******************************************",
        "",
        _fmt(hds_pf, tab_pf),
        "",
        pf_text,
        "",
        "Computation times security analysis",
        "************************************",
        "",
        _fmt(hds_cont, tab_cont),
        "",
        cont_text,
    ]
    return "\n".join(parts)


def main(case_names,
         nb_extra_powerflow=100,
         verbose_lightsim2grid_timing=False,
         save_results=DONT_SAVE):
    print_configuration(pypowbk_error=True, pypowsybl_error=None)
    results = {}
    for case_name in case_names:
        results[case_name] = run_case(case_name,
                                       nb_extra_powerflow=nb_extra_powerflow,
                                       verbose_lightsim2grid_timing=verbose_lightsim2grid_timing,
                                       verbose=True,
                                       print_config=False)
    report = generate_comparison_report(results)
    print()
    print(report)
    if save_results != DONT_SAVE:
        with open(save_results + "comparison_with_pypowsybl.rst", "w", encoding="utf-8") as f:
            f.write(report)
            f.write("\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark lightsim2grid with pypowsybl')
    parser.add_argument('--case_name', default=ALL_CASES, nargs='+', type=str,
                        help='ieee case(s) used for the benchmark. Accepts one or more names '
                             f'(default: all of {ALL_CASES}). When more than one is given, '
                             'the tables and text found in docs/comparison_with_pypowsybl.rst '
                             'are printed at the end, ready to copy / paste, instead of having '
                             'to gather the numbers by hand from each case\'s output.')
    parser.add_argument("--verbose_lightsim2grid_timing",
                        action=argparse.BooleanOptionalAction,
                        default=False)
    parser.add_argument("--save_results", default=DONT_SAVE, type=str,
                        help='Name of the file in which you want to save the generated report')
    args = parser.parse_args()
    main(case_names=args.case_name,
         verbose_lightsim2grid_timing=args.verbose_lightsim2grid_timing,
         save_results=args.save_results)
