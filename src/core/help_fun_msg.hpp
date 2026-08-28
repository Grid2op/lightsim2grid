// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// containst he help messages for the functions and classes defined in "main" to avoid "pollute" the "main.cpp" file
// with all of this and keep it cleaner

#ifndef HELP_FUN_MSG_H
#define HELP_FUN_MSG_H

#include "ls2g_api.hpp"

#include<string>


namespace ls2g {

struct LS2G_API DocSolver
{
    // generic functions
    static const std::string get_J_python;
    static const std::string get_Va;
    static const std::string get_Vm;
    static const std::string get_V;
    static const std::string get_error;
    static const std::string get_nb_iter;
    static const std::string reset;
    static const std::string converged;
    static const std::string compute_pf;
    static const std::string get_timers;

    // TimerJac: returned by get_timers_jacobian()
    static const std::string TimerJac;
    static const std::string timer_Fx;
    static const std::string timer_solve;
    static const std::string timer_factor;
    static const std::string timer_refactor;
    static const std::string timer_initialize;
    static const std::string timer_check;
    static const std::string timer_dSbus;
    static const std::string timer_fillJ;
    static const std::string timer_Va_Vm;
    static const std::string timer_pre_proc;
    static const std::string timer_scale;
    static const std::string timer_mismatch;
    static const std::string timer_total_nr;

    // LinearSolverStats: returned by get_linear_solver_stats() (and the
    // get_linear_solver_stats_bp/bpp FDPF variants)
    static const std::string LinearSolverStats;
    static const std::string nb_reset;
    static const std::string nb_analyze;
    static const std::string nb_factorize;
    static const std::string nb_refactorize;
    static const std::string nb_refactorize_failed;
    static const std::string nb_fallback_factorize;
    static const std::string nb_fallback_factorize_failed;
    static const std::string nb_solve;
    static const std::string get_linear_solver_stats;
    static const std::string get_linear_solver_stats_bp;
    static const std::string get_linear_solver_stats_bpp;

    // NR-only: scaling / refactor policy accessors (bind_nr_algo_policies)
    static const std::string get_scaling_policy_type;
    static const std::string set_scaling_policy;
    static const std::string get_refactor_policy;
    static const std::string set_refactor_policy;
    static const std::string get_max_dVa;
    static const std::string set_max_dVa;
    static const std::string get_max_dVm;
    static const std::string set_max_dVm;
    static const std::string get_ls_c;
    static const std::string set_ls_c;
    static const std::string get_ls_rho;
    static const std::string set_ls_rho;
    static const std::string get_ls_max_iter;
    static const std::string set_ls_max_iter;
    static const std::string get_iw_mu_min;
    static const std::string set_iw_mu_min;
    static const std::string get_iw_mu_max;
    static const std::string set_iw_mu_max;
    static const std::string get_refactor_every_n;
    static const std::string set_refactor_every_n;
    static const std::string get_config;
    static const std::string set_config;
    static const std::string get_theta_to_J_col;
    static const std::string get_vm_to_J_col;
    static const std::string get_q_to_J_col;

    // solver description
    static const std::string NR_SparseLU;
    static const std::string NRSing_SparseLU;
    static const std::string DC_SparseLU;
    static const std::string FDPF_XB_SparseLU;
    static const std::string FDPF_BX_SparseLU;

    static const std::string NR_KLU;
    static const std::string NRSing_KLU;
    static const std::string DC_KLU;
    static const std::string FDPF_XB_KLU;
    static const std::string FDPF_BX_KLU;
    static const std::string NRRefactorRetry_KLU;

    static const std::string NR_NICSLU;
    static const std::string NRSing_NICSLU;
    static const std::string DC_NICSLU;
    static const std::string FDPF_XB_NICSLU;
    static const std::string FDPF_BX_NICSLU;
    static const std::string NRRefactorRetry_NICSLU;

    static const std::string NR_CKTSO;
    static const std::string NRSing_CKTSO;
    static const std::string DC_CKTSO;
    static const std::string FDPF_XB_CKTSO;
    static const std::string FDPF_BX_CKTSO;
    static const std::string NRRefactorRetry_CKTSO;

    static const std::string GaussSeidelAlgo;
    static const std::string GaussSeidelSynchAlgo;

    // function to select the solver
    static const std::string AlgorithmSelector;
    static const std::string get_type;
    static const std::string chooseSolver_get_J_python;
    static const std::string get_computation_time;

};

struct LS2G_API DocIterator
{
    // generic functions
    static const std::string only_avail_res;
    static const std::string id;
    static const std::string name;
    static const std::string sub_id;
    static const std::string pos_topo_vect;
    static const std::string connected;
    static const std::string bus_id;
    static const std::string get_bus_id;
    static const std::string target_p_mw;
    static const std::string target_vm_pu;
    static const std::string target_q_mvar;
    static const std::string has_res;
    static const std::string res_p_mw;
    static const std::string res_q_mvar;
    static const std::string res_theta_deg;
    static const std::string res_v_kv;
    static const std::string min_p_mw;
    static const std::string max_p_mw;
    static const std::string min_q_mvar;
    static const std::string max_q_mvar;
    // shared by any element that can either regulate a bus voltage (PV-like) or
    // apply a fixed reactive setpoint (PQ-like): currently GenInfo and
    // ConverterStationInfo (VSC only)
    static const std::string voltage_regulator_on;
    static const std::string regulated_bus_id;
    static const std::string line_model;
    static const std::string r_pu;
    static const std::string x_pu;
    static const std::string h_pu;

    // specific to substations
    static const std::string SubstationContainer;
    static const std::string SubstationInfo;
    static const std::string nb_max_busbars;
    static const std::string vn_kv;

    // specific to generators
    static const std::string GeneratorContainer;
    static const std::string GenInfo;
    static const std::string is_slack;
    static const std::string slack_weight;

    // specific to sgens
    static const std::string SGenContainer;
    static const std::string SGenInfo;

    // specific to SVC (Static Var Compensators)
    static const std::string SvcContainer;
    static const std::string SvcInfo;
    static const std::string RegulationMode;
    static const std::string regulation_mode;
    static const std::string svc_target_vm_pu;
    static const std::string svc_target_q_mvar;
    static const std::string slope_pu;
    static const std::string b_min;
    static const std::string b_max;
    static const std::string svc_regulated_bus_id;

    // specific to loads (and storage units)
    static const std::string LoadContainer;
    static const std::string LoadInfo;

    // specific to shunts
    static const std::string ShuntContainer;
    static const std::string ShuntInfo;

    // specific to transformers
    static const std::string TrafoContainer;
    static const std::string TrafoInfo;
    static const std::string bus_hv_id;
    static const std::string bus_lv_id;
    static const std::string is_tap_hv_side;
    static const std::string ratio;
    static const std::string shift_rad;
    static const std::string res_p_hv_mw;
    static const std::string res_q_hv_mvar;
    static const std::string res_v_hv_kv;
    static const std::string res_a_hv_ka;
    static const std::string res_p_lv_mw;
    static const std::string res_q_lv_mvar;
    static const std::string res_v_lv_kv;
    static const std::string res_a_lv_ka;
    static const std::string res_theta_hv_deg;
    static const std::string res_theta_lv_deg;

    // specific to powerlines (also used, where noted, for hvdc lines' bus ids)
    static const std::string LineContainer;
    static const std::string LineInfo;
    static const std::string get_bus_id_side_1;
    static const std::string get_bus_id_side_2;
    static const std::string res_p_1_mw;
    static const std::string res_q_1_mvar;
    static const std::string res_v_1_kv;
    static const std::string res_a_1_ka;
    static const std::string res_p_2_mw;
    static const std::string res_q_2_mvar;
    static const std::string res_v_2_kv;
    static const std::string res_a_2_ka;
    static const std::string res_theta_1_deg;
    static const std::string res_theta_2_deg;

    // shared by LineInfo and TrafoInfo (the two-port pi-model admittance
    // matrix, see `line_model`)
    static const std::string yac_11;
    static const std::string yac_12;
    static const std::string yac_21;
    static const std::string yac_22;
    static const std::string yac_eff_11;
    static const std::string yac_eff_12;
    static const std::string yac_eff_21;
    static const std::string yac_eff_22;
    static const std::string get_yac_eff_11;
    static const std::string get_yac_eff_12;
    static const std::string get_yac_eff_21;
    static const std::string get_yac_eff_22;
    static const std::string ydc_11;
    static const std::string ydc_12;
    static const std::string ydc_21;
    static const std::string ydc_22;

    // specific to hvdc lines (kept under the "dc line" name for backward compatibility)
    static const std::string dc_line_formula;
    static const std::string DCLineContainer;
    static const std::string DCLineInfo;
    static const std::string target_p_1_mw_dcline;
    static const std::string target_p_2_mw_dcline;
    static const std::string target_vm_1_pu_dcline;
    static const std::string target_vm_2_pu_dcline;
    static const std::string loss_pct;
    static const std::string loss_mw;
    static const std::string res_p_1_mw_dcline;
    static const std::string res_p_2_mw_dcline;
    static const std::string res_q_1_mvar_dcline;
    static const std::string res_q_2_mvar_dcline;
    static const std::string res_v_1_kv_dcline;
    static const std::string res_v_2_kv_dcline;
    static const std::string res_theta_1_deg_dcline;
    static const std::string res_theta_2_deg_dcline;

    // IIDM model (powsybl-style hvdc lines)
    static const std::string converters_mode;
    static const std::string p_setpoint_mw;
    static const std::string r_ohm;
    static const std::string nominal_v_kv;

    // angle-droop ("AC emulation")
    static const std::string droop_enabled;
    static const std::string droop_p0_mw;
    static const std::string droop_k_mw_per_rad;
    static const std::string pmax_1to2_mw;
    static const std::string pmax_2to1_mw;
    static const std::string status_droop;

    // hvdc converter stations (one per side of a hvdc line)
    static const std::string station_side_1;
    static const std::string station_side_2;
    static const std::string ConverterStationInfo;
    static const std::string converter_type;
    static const std::string loss_factor;
    static const std::string power_factor;

    // per-container specific variants of the shared fields above, cross-referencing the
    // LSGrid getter/setter that reads / writes them (added so a bare `bus_id` etc. cannot be
    // read from the field alone, the corresponding `LSGrid` method is always one click away)
    static const std::string load_bus_id;
    static const std::string gen_bus_id;
    static const std::string shunt_bus_id;
    static const std::string sgen_bus_id;
    static const std::string storage_bus_id;
    static const std::string svc_bus_id;
    static const std::string line_bus1_id;
    static const std::string line_bus2_id;
    static const std::string hvdc_bus1_id;
    static const std::string hvdc_bus2_id;

    static const std::string load_connected;
    static const std::string gen_connected;
    static const std::string shunt_connected;
    static const std::string sgen_connected;
    static const std::string storage_connected;
    static const std::string svc_connected;
    static const std::string line_connected_global;
    static const std::string trafo_connected_global;
    static const std::string hvdc_connected_global;
    static const std::string line_connected1;
    static const std::string line_connected2;
    static const std::string trafo_connected1;
    static const std::string trafo_connected2;
    static const std::string hvdc_connected1;
    static const std::string hvdc_connected2;

    static const std::string load_target_p_mw;
    static const std::string gen_target_p_mw;
    static const std::string shunt_target_p_mw;
    static const std::string sgen_target_p_mw;
    static const std::string storage_target_p_mw;

    static const std::string load_target_q_mvar;
    static const std::string shunt_target_q_mvar;
    static const std::string sgen_target_q_mvar;
    static const std::string storage_target_q_mvar;

    static const std::string gen_target_vm_pu;

    static const std::string load_pos_topo_vect;
    static const std::string gen_pos_topo_vect;
    static const std::string shunt_pos_topo_vect;
    static const std::string sgen_pos_topo_vect;
    static const std::string storage_pos_topo_vect;
    static const std::string svc_pos_topo_vect;
    static const std::string line_pos1_topo_vect;
    static const std::string line_pos2_topo_vect;
    static const std::string trafo_pos1_topo_vect;
    static const std::string trafo_pos2_topo_vect;
    static const std::string hvdc_pos1_topo_vect;
    static const std::string hvdc_pos2_topo_vect;

    static const std::string load_sub_id;
    static const std::string gen_sub_id;
    static const std::string shunt_sub_id;
    static const std::string sgen_sub_id;
    static const std::string storage_sub_id;
    static const std::string svc_sub_id;
    static const std::string line_sub1_id;
    static const std::string line_sub2_id;
    static const std::string trafo_sub1_id;
    static const std::string trafo_sub2_id;
    static const std::string hvdc_sub1_id;
    static const std::string hvdc_sub2_id;

    static const std::string substation_name;
};

struct LS2G_API DocLSGrid
{
    static const std::string _internal_do_not_use;
    static const std::string J_description;
    
    static const std::string LSGrid;

    static const std::string check_grid;

    static const std::string change_algorithm;
    static const std::string available_algorithm_names;
    static const std::string available_default_algorithms;
    static const std::string get_computation_time;
    static const std::string get_dc_computation_time;
    static const std::string get_algo_type;
    static const std::string get_dc_algo_type;
    static const std::string get_algo;
    static const std::string get_dc_algo;

    // accessor
    static const std::string get_lines;
    static const std::string get_trafos;
    static const std::string get_generators;
    static const std::string get_static_generators;
    static const std::string get_shunts;
    static const std::string get_storages;
    static const std::string get_loads;
    static const std::string get_dclines;
    static const std::string get_substations;

    // solver-family change tracking (see DocMisc::AlgoControl)
    static const std::string get_ac_algo_controler;
    static const std::string get_dc_algo_controler;

    // internal bookkeeping properties (still user-visible, eg for pickling / debugging)
    static const std::string _ls_to_orig;
    static const std::string _orig_to_ls;
    static const std::string _max_nb_bus_per_sub;
    static const std::string _init_kwargs;
    static const std::string _bus_fusion_rep;

    // algo config (scaling/refactor policy params), see DocMisc::AlgoConfig
    static const std::string get_ac_algo_config;
    static const std::string set_ac_algo_config;
    static const std::string get_dc_algo_config;
    static const std::string set_dc_algo_config;
    static const std::string change_algorithm_by_name;

    // per-bus operating voltage limits
    static const std::string set_bus_voltage_limits;
    static const std::string get_bus_vmin_kv;
    static const std::string get_bus_vmax_kv;

    static const std::string get_svcs;

    // slack / pv-pq bookkeeping
    static const std::string turnedoff_no_pv;
    static const std::string turnedoff_pv;
    static const std::string set_reference_slack_bus;
    static const std::string get_reference_slack_bus;
    static const std::string set_ignore_status_global;
    static const std::string set_synch_status_both_side;

    // names
    static const std::string get_line_names;
    static const std::string get_trafo_names;
    static const std::string set_line_current_limit_side1;
    static const std::string set_line_current_limit_side2;
    static const std::string set_trafo_current_limit_side1;
    static const std::string set_trafo_current_limit_side2;

    // half-open (per-side) connect / disconnect
    static const std::string deactivate_powerline_side1;
    static const std::string deactivate_powerline_side2;
    static const std::string reactivate_powerline_side1;
    static const std::string reactivate_powerline_side2;
    static const std::string deactivate_trafo_side1;
    static const std::string deactivate_trafo_side2;
    static const std::string reactivate_trafo_side1;
    static const std::string reactivate_trafo_side2;
    static const std::string get_lines_status_side1;
    static const std::string get_lines_status_side2;
    static const std::string get_trafo_status_side1;
    static const std::string get_trafo_status_side2;
    static const std::string deactivate_dcline_side1;
    static const std::string deactivate_dcline_side2;

    // transformer tap / phase-shift
    static const std::string change_shift_trafo;
    static const std::string change_shift_trafo_deg;
    static const std::string set_trafo_shift_dependent_rx;

    // remote voltage control (generators)
    static const std::string set_gen_regulated_bus;

    // hvdc angle-droop regime
    static const std::string set_status_droop_hvdc;
    static const std::string get_status_droop_hvdc;
    static const std::string get_status_droop_hvdc_vect;

    // solver-internal state, for external solvers re-deriving the NR system
    // (distributed slack / VoltageControl / bordered-block bookkeeping)
    static const std::string get_slack_col_solver;
    static const std::string get_slack_absorbed_solver;
    static const std::string get_controller_q_solver;
    static const std::string get_controller_kind_solver;
    static const std::string get_controller_elem_id_solver;
    static const std::string get_controller_q_col_solver;
    static const std::string get_p_buses_solver;
    static const std::string get_p_rows_solver;
    static const std::string get_q_buses_solver;
    static const std::string get_q_rows_solver;
    static const std::string get_theta_buses_solver;
    static const std::string get_theta_cols_solver;
    static const std::string get_vm_buses_solver;
    static const std::string get_vm_cols_solver;
    static const std::string get_hvdc_droop_data_solver;

    // misc grid-level bookkeeping (previously "TODO")
    static const std::string copy;
    static const std::string timer_last_ac_pf;
    static const std::string timer_last_dc_pf;
    static const std::string get_turnedoff_gen_pv;
    static const std::string update_slack_weights;
    static const std::string update_slack_weights_by_id;
    static const std::string assign_slack_to_most_connected;
    static const std::string consider_only_main_component;
    static const std::string get_ignore_status_global;
    static const std::string get_synch_status_both_side;
    static const std::string set_line_names;
    static const std::string set_dcline_names;
    static const std::string set_trafo_names;
    static const std::string set_gen_names;
    static const std::string set_load_names;
    static const std::string set_storage_names;
    static const std::string set_sgen_names;
    static const std::string set_shunt_names;
    static const std::string set_svc_names;
    static const std::string change_ratio_trafo;

    // retrieve the results
    static const std::string get_J_python;
    static const std::string get_Va;
    static const std::string get_Vm;
    static const std::string get_V;
    static const std::string get_J_python_solver;
    static const std::string get_Va_solver;
    static const std::string get_Vm_solver;
    static const std::string get_V_solver;
    static const std::string id_me_to_ac_solver;
    static const std::string id_ac_solver_to_me;
    static const std::string id_me_to_dc_solver;
    static const std::string id_dc_solver_to_me;
    static const std::string total_bus;
    static const std::string nb_connected_bus;

    // TODO doc: make more precise when things are copied and when things are not
    static const std::string get_pv;
    static const std::string get_pq;
    static const std::string get_slack_ids;
    static const std::string get_slack_ids_dc;
    static const std::string get_slack_weights;
    static const std::string get_pv_solver;
    static const std::string get_pq_solver;
    static const std::string get_ac_pv_solver;
    static const std::string get_dc_pv_solver;
    static const std::string get_ac_pq_solver;
    static const std::string get_dc_pq_solver;
    static const std::string get_slack_ids_solver;
    static const std::string get_slack_ids_dc_solver;
    static const std::string get_slack_weights_solver;
    static const std::string get_ac_slack_weights_solver;
    static const std::string get_dc_slack_weights_solver;
    static const std::string get_Ybus;
    static const std::string get_dcYbus;
    static const std::string get_Sbus;
    static const std::string get_dcSbus;
    static const std::string get_Ybus_solver;
    static const std::string get_dcYbus_solver;
    static const std::string get_Sbus_solver;
    static const std::string get_dcSbus_solver;
    static const std::string get_ptdf;
    static const std::string get_ptdf_solver;
    static const std::string get_lodf;
    static const std::string get_Bf;
    static const std::string get_Bf_solver;

    static const std::string check_solution;

    static const std::string deactivate_result_computation;
    static const std::string reactivate_result_computation;
    static const std::string ac_pf;
    static const std::string dc_pf;

    // grid-level scalars / bus infrastructure (no per-element python class to reference)
    static const std::string get_bus_vn_kv;
    static const std::string get_bus_status;
    static const std::string set_init_vm_pu;
    static const std::string get_init_vm_pu;
    static const std::string set_sn_mva;
    static const std::string get_sn_mva;
    static const std::string set_n_sub;
    static const std::string get_n_sub;
    static const std::string set_max_nb_bus_per_sub;

    // solver-state bookkeeping
    static const std::string unset_changes;
    static const std::string tell_solver_need_reset;
    static const std::string allow_ac_cache_reuse;
    static const std::string allow_dc_cache_reuse;
    static const std::string allow_cache_reuse;
    static const std::string get_allow_ac_cache_reuse;
    static const std::string get_allow_dc_cache_reuse;
    static const std::string get_allow_cache_reuse;
    static const std::string prevent_ac_cache_reuse;
    static const std::string prevent_dc_cache_reuse;
    static const std::string prevent_cache_reuse;

    // slack designation
    static const std::string add_gen_slackbus;
    static const std::string remove_gen_slackbus;

    // substation names (bulk)
    static const std::string set_substation_names;
    static const std::string get_substation_names;

    // per-element bus getters
    static const std::string get_bus_load;
    static const std::string get_bus_gen;
    static const std::string get_bus_shunt;
    static const std::string get_bus_sgen;
    static const std::string get_bus_storage;
    static const std::string get_bus_svc;
    static const std::string get_bus1_powerline;
    static const std::string get_bus2_powerline;
    static const std::string get_bus1_trafo;
    static const std::string get_bus2_trafo;
    static const std::string get_bus1_dcline;
    static const std::string get_bus2_dcline;

    // per-element bus setters (topology change)
    static const std::string change_bus_load;
    static const std::string change_bus_gen;
    static const std::string change_bus_svc;
    static const std::string change_bus_shunt;
    static const std::string change_bus_sgen;
    static const std::string change_bus_storage;
    static const std::string change_bus1_powerline;
    static const std::string change_bus2_powerline;
    static const std::string change_bus1_trafo;
    static const std::string change_bus2_trafo;
    static const std::string change_bus1_dcline;
    static const std::string change_bus2_dcline;

    // activation status
    static const std::string deactivate_powerline;
    static const std::string reactivate_powerline;
    static const std::string deactivate_trafo;
    static const std::string reactivate_trafo;
    static const std::string deactivate_load;
    static const std::string reactivate_load;
    static const std::string deactivate_gen;
    static const std::string reactivate_gen;
    static const std::string deactivate_svc;
    static const std::string reactivate_svc;
    static const std::string deactivate_shunt;
    static const std::string reactivate_shunt;
    static const std::string deactivate_sgen;
    static const std::string reactivate_sgen;
    static const std::string deactivate_storage;
    static const std::string reactivate_storage;
    static const std::string deactivate_dcline;
    static const std::string reactivate_dcline;

    // setpoint setters
    static const std::string change_p_load;
    static const std::string change_q_load;
    static const std::string change_p_gen;
    static const std::string change_v_gen;
    static const std::string change_p_shunt;
    static const std::string change_q_shunt;
    static const std::string change_p_sgen;
    static const std::string change_q_sgen;
    static const std::string change_p_storage;
    static const std::string change_q_storage;
    static const std::string change_p_dcline;
    static const std::string change_v1_dcline;
    static const std::string change_v2_dcline;

    // bulk vector getters (whole-container equivalents of the per-element attributes above)
    static const std::string get_loads_res;
    static const std::string get_shunts_res;
    static const std::string get_gen_res;
    static const std::string get_storages_res;
    static const std::string get_sgens_res;
    static const std::string get_loads_status;
    static const std::string get_shunts_status;
    static const std::string get_gen_status;
    static const std::string get_storages_status;
    static const std::string get_sgens_status;
    static const std::string get_lines_status;
    static const std::string get_trafo_status;
    static const std::string get_line_res1;
    static const std::string get_line_res2;
    static const std::string get_trafo_res1;
    static const std::string get_trafo_res2;
    static const std::string get_gen_theta;
    static const std::string get_load_theta;
    static const std::string get_shunt_theta;
    static const std::string get_storage_theta;
    static const std::string get_line_theta1;
    static const std::string get_line_theta2;
    static const std::string get_trafo_theta1;
    static const std::string get_trafo_theta2;
    static const std::string get_all_shunt_buses;
    static const std::string get_loads_res_full;
    static const std::string get_shunts_res_full;
    static const std::string get_gen_res_full;
    static const std::string get_storages_res_full;
    static const std::string get_sgens_res_full;
    static const std::string get_line_res1_full;
    static const std::string get_line_res2_full;
    static const std::string get_trafo_res1_full;
    static const std::string get_trafo_res2_full;
    static const std::string get_dcline_res1_full;
    static const std::string get_dcline_res2_full;
    static const std::string get_shunt_target_p;
    static const std::string get_load_target_p;
    static const std::string get_gen_target_p;
    static const std::string get_sgen_target_p;
    static const std::string get_storage_target_p;

    // bulk vectorized setters (grid2op-backend fast path)
    static const std::string update_gens_p;
    static const std::string update_sgens_p;
    static const std::string update_gens_v;
    static const std::string update_loads_p;
    static const std::string update_loads_q;
    static const std::string update_storages_p;
    static const std::string update_topo;

    // structural position / substation-id setters (one-time loader setup)
    static const std::string set_load_pos_topo_vect;
    static const std::string set_gen_pos_topo_vect;
    static const std::string set_line_pos1_topo_vect;
    static const std::string set_line_pos2_topo_vect;
    static const std::string set_trafo_pos1_topo_vect;
    static const std::string set_trafo_pos2_topo_vect;
    static const std::string set_storage_pos_topo_vect;
    static const std::string set_load_to_subid;
    static const std::string set_gen_to_subid;
    static const std::string set_shunt_to_subid;
    static const std::string set_line_to_sub1_id;
    static const std::string set_line_to_sub2_id;
    static const std::string set_trafo_to_sub1_id;
    static const std::string set_trafo_to_sub2_id;
    static const std::string set_storage_to_subid;

    // bulk container constructors
    static const std::string init_powerlines;
    static const std::string init_powerlines_full;
    static const std::string init_shunt;
    static const std::string init_trafo_pandapower;
    static const std::string init_trafo;
    static const std::string init_generators;
    static const std::string init_generators_full;
    static const std::string init_loads;
    static const std::string init_storages;
    static const std::string init_sgens;
    static const std::string init_dclines;
    static const std::string init_hvdc_lines;
    static const std::string init_svcs;
};

struct LS2G_API DocTimeSeries
{
    static const std::string TimeSeries;
    static const std::string total_time;
    static const std::string solver_time;
    static const std::string preprocessing_time;
    static const std::string amps_computation_time;
    static const std::string nb_solved;
    static const std::string nb_converged;
    static const std::string converged_mask;
    static const std::string get_status;

    static const std::string compute_Vs;
    static const std::string compute_flows;
    static const std::string compute_power_flows;

    static const std::string get_flows;
    static const std::string get_power_flows;
    static const std::string get_voltages;
    static const std::string get_sbuses;
    static const std::string clear;

    // shared by TimeSeriesCPP and InjectionSweepCPP (see DocInjectionSweep)
    static const std::string nb_thread;
    static const std::string thread_init_time;
    static const std::string init_from_n_powerflow;
};

// InjectionSweepCPP has the very same interface as TimeSeriesCPP, so it reuses every
// DocTimeSeries member above; only what differs is redefined here.
struct LS2G_API DocInjectionSweep
{
    static const std::string InjectionSweep;
    static const std::string init_from_n_powerflow;
};

struct LS2G_API DocContingencyAnalysis
{
    static const std::string ContingencyAnalysis;

    static const std::string violation_threshold;

    static const std::string preprocessing_time;
    static const std::string modif_Ybus_time;

    static const std::string add_all_n1;
    static const std::string add_n1;
    static const std::string add_nk;
    static const std::string add_multiple_n1;
    static const std::string clear;
    static const std::string remove_n1;
    static const std::string remove_nk;
    static const std::string remove_multiple_n1;

    static const std::string my_defaults_vect;

    static const std::string compute;
    static const std::string compute_flows;
    static const std::string compute_power_flows;

    static const std::string get_flows;
    static const std::string get_voltages;
    static const std::string get_power_flows;

    // limit-violation results (batch_algorithm/LimitViolation.hpp)
    static const std::string ViolationElementType;
    static const std::string LimitViolationType;
    static const std::string LimitViolation;
    static const std::string element_type;
    static const std::string element_id;
    static const std::string side;
    static const std::string violation_type;
    static const std::string value;
    static const std::string limit;
    static const std::string violation_name;
};

struct LS2G_API DocMisc
{
    // PandaPowerConverter: pandapower per-unit data conversion, used by the
    // pandapower grid loader (network/from_pandapower)
    static const std::string PandaPowerConverter;
    static const std::string set_f_hz;
    static const std::string set_sn_mva;
    static const std::string get_line_param_legacy;
    static const std::string get_line_param;
    static const std::string get_trafo_param_pp3;
    static const std::string get_trafo_param_pp2;

    // AlgoControl: per-solver-family (AC / DC) change-tracking flags, read by
    // LSGrid.get_ac_algo_controler() / get_dc_algo_controler()
    static const std::string AlgoControl;
    static const std::string has_dimension_changed;
    static const std::string has_pv_changed;
    static const std::string has_pq_changed;
    static const std::string has_slack_participate_changed;
    static const std::string need_reset_solver;
    static const std::string need_recompute_sbus;
    static const std::string need_recompute_ybus;
    static const std::string ybus_change_sparsity_pattern;
    static const std::string has_slack_weight_changed;
    static const std::string has_v_changed;
    static const std::string has_ybus_some_coeffs_zero;
    static const std::string has_one_el_changed_bus;

    // AlgoConfig: serializable NR scaling/refactor policy configuration blob,
    // read/set via get_config()/set_config() on any NR-based solver
    static const std::string AlgoConfig;
    static const std::string int_params;
    static const std::string real_params;
};


} // namespace ls2g

#endif  // HELP_FUN_MSG_H
