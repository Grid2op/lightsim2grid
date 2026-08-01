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
    static const std::string bus_1_id;
    static const std::string bus_2_id;
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
    static const std::string get_slack_ids_solver;
    static const std::string get_slack_ids_dc_solver;
    static const std::string get_slack_weights_solver;
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
};

struct LS2G_API DocTimeSeries
{
    static const std::string TimeSeries;
    static const std::string total_time;
    static const std::string solver_time;
    static const std::string preprocessing_time;
    static const std::string amps_computation_time;
    static const std::string nb_solved;
    static const std::string get_status;

    static const std::string compute_Vs;
    static const std::string compute_flows;
    static const std::string compute_power_flows;

    static const std::string get_flows;
    static const std::string get_power_flows;
    static const std::string get_voltages;
    static const std::string get_sbuses;
    static const std::string clear;
};

struct LS2G_API DocContingencyAnalysis
{
    static const std::string ContingencyAnalysis;

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
};


} // namespace ls2g

#endif  // HELP_FUN_MSG_H
