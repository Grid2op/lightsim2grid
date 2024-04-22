// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
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

#include<string>

struct DocSolver
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
    static const std::string SparseLUSolver;
    static const std::string SparseLUSolverSingleSlack;
    static const std::string DCSolver;
    static const std::string FDPF_XB_SparseLUSolver;
    static const std::string FDPF_BX_SparseLUSolver;

    static const std::string KLUSolver;
    static const std::string KLUSolverSingleSlack;
    static const std::string KLUDCSolver;
    static const std::string FDPF_XB_KLUSolver;
    static const std::string FDPF_BX_KLUSolver;

    static const std::string NICSLUSolver;
    static const std::string NICSLUSolverSingleSlack;
    static const std::string NICSLUDCSolver;
    static const std::string FDPF_XB_NICSLUSolver;
    static const std::string FDPF_BX_NICSLUSolver;

    static const std::string CKTSOSolver;
    static const std::string CKTSOSolverSingleSlack;
    static const std::string CKTSODCSolver;
    static const std::string FDPF_XB_CKTSOSolver;
    static const std::string FDPF_BX_CKTSOSolver;

    static const std::string GaussSeidelSolver;
    static const std::string GaussSeidelSynchSolver;

    // function to select the solver
    static const std::string AnySolver;
    static const std::string get_type;
    static const std::string chooseSolver_get_J_python;
    static const std::string get_computation_time;

};

struct DocIterator
{
    // generic functions
    static const std::string only_avail_res;
    static const std::string id;
    static const std::string name;
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
    static const std::string line_model;
    static const std::string r_pu;
    static const std::string x_pu;
    static const std::string h_pu;

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

    // specific to powerlines
    static const std::string LineContainer;
    static const std::string LineInfo;
    static const std::string bus_or_id;
    static const std::string bus_ex_id;
    static const std::string res_p_or_mw;
    static const std::string res_q_or_mvar;
    static const std::string res_v_or_kv;
    static const std::string res_a_or_ka;
    static const std::string res_p_ex_mw;
    static const std::string res_q_ex_mvar;
    static const std::string res_v_ex_kv;
    static const std::string res_a_ex_ka;
    static const std::string res_theta_or_deg;
    static const std::string res_theta_ex_deg;

    // specific to dc lines
    static const std::string dc_line_formula;
    static const std::string DCLineContainer;
    static const std::string DCLineInfo;
    static const std::string target_p_or_mw;
    static const std::string target_vm_or_pu;
    static const std::string target_vm_ex_pu;
    static const std::string loss_pct;
    static const std::string loss_mw;
    static const std::string res_p_or_mw_dcline;
    static const std::string res_p_ex_mw_dcline;
    static const std::string res_q_or_mvar_dcline;
    static const std::string res_q_ex_mvar_dcline;
    static const std::string res_v_or_kv_dcline;
    static const std::string res_v_ex_kv_dcline;
    static const std::string res_theta_or_deg_dcline;
    static const std::string res_theta_ex_deg_dcline;
    static const std::string gen_or;
    static const std::string gen_ex;
    
};

struct DocGridModel
{
    static const std::string _internal_do_not_use;
    static const std::string J_description;
    
    static const std::string GridModel;

    static const std::string change_solver;
    static const std::string available_solvers;
    static const std::string get_computation_time;
    static const std::string get_dc_computation_time;
    static const std::string get_solver_type;
    static const std::string get_dc_solver_type;
    static const std::string get_solver;
    static const std::string get_dc_solver;

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
    static const std::string id_me_to_ac_solver;
    static const std::string id_ac_solver_to_me;
    static const std::string id_me_to_dc_solver;
    static const std::string id_dc_solver_to_me;
    static const std::string total_bus;
    static const std::string nb_bus;
    static const std::string get_pv;
    static const std::string get_pq;
    static const std::string get_slack_ids;
    static const std::string get_slack_weights;
    static const std::string get_Ybus;
    static const std::string get_dcYbus;
    static const std::string get_Sbus;

    static const std::string check_solution;

    static const std::string deactivate_result_computation;
    static const std::string reactivate_result_computation;
    static const std::string ac_pf;
    static const std::string dc_pf;
};

struct DocComputers
{
    static const std::string Computers;
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

struct DocSecurityAnalysis
{
    static const std::string SecurityAnalysis;

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

#endif  // HELP_FUN_MSG_H