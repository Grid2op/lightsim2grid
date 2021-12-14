// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// containst he help messages for the functions and classes defined in "main" to avoid "pollute" the "main.cpp" file
// with all of this and keep it cleaner

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

    static const std::string KLUSolver;
    static const std::string KLUSolverSingleSlack;
    static const std::string KLUDCSolver;

    static const std::string NICSLUSolver;
    static const std::string NICSLUSolverSingleSlack;
    static const std::string NICSLUDCSolver;

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
    static const std::string id;
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
    static const std::string DataGen;
    static const std::string GenInfo;
    static const std::string is_slack;
    static const std::string slack_weight;

    // specific to sgens
    static const std::string DataSGen;
    static const std::string SGenInfo;

    // specific to loads (and storage units)
    static const std::string DataLoad;
    static const std::string LoadInfo;

    // specific to shunts
    static const std::string DataShunt;
    static const std::string ShuntInfo;

    // specific to transformers
    static const std::string DataTrafo;
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
    static const std::string DataLine;
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

};

struct DocGridModel
{
    static const std::string _internal_do_not_use;
    static const std::string J_description;
    
    static const std::string GridModel;

    static const std::string change_solver;
    static const std::string available_solvers;
    static const std::string get_computation_time;
    static const std::string get_solver_type;
    static const std::string get_solver;

    // accessor
    static const std::string get_lines;
    static const std::string get_trafos;
    static const std::string get_generators;
    static const std::string get_static_generators;
    static const std::string get_shunts;
    static const std::string get_storages;
    static const std::string get_loads;
    
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
