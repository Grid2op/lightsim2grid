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
    static const std::string has_res;
    static const std::string res_p_mw;
    static const std::string res_q_mvar;
    static const std::string res_theta_deg;
    static const std::string res_v_kv;

    // specifi to generators
    static const std::string DataGen;
    static const std::string GenInfo;
    static const std::string is_slack;
    static const std::string slack_weight;
    static const std::string min_q_mvar;
    static const std::string max_q_mvar;

};