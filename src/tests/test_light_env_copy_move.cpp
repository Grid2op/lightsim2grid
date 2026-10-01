// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// LightEnv copies and moves. The python tests (test_LightEnv.py) cover the copies against a
// grid2op env; the moves only exist on the C++ side, they are pinned here: nothing is copied
// (the observation buffers keep their address), the observation follows the env it belongs
// to, and a moved-from env refuses to be used instead of dereferencing a null grid.

#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "light_env/light_env.hpp"

using ls2g::AlgorithmType;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::LightEnv;
using ls2g::Protections;
using ls2g::RealVect;
using ls2g::real_type;

namespace {

using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const int N_SUB = 4;
const int N_LINE = 4;
const int N_TS = 6;

// a 4-substation feeder 0-1-2-3 with a second line between 1 and 2, a slack generator at 0,
// a PV generator and a load at 2, a load at 3
LSGrid make_grid()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    const RealVect bus_vn_kv = RealVect::Constant(N_SUB, 138.);
    grid.init_bus(static_cast<unsigned int>(N_SUB), 1u, bus_vn_kv, 0, 0);

    RealVect branch_r = RealVect::Constant(N_LINE, 0.01);
    RealVect branch_x = RealVect::Constant(N_LINE, 0.1);
    const CplxVect branch_h = CplxVect::Zero(N_LINE);
    Eigen::VectorXi from_id(N_LINE), to_id(N_LINE);
    from_id << 0, 1, 2, 1;
    to_id << 1, 2, 3, 2;
    grid.init_powerlines(branch_r, branch_x, branch_h, from_id, to_id);

    RealVect load_p(2), load_q(2);
    load_p << 10., 50.;
    load_q << 2., 10.;
    Eigen::VectorXi load_bus(2);
    load_bus << 2, 3;
    grid.init_loads(load_p, load_q, load_bus);

    RealVect gen_p(2), gen_v(2);
    gen_p << 0., 20.;
    gen_v << 1.04, 1.02;
    Eigen::VectorXi gen_bus(2);
    gen_bus << 0, 2;
    grid.init_generators(gen_p, gen_v, RealVect::Constant(2, -1000.), RealVect::Constant(2, 1000.), gen_bus);
    grid.add_gen_slackbus(0, 1.);

    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    return grid;
}

// an env ready to step: loads growing a little every step, protections that never trip
LightEnv make_env()
{
    LightEnv env(make_grid());
    RealMat load_p(N_TS, 2), load_q(N_TS, 2);
    for(int ts = 0; ts < N_TS; ++ts){
        load_p(ts, 0) = 10. + ts;
        load_p(ts, 1) = 50. + 2. * ts;
        load_q(ts, 0) = 2.;
        load_q(ts, 1) = 10.;
    }
    const real_type nan = std::numeric_limits<real_type>::quiet_NaN();
    env.assign_time_series(load_p, load_q,
                           RealMat::Constant(N_TS, 2, nan), RealMat::Constant(N_TS, 2, nan),
                           RealMat(N_TS, 0), RealMat(N_TS, 0), RealMat(N_TS, 0),
                           RealMat(N_TS, 0), RealMat(N_TS, 0));
    Protections protections;
    RealVect th_lim = RealVect::Constant(N_LINE, 1e9);
    Eigen::VectorXi max_overflow = Eigen::VectorXi::Constant(N_LINE, 99);
    protections.set_thermal_limit_or(th_lim);
    protections.set_thermal_limit_ex(th_lim);
    protections.set_max_line_time_step_overflow(max_overflow);
    env.assign_protections(protections);
    return env;
}

void require_same_state(const LightEnv & a, const LightEnv & b)
{
    REQUIRE(a.get_current_step() == b.get_current_step());
    REQUIRE(a.get_obs().get_current_step() == b.get_obs().get_current_step());
    REQUIRE(a.get_obs().get_p_or() == b.get_obs().get_p_or());
    REQUIRE(a.get_obs().get_a_ex() == b.get_obs().get_a_ex());
    REQUIRE(a.get_obs().get_load_p() == b.get_obs().get_load_p());
    REQUIRE(a.get_obs().get_gen_p() == b.get_obs().get_gen_p());
    REQUIRE(a.get_obs().get_rho() == b.get_obs().get_rho());
}

void require_moved_from(LightEnv & env)
{
    REQUIRE_THROWS_AS(env.reset(), std::logic_error);
    REQUIRE_THROWS_AS(env.step(0), std::logic_error);
    REQUIRE_THROWS_AS(env.get_grid(), std::logic_error);
    REQUIRE_THROWS_AS(env.get_obs().get_load_p(), std::logic_error);
    REQUIRE_THROWS_AS(env.get_obs().get_gen_p(), std::logic_error);
}

}  // namespace

static_assert(std::is_nothrow_move_constructible<LightEnv>::value, "LightEnv should be nothrow move constructible");
static_assert(std::is_nothrow_move_assignable<LightEnv>::value, "LightEnv should be nothrow move assignable");
static_assert(std::is_copy_constructible<LightEnv>::value, "LightEnv should be copy constructible");
static_assert(std::is_copy_assignable<LightEnv>::value, "LightEnv should be copy assignable");

TEST_CASE("A LightEnv copy is an independent env at the same step", "[light_env]")
{
    LightEnv env = make_env();
    env.reset();
    env.step(0);

    LightEnv cpy(env);
    require_same_state(cpy, env);
    REQUIRE(cpy.get_obs().get_p_or().data() != env.get_obs().get_p_or().data());
    REQUIRE(&cpy.get_grid() != &env.get_grid());

    // the copy steps alone, then the original catches up: same result
    const RealVect p_or_before = env.get_obs().get_p_or();
    cpy.step(0);
    REQUIRE(env.get_current_step() == 1);
    REQUIRE(env.get_obs().get_p_or() == p_or_before);
    env.step(0);
    require_same_state(cpy, env);

    // copy assignment, onto an env at another step
    LightEnv other = make_env();
    other.reset();
    other = env;
    require_same_state(other, env);
    other.step(0);
    REQUIRE(env.get_current_step() == 2);
    REQUIRE(other.get_current_step() == 3);
}

TEST_CASE("A LightEnv move copies nothing and its observation follows it", "[light_env]")
{
    LightEnv env = make_env();
    env.reset();
    env.step(0);
    const LightEnv reference(env);

    const real_type * p_or_buffer = env.get_obs().get_p_or().data();
    const real_type * load_p_buffer = env.get_obs().get_load_p().data();
    const LSGrid * grid = &env.get_grid();

    LightEnv moved(std::move(env));
    // the same memory, now owned (and viewed) by `moved`
    REQUIRE(moved.get_obs().get_p_or().data() == p_or_buffer);
    REQUIRE(moved.get_obs().get_load_p().data() == load_p_buffer);
    REQUIRE(&moved.get_grid() == grid);
    require_same_state(moved, reference);
    require_moved_from(env);

    // it keeps stepping as the original would have
    LightEnv reference_next(reference);
    moved.step(0);
    reference_next.step(0);
    require_same_state(moved, reference_next);

    // move assignment
    LightEnv target = make_env();
    target.reset();
    target = std::move(moved);
    REQUIRE(target.get_obs().get_p_or().data() == p_or_buffer);
    REQUIRE(&target.get_grid() == grid);
    require_same_state(target, reference_next);
    require_moved_from(moved);

    // a self move assignment is a no-op
    LightEnv & alias = target;
    target = std::move(alias);
    require_same_state(target, reference_next);
}

TEST_CASE("A moved-from LightEnv is usable again once assigned to", "[light_env]")
{
    LightEnv env = make_env();
    env.reset();
    LightEnv moved(std::move(env));
    require_moved_from(env);

    env = moved;  // copy assignment
    require_same_state(env, moved);
    env.step(0);
    REQUIRE(env.get_current_step() == 1);
    env.reset();
    REQUIRE(env.get_current_step() == 0);

    LightEnv again(std::move(env));
    env = std::move(again);  // move assignment
    env.step(0);
    REQUIRE(env.get_current_step() == 1);
}

TEST_CASE("LightEnvs in a growing std::vector keep their own observation", "[light_env]")
{
    std::vector<LightEnv> envs;
    for(int i = 0; i < 5; ++i){
        envs.push_back(make_env());  // reallocations move (noexcept) the envs already there
        envs.back().reset();
        for(int ts = 0; ts < i; ++ts) envs.back().step(0);
    }
    for(int i = 0; i < 5; ++i){
        REQUIRE(envs[i].get_current_step() == i);
        REQUIRE(envs[i].get_obs().get_current_step() == i);
    }
}
