// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Tests of the outer-loop driver of the NROuter_* algorithms (OuterLoopAlgo): the order in
// which OpenLoadFlow's engine runs, re-runs and stops the loops, and the one-analyze
// premise. The loops here are scripted stand-ins that say STABLE / UNSTABLE / FAILED on cue
// and log every call; the real loops have their own tests. C++14 only (project policy).

#include <complex>
#include <memory>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "powerflow_algorithm/outer_loop/BaseOuterLoop.hpp"

using Catch::Approx;
using ls2g::BaseOuterLoop;
using ls2g::CplxVect;
using ls2g::ErrorType;
using ls2g::LSGrid;
using ls2g::LimitViolation;
using ls2g::OuterContext;
using ls2g::OuterLoopStats;
using ls2g::OuterLoopStatus;
using ls2g::RealVect;
using ls2g::real_type;

namespace {

// radial feeder 0-1-2, slack generator at bus 0, 50 MW / 10 MVar load at bus 2; with
// `pv_gen`, a second generator (10 MW, PV at 1.01 pu) at bus 1
LSGrid make_grid(bool pv_gen = false)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(3, 1, RealVect::Constant(3, 138.), 0, 0);
    Eigen::VectorXi from_id(2), to_id(2);
    from_id << 0, 1;
    to_id << 1, 2;
    grid.init_powerlines(RealVect::Constant(2, 0.01), RealVect::Constant(2, 0.1),
                         CplxVect::Zero(2), from_id, to_id);
    RealVect load_p(1), load_q(1);
    load_p << 50.;
    load_q << 10.;
    Eigen::VectorXi load_bus(1);
    load_bus << 2;
    grid.init_loads(load_p, load_q, load_bus);
    const int n_gen = pv_gen ? 2 : 1;
    RealVect gen_p(n_gen), gen_v(n_gen), gen_min_q(n_gen), gen_max_q(n_gen);
    Eigen::VectorXi gen_bus(n_gen);
    gen_p(0) = 0.;
    gen_v(0) = 1.02;
    gen_min_q(0) = -1000.;
    gen_max_q(0) = 1000.;
    gen_bus(0) = 0;
    if (pv_gen) {
        gen_p(1) = 10.;
        gen_v(1) = 1.01;
        gen_min_q(1) = -100.;
        gen_max_q(1) = 100.;
        gen_bus(1) = 1;
    }
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    return grid;
}

// a loop that answers its check() calls from a script (its last entry repeats), logs
// "<name>.<call>" for every call, and on UNSTABLE nudges every injection so that the next
// Newton solve has something to do
class ScriptedLoop final : public BaseOuterLoop
{
    public:
        ScriptedLoop(std::string name, std::vector<OuterLoopStatus> script,
                     std::shared_ptr<std::vector<std::string> > log,
                     bool needed = true, bool fixes_unrealistic = false):
            name_(std::move(name)), script_(std::move(script)), log_(std::move(log)),
            needed_(needed), fixes_unrealistic_(fixes_unrealistic) {}

    protected:
        std::string _name() const override { return name_; }
        bool _is_needed(const OuterContext &) const override { return needed_; }
        void _initialize(OuterContext &) override {
            calls_ = 0;
            log_->push_back(name_ + ".initialize");
        }
        void _detect(const OuterContext &, std::vector<LimitViolation> &) const override {}
        OuterLoopStatus _check(OuterContext & ctx) override {
            log_->push_back(name_ + ".check" + std::to_string(ctx.iteration));
            const std::size_t pos = std::min(calls_, script_.size() - 1);
            ++calls_;
            const OuterLoopStatus res = script_[pos];
            if (res == OuterLoopStatus::UNSTABLE) *ctx.state->Sbus *= 1.01;
            return res;
        }
        void _cleanup(OuterContext &) override { log_->push_back(name_ + ".cleanup"); }
        bool _can_fix_unrealistic_state() const override { return fixes_unrealistic_; }
        std::unique_ptr<BaseOuterLoop> _clone() const override {
            return std::unique_ptr<BaseOuterLoop>(new ScriptedLoop(*this));
        }

    private:
        std::string name_;
        std::vector<OuterLoopStatus> script_;
        std::shared_ptr<std::vector<std::string> > log_;
        bool needed_;
        bool fixes_unrealistic_;
        std::size_t calls_ = 0;
};

// a loop that only reserves a switchable bus (a Vm unknown + a Q equation) for every PV bus
class DeclaringLoop final : public BaseOuterLoop
{
    protected:
        std::string _name() const override { return "Declaring"; }
        void _declare(const OuterContext & ctx, ls2g::OuterDeclaration &) const override {
            const ls2g::SolverBusIdVect & pv = ctx.grid->get_ac_pv_solver();
            for (std::size_t i = 0; i < pv.size(); ++i) ctx.controls->reserve_bus_voltage(pv[i].cast_int());
        }
        void _detect(const OuterContext &, std::vector<LimitViolation> &) const override {}
        OuterLoopStatus _check(OuterContext &) override { return OuterLoopStatus::STABLE; }
        std::unique_ptr<BaseOuterLoop> _clone() const override {
            return std::unique_ptr<BaseOuterLoop>(new DeclaringLoop(*this));
        }
};

const OuterLoopStatus S = OuterLoopStatus::STABLE;
const OuterLoopStatus U = OuterLoopStatus::UNSTABLE;
const OuterLoopStatus F = OuterLoopStatus::FAILED;

CplxVect solve(LSGrid & grid)
{
    return grid.ac_pf(CplxVect::Constant(grid.total_bus(), ls2g::cplx_type(1., 0.)), 30, 1e-8);
}

}  // namespace

TEST_CASE("an outer-loop solve without any loop is the single-slack Newton's", "[outer_loop]")
{
    LSGrid ref = make_grid();
    ref.change_algorithm("NRSing_SparseLU");
    const CplxVect V_ref = solve(ref);
    REQUIRE(V_ref.size() == 3);

    LSGrid grid = make_grid();
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    const CplxVect V = solve(grid);
    REQUIRE(V.size() == 3);
    for (int i = 0; i < 3; ++i) REQUIRE(V(i) == V_ref(i));

    const OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::STABLE);
    CHECK(stats.nb_outer_iterations == 0);
    CHECK(stats.nr_iterations.size() == 1);
}

TEST_CASE("the loops run in passes, as OpenLoadFlow's engine runs them", "[outer_loop]")
{
    auto log = std::make_shared<std::vector<std::string> >();
    LSGrid grid = make_grid();
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    // A changes something twice, B once
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("A", std::vector<OuterLoopStatus>{U, U, S}, log));
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("B", std::vector<OuterLoopStatus>{U, S}, log));
    REQUIRE(solve(grid).size() == 3);

    // pass 1: A until stable, then B until stable (B is now the last unstable loop);
    // pass 2: A is checked again (stable), and the pass stops before B -- every other loop
    // has been checked since B's change -- and adds no iteration, so there is no pass 3
    const std::vector<std::string> expected{
        "A.initialize", "B.initialize",
        "A.check0", "A.check1", "A.check2", "B.check0", "B.check1",
        "A.check2",
        "B.cleanup", "A.cleanup"};
    CHECK(*log == expected);

    const OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::STABLE);
    CHECK(stats.nb_outer_iterations == 3);
    CHECK(stats.nb_passes == 2);
    REQUIRE(stats.loop_iterations.size() == 2);
    CHECK(stats.loop_iterations[0].second == 2);
    CHECK(stats.loop_iterations[1].second == 1);
    CHECK(stats.nr_iterations.size() == 4);  // the first solve and one per change
}

TEST_CASE("a whole outer-loop solve is one symbolic analysis", "[outer_loop]")
{
    auto log = std::make_shared<std::vector<std::string> >();
    LSGrid grid = make_grid();
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("A", std::vector<OuterLoopStatus>{U, U, U, S}, log));
    REQUIRE(solve(grid).size() == 3);
    REQUIRE(solve(grid).size() == 3);  // and a second solve on the same grid reuses it
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
}

TEST_CASE("the outer iterations are capped", "[outer_loop]")
{
    auto log = std::make_shared<std::vector<std::string> >();
    LSGrid grid = make_grid();
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("A", std::vector<OuterLoopStatus>{U}, log));
    CHECK(solve(grid).size() == 0);
    const OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::UNSTABLE);
    CHECK(stats.nb_outer_iterations == 30);
    CHECK(grid.get_algo().get_error() == ErrorType::TooManyIterations);
}

TEST_CASE("a failed loop stops the solve", "[outer_loop]")
{
    auto log = std::make_shared<std::vector<std::string> >();
    LSGrid grid = make_grid();
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("A", std::vector<OuterLoopStatus>{F}, log));
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("B", std::vector<OuterLoopStatus>{U, S}, log));
    CHECK(solve(grid).size() == 0);
    const OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::FAILED);
    CHECK(stats.failed_loop == "A");
    CHECK(grid.get_algo().get_error() == ErrorType::OuterLoopFailed);
    const std::vector<std::string> expected{
        "A.initialize", "B.initialize", "A.check0", "B.cleanup", "A.cleanup"};
    CHECK(*log == expected);
}

TEST_CASE("a loop with nothing to do is left out", "[outer_loop]")
{
    auto log = std::make_shared<std::vector<std::string> >();
    LSGrid grid = make_grid();
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("A", std::vector<OuterLoopStatus>{U, S}, log, false));
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("B", std::vector<OuterLoopStatus>{S}, log));
    REQUIRE(solve(grid).size() == 3);
    const std::vector<std::string> expected{"B.initialize", "B.check0", "B.cleanup"};
    CHECK(*log == expected);
}

TEST_CASE("the unrealistic-voltage check", "[outer_loop]")
{
    SECTION("every solve is checked when no loop can fix it") {
        LSGrid grid = make_grid();
        grid.change_algorithm("NROuter_SparseLU");
        grid.clear_outer_loops();
        // every bus is "unrealistic" in a band that excludes 1 pu, and is checked whatever
        // its nominal voltage
        ls2g::AlgoConfig cfg = grid.get_ac_algo_config();
        REQUIRE(cfg.real_params.size() >= 9);
        cfg.real_params[6] = 1.05;
        cfg.real_params[7] = 1.2;
        cfg.real_params[8] = 0.;
        grid.set_ac_algo_config(cfg);
        CHECK(solve(grid).size() == 0);
        CHECK(grid.get_algo().get_error() == ErrorType::UnrealisticState);
        CHECK(grid.get_algo().get_outer_loop_stats().unrealistic_state);
    }
    SECTION("only buses at or above the nominal voltage threshold are checked") {
        LSGrid grid = make_grid();
        grid.change_algorithm("NROuter_SparseLU");
        grid.clear_outer_loops();
        ls2g::AlgoConfig cfg = grid.get_ac_algo_config();
        cfg.real_params[6] = 1.05;  // the grid is 138 kV, below the default 180 kV threshold
        grid.set_ac_algo_config(cfg);
        CHECK(solve(grid).size() == 3);
    }
    SECTION("in robust mode, a loop able to fix it postpones the check") {
        auto log = std::make_shared<std::vector<std::string> >();
        LSGrid grid = make_grid();
        grid.change_algorithm("NROuter_SparseLU");
        grid.clear_outer_loops();
        grid.add_outer_loop(std::make_shared<ScriptedLoop>("A", std::vector<OuterLoopStatus>{U, S}, log, true, true));
        ls2g::AlgoConfig cfg = grid.get_ac_algo_config();
        cfg.real_params[6] = 1.05;
        cfg.real_params[8] = 0.;
        grid.set_ac_algo_config(cfg);
        // the loop still ran (the first solve was not rejected), then the check failed
        CHECK(solve(grid).size() == 0);
        CHECK(grid.get_algo().get_error() == ErrorType::UnrealisticState);
        const std::vector<std::string> expected{"A.initialize", "A.check0", "A.check1", "A.cleanup"};
        CHECK(*log == expected);
    }
}

TEST_CASE("a copied grid runs its own loops", "[outer_loop]")
{
    auto log = std::make_shared<std::vector<std::string> >();
    LSGrid grid = make_grid();
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("A", std::vector<OuterLoopStatus>{S}, log));
    LSGrid copy(grid);
    REQUIRE(copy.get_outer_loops().size() == 1);
    CHECK(copy.get_outer_loops()[0] != grid.get_outer_loops()[0]);
    CHECK(copy.get_algo().get_name() == "NROuter_SparseLU");
    REQUIRE(solve(copy).size() == 3);
}

TEST_CASE("the same loop cannot be listed twice", "[outer_loop]")
{
    auto log = std::make_shared<std::vector<std::string> >();
    LSGrid grid = make_grid();
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("A", std::vector<OuterLoopStatus>{S}, log));
    CHECK_THROWS(grid.add_outer_loop(std::make_shared<ScriptedLoop>("A", std::vector<OuterLoopStatus>{S}, log)));
}

TEST_CASE("an algorithm without outer loops ignores the list", "[outer_loop]")
{
    auto log = std::make_shared<std::vector<std::string> >();
    LSGrid grid = make_grid();
    grid.change_algorithm("NRSing_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ScriptedLoop>("A", std::vector<OuterLoopStatus>{U}, log));
    REQUIRE(solve(grid).size() == 3);
    CHECK(log->empty());
}

TEST_CASE("what a loop declares is reserved in the Jacobian", "[outer_loop]")
{
    LSGrid grid = make_grid(true);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    REQUIRE(solve(grid).size() == 3);
    const int pv_bus = grid.id_me_to_ac_solver()[1].cast_int();
    const Eigen::Index n_without = grid.get_algo().get_J().rows();
    CHECK(grid.get_algo().get_vm_to_J_col_python()(pv_bus) == -1);

    grid.add_outer_loop(std::make_shared<DeclaringLoop>());
    REQUIRE(solve(grid).size() == 3);
    // the PV bus got a Vm column and a Q row: one more of each
    CHECK(grid.get_algo().get_J().rows() == n_without + 1);
    CHECK(grid.get_algo().get_vm_to_J_col_python()(pv_bus) >= 0);
}

TEST_CASE("a rejected configuration changes nothing", "[outer_loop]")
{
    LSGrid grid = make_grid();
    grid.change_algorithm("NROuter_SparseLU");
    const ls2g::AlgoConfig before = grid.get_ac_algo_config();
    ls2g::AlgoConfig cfg = before;
    cfg.int_params[4] = 12;   // a valid driver part ...
    cfg.int_params[3] = 0;    // ... with an invalid Newton part (refactor_every_n)
    CHECK_THROWS(grid.set_ac_algo_config(cfg));
    CHECK(grid.get_ac_algo_config().int_params == before.int_params);
    CHECK(grid.get_ac_algo_config().real_params == before.real_params);
}
