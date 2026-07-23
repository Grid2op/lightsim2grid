// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Registration-safety tests for the solver-plugin mechanism. The concern these
// guard is specifically the plugin boundary: a plugin registers solvers from
// code that lightsim2grid loads at runtime, and the two failure modes that
// mechanism has -- an ABI-tag mismatch and a duplicate solver name -- must
// (a) never leave the registry half-updated, and (b) never let a C++ exception
// escape the exported entry point (where, running under dlopen's C frames, it
// would be an uncatchable std::terminate). The errbuf-handling tests also run
// under valgrind in the cpp_unit_tests workflow, which is what actually proves
// the bounded copy never writes out of bounds.

#include <catch2/catch_test_macros.hpp>

#include <cstring>
#include <memory>
#include <string>

#include "AlgorithmRegistry.hpp"
#include "Ls2gAbiTag.hpp"

namespace {

ls2g::AlgorithmRegistry::Factory null_factory() {
    // Registration must never need to *build* a solver, so the factory can be
    // a no-op: these tests only exercise the registration bookkeeping.
    return []() -> std::unique_ptr<ls2g::BaseAlgo> { return nullptr; };
}

// Referenced from the captureless registration lambda in the overflow test
// below (a function pointer cannot capture, but may read a static).
std::string g_long_solver_name;

} // namespace

TEST_CASE("register_all commits every name on success", "[plugin]")
{
    ls2g::AlgorithmRegistry::FactoryMap batch;
    batch.emplace("__plugin_ok_A__", null_factory());
    batch.emplace("__plugin_ok_B__", null_factory());

    ls2g::AlgorithmRegistry::instance().register_all(std::move(batch));

    REQUIRE(ls2g::AlgorithmRegistry::instance().is_registered("__plugin_ok_A__"));
    REQUIRE(ls2g::AlgorithmRegistry::instance().is_registered("__plugin_ok_B__"));
}

TEST_CASE("register_all is atomic: a colliding name rolls back the whole batch", "[plugin]")
{
    // Pre-register a name the batch will collide with.
    ls2g::AlgorithmRegistry::instance().register_solver("__plugin_existing__", null_factory());

    ls2g::AlgorithmRegistry::FactoryMap batch;
    batch.emplace("__plugin_fresh_should_not_survive__", null_factory());
    batch.emplace("__plugin_existing__", null_factory());  // collision

    REQUIRE_THROWS_AS(
        ls2g::AlgorithmRegistry::instance().register_all(std::move(batch)),
        std::invalid_argument);

    // The fresh name from the same (rejected) batch must NOT have been added.
    REQUIRE_FALSE(ls2g::AlgorithmRegistry::instance().is_registered("__plugin_fresh_should_not_survive__"));
    // ... and the pre-existing one is untouched.
    REQUIRE(ls2g::AlgorithmRegistry::instance().is_registered("__plugin_existing__"));
}

TEST_CASE("register_all rejects a mismatched ABI tag and commits nothing", "[plugin]")
{
    ls2g::Ls2gAbiTag mismatched = ls2g::ls2g_current_abi_tag();
    mismatched.eigen_max_align_bytes += 16;  // guaranteed to differ from core's tag

    ls2g::AlgorithmRegistry::FactoryMap batch;
    batch.emplace("__plugin_abi_reject__", null_factory());

    REQUIRE_THROWS_AS(
        ls2g::AlgorithmRegistry::instance().register_all(std::move(batch), mismatched),
        std::runtime_error);

    REQUIRE_FALSE(ls2g::AlgorithmRegistry::instance().is_registered("__plugin_abi_reject__"));
}

TEST_CASE("PluginRegistrar rejects a name the same plugin stages twice", "[plugin]")
{
    ls2g::PluginRegistrar reg;
    reg.add("__plugin_intra_dup__", null_factory());
    REQUIRE_THROWS_AS(reg.add("__plugin_intra_dup__", null_factory()), std::invalid_argument);
}

TEST_CASE("plugin entry helper: success returns 0 and registers", "[plugin]")
{
    char errbuf[256] = {'x'};  // deliberately not empty to confirm it stays untouched-but-terminated
    int rc = ls2g::detail::ls2g_run_plugin_entry(
        &ls2g::AlgorithmRegistry::instance(), errbuf, sizeof(errbuf),
        +[](ls2g::PluginRegistrar& r) { r.add("__plugin_entry_ok__", null_factory()); },
        ls2g::ls2g_current_abi_tag());

    REQUIRE(rc == 0);
    REQUIRE(ls2g::AlgorithmRegistry::instance().is_registered("__plugin_entry_ok__"));
}

TEST_CASE("plugin entry helper: a registration failure is a return code, never an exception", "[plugin]")
{
    ls2g::AlgorithmRegistry::instance().register_solver("__plugin_entry_dup__", null_factory());

    char errbuf[256] = {'\0'};
    int rc = ls2g::detail::ls2g_run_plugin_entry(
        &ls2g::AlgorithmRegistry::instance(), errbuf, sizeof(errbuf),
        +[](ls2g::PluginRegistrar& r) { r.add("__plugin_entry_dup__", null_factory()); },
        ls2g::ls2g_current_abi_tag());

    REQUIRE(rc != 0);
    REQUIRE(std::strlen(errbuf) > 0);  // a diagnostic was written
}

TEST_CASE("plugin entry helper: a short error buffer is never overflowed", "[plugin]")
{
    // Force a (long) error message, then hand the helper a tiny buffer.
    ls2g::AlgorithmRegistry::instance().register_solver("__plugin_short_buf_dup__", null_factory());

    constexpr std::size_t N = 8;
    char errbuf[N];
    std::memset(errbuf, 'Z', sizeof(errbuf));  // no NUL anywhere to start with

    int rc = ls2g::detail::ls2g_run_plugin_entry(
        &ls2g::AlgorithmRegistry::instance(), errbuf, N,
        +[](ls2g::PluginRegistrar& r) { r.add("__plugin_short_buf_dup__", null_factory()); },
        ls2g::ls2g_current_abi_tag());

    REQUIRE(rc != 0);
    REQUIRE(errbuf[N - 1] == '\0');          // always NUL-terminated within the bound
    REQUIRE(std::strlen(errbuf) <= N - 1);   // never wrote past the buffer
}

TEST_CASE("plugin entry helper: a solver name longer than the buffer cannot overflow it", "[plugin]")
{
    // The scenario from the review: a plugin registers a name so long that the
    // resulting error message far exceeds the host's fixed 512-byte diagnostic
    // buffer. The bounded write must truncate it, never overflow. Under valgrind
    // (cpp_unit_tests workflow) this also proves there is no out-of-bounds write.
    g_long_solver_name = std::string(1000, 'x');  // >> 512
    ls2g::AlgorithmRegistry::instance().register_solver(g_long_solver_name, null_factory());

    constexpr std::size_t host_buf_size = 512;  // matches load_algorithm_plugin_impl
    char errbuf[host_buf_size];
    std::memset(errbuf, 'Z', sizeof(errbuf));  // start with no NUL anywhere

    int rc = ls2g::detail::ls2g_run_plugin_entry(
        &ls2g::AlgorithmRegistry::instance(), errbuf, host_buf_size,
        // staging the (already-registered) long name triggers the long
        // "... is already registered" message, which contains the 1000-char name
        +[](ls2g::PluginRegistrar& r) { r.add(g_long_solver_name, null_factory()); },
        ls2g::ls2g_current_abi_tag());

    REQUIRE(rc != 0);
    REQUIRE(errbuf[host_buf_size - 1] == '\0');          // terminated within bounds
    REQUIRE(std::strlen(errbuf) <= host_buf_size - 1);   // never wrote past the buffer
}

TEST_CASE("plugin entry helper: a null / zero-length buffer is handled safely", "[plugin]")
{
    ls2g::AlgorithmRegistry::instance().register_solver("__plugin_nullbuf_dup__", null_factory());

    // null buffer on an error path -> no write, no crash.
    int rc1 = ls2g::detail::ls2g_run_plugin_entry(
        &ls2g::AlgorithmRegistry::instance(), nullptr, 0,
        +[](ls2g::PluginRegistrar& r) { r.add("__plugin_nullbuf_dup__", null_factory()); },
        ls2g::ls2g_current_abi_tag());
    REQUIRE(rc1 != 0);

    // zero-length (non-null) buffer: the single byte must not be written.
    char one = 'Q';
    int rc2 = ls2g::detail::ls2g_run_plugin_entry(
        &ls2g::AlgorithmRegistry::instance(), &one, 0,
        +[](ls2g::PluginRegistrar& r) { r.add("__plugin_nullbuf_dup__", null_factory()); },
        ls2g::ls2g_current_abi_tag());
    REQUIRE(rc2 != 0);
    REQUIRE(one == 'Q');
}

TEST_CASE("plugin entry helper: a null registry pointer is reported, not dereferenced", "[plugin]")
{
    char errbuf[128] = {'\0'};
    int rc = ls2g::detail::ls2g_run_plugin_entry(
        nullptr, errbuf, sizeof(errbuf),
        +[](ls2g::PluginRegistrar& r) { r.add("__plugin_never__", null_factory()); },
        ls2g::ls2g_current_abi_tag());
    REQUIRE(rc == 2);
    REQUIRE(std::strlen(errbuf) > 0);
}
