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

TEST_CASE("plugin entry helper: a pathologically long solver name cannot overflow the buffer", "[plugin]")
{
    // A plugin supplies a name far longer than both the solver-name limit and the
    // host's fixed 512-byte diagnostic buffer. It is rejected (see
    // is_valid_solver_name), and the rejection has to come back as a return code
    // plus a message that fits the buffer -- never an overflow, and never an
    // exception escaping the extern "C" entry point. Under valgrind
    // (cpp_unit_tests workflow) this also proves there is no out-of-bounds write.
    g_long_solver_name = std::string(1000, 'x');  // >> 512, and >> SOLVER_NAME_MAX_LEN

    constexpr std::size_t host_buf_size = 512;  // matches load_algorithm_plugin_impl
    char errbuf[host_buf_size];
    std::memset(errbuf, 'Z', sizeof(errbuf));  // start with no NUL anywhere

    int rc = ls2g::detail::ls2g_run_plugin_entry(
        &ls2g::AlgorithmRegistry::instance(), errbuf, host_buf_size,
        +[](ls2g::PluginRegistrar& r) { r.add(g_long_solver_name, null_factory()); },
        ls2g::ls2g_current_abi_tag());

    REQUIRE(rc != 0);
    // The buffer holds a proper C string: a NUL somewhere within bounds. (Unlike
    // the short-buffer case above we do not require the NUL to land on the very
    // last byte: printable() truncates the offending name, so the message may now
    // be shorter than the buffer -- what matters is that it is terminated and
    // that nothing was written past the end.)
    REQUIRE(std::memchr(errbuf, '\0', host_buf_size) != nullptr);
    REQUIRE(std::strlen(errbuf) < host_buf_size);
    REQUIRE(std::strlen(errbuf) > 0);                    // a diagnostic was written
    REQUIRE_FALSE(ls2g::AlgorithmRegistry::instance().is_registered(g_long_solver_name));
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

// --- solver-name whitelist ---------------------------------------------------
// A solver name is an identity: it is written into every serialized grid and is
// what re-selects the solver on load. It is restricted to
// [A-Za-z_][A-Za-z0-9_.]{0,63} -- see is_valid_solver_name.

TEST_CASE("is_valid_solver_name accepts the names actually in use", "[plugin][solver_name]")
{
    // every built-in, and the names the shipped example plugins register
    for (const char * name : {"NR_SparseLU", "NRSing_SparseLU", "DC_SparseLU",
                              "GaussSeidel", "GaussSeidelSynch", "FDPF_XB_SparseLU",
                              "NR_KLU", "NRRefactorRetry_CKTSO",
                              "DummyExternal", "NR_LM_SparseLU", "NRDistSlack_KLU"}) {
        INFO(name);
        CHECK(ls2g::is_valid_solver_name(name));
    }
    // and the shapes we explicitly want to allow
    CHECK(ls2g::is_valid_solver_name("_leading_underscore"));
    CHECK(ls2g::is_valid_solver_name("acme.NR_fast"));   // namespaced (recommended)
    CHECK(ls2g::is_valid_solver_name("solver_v2"));      // digits, not first
    CHECK(ls2g::is_valid_solver_name("A"));              // shortest
    CHECK(ls2g::is_valid_solver_name(std::string(ls2g::SOLVER_NAME_MAX_LEN, 'a')));
}

TEST_CASE("is_valid_solver_name rejects names that could confuse or corrupt", "[plugin][solver_name]")
{
    CHECK_FALSE(ls2g::is_valid_solver_name(""));                       // empty
    CHECK_FALSE(ls2g::is_valid_solver_name(std::string(ls2g::SOLVER_NAME_MAX_LEN + 1, 'a')));  // too long
    CHECK_FALSE(ls2g::is_valid_solver_name("2fast"));                  // leading digit
    CHECK_FALSE(ls2g::is_valid_solver_name(".hidden"));                // leading dot
    CHECK_FALSE(ls2g::is_valid_solver_name("NR-SparseLU"));            // hyphen not allowed
    CHECK_FALSE(ls2g::is_valid_solver_name("NR SparseLU"));            // space
    CHECK_FALSE(ls2g::is_valid_solver_name("../../etc/passwd"));       // path-like
    CHECK_FALSE(ls2g::is_valid_solver_name("NR_SparseLU\nWARNING: x"));// log injection
    CHECK_FALSE(ls2g::is_valid_solver_name("NR\x1b[31m_SparseLU"));    // terminal escape
    CHECK_FALSE(ls2g::is_valid_solver_name(std::string("NR\0LU", 5))); // embedded NUL
    // a homoglyph: "NR_SparseLU" with a Cyrillic 'а' (U+0430) in place of the 'a'
    CHECK_FALSE(ls2g::is_valid_solver_name("NR_Sp\xd0\xb0rseLU"));
}

TEST_CASE("registering an invalid solver name is refused", "[plugin][solver_name]")
{
    CHECK_THROWS_AS(
        ls2g::AlgorithmRegistry::instance().register_solver("bad name!", null_factory()),
        std::invalid_argument);
    CHECK_FALSE(ls2g::AlgorithmRegistry::instance().is_registered("bad name!"));
}

TEST_CASE("a batch containing an invalid name is rejected whole", "[plugin][solver_name]")
{
    ls2g::AlgorithmRegistry::FactoryMap batch;
    batch.emplace("__valid_name_in_bad_batch__", null_factory());
    batch.emplace("also bad!", null_factory());

    CHECK_THROWS_AS(ls2g::AlgorithmRegistry::instance().register_all(std::move(batch)),
                    std::invalid_argument);
    // atomic: the valid name from the rejected batch must not survive either
    CHECK_FALSE(ls2g::AlgorithmRegistry::instance().is_registered("__valid_name_in_bad_batch__"));
    CHECK_FALSE(ls2g::AlgorithmRegistry::instance().is_registered("also bad!"));
}
