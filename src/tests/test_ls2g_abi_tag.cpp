// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Ls2gAbiTag is the mechanism that catches, at plugin-registration / bindings
// import time, the exact bug this repository hit once already: core and a
// separately-compiled consumer (a solver plugin, or lightsim2grid_cpp itself)
// disagreeing on -march=native (and therefore on Eigen's alignment), which
// otherwise silently corrupts the heap the first time an Eigen object crosses
// the boundary. These tests can't reproduce a real cross-flag build inside a
// single binary, so they exercise the comparison/rejection logic directly by
// constructing a deliberately-mismatched tag.

#include <catch2/catch_test_macros.hpp>

#include "AlgorithmRegistry.hpp"
#include "Ls2gAbiTag.hpp"

TEST_CASE("Ls2gAbiTag equality reflects every field", "[ls2g_abi_tag]")
{
    ls2g::Ls2gAbiTag a = ls2g::ls2g_current_abi_tag();
    ls2g::Ls2gAbiTag b = a;
    REQUIRE(a == b);

    b = a;
    b.eigen_max_align_bytes += 16;
    REQUIRE(a != b);

    b = a;
    b.vectorize_avx2 = !b.vectorize_avx2;
    REQUIRE(a != b);

    b = a;
    b.tag_version += 1;
    REQUIRE(a != b);

    // Eigen version is compared in full: two independently-compiled TUs
    // built against different Eigen copies (even just a patch bump apart)
    // are treated as a mismatch -- see the comment on operator== for why.
    b = a;
    b.eigen_patch_version += 1;
    REQUIRE(a != b);
}

TEST_CASE("register_solver rejects a solver registered with a mismatched ABI tag", "[ls2g_abi_tag]")
{
    ls2g::Ls2gAbiTag mismatched = ls2g::ls2g_current_abi_tag();
    mismatched.eigen_max_align_bytes += 16;  // guaranteed to differ from core's own tag

    // The factory is never invoked -- the check must reject the registration
    // before any solver object (and any Eigen member it holds) is ever built.
    REQUIRE_THROWS_AS(
        ls2g::AlgorithmRegistry::instance().register_solver(
            "__test_mismatched_abi_tag_plugin__",
            []() -> std::unique_ptr<ls2g::BaseAlgo> { return nullptr; },
            mismatched),
        std::runtime_error);

    REQUIRE_FALSE(ls2g::AlgorithmRegistry::instance().is_registered("__test_mismatched_abi_tag_plugin__"));
}

TEST_CASE("register_solver accepts a solver registered with a matching ABI tag", "[ls2g_abi_tag]")
{
    bool called = false;
    ls2g::AlgorithmRegistry::instance().register_solver(
        "__test_matching_abi_tag_plugin__",
        [&called]() -> std::unique_ptr<ls2g::BaseAlgo> { called = true; return std::make_unique<ls2g::BaseAlgo>(); });
        // caller_tag defaults to ls2g_current_abi_tag(), i.e. this TU's own --
        // which is core's own TU here, so it always matches core_abi_tag().

    REQUIRE(ls2g::AlgorithmRegistry::instance().is_registered("__test_matching_abi_tag_plugin__"));
    auto algo = ls2g::AlgorithmRegistry::instance().make("__test_matching_abi_tag_plugin__");
    REQUIRE(called);
    REQUIRE(algo != nullptr);
}
