// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Deterministic corruption sweep: the C++ twin of the python sweep in
// lightsim2grid/tests/test_binary_serialization.py, but running on the
// archive layer alone -- so it is cheap enough to run at every byte offset,
// and any memory error it provokes becomes a hard failure under the valgrind
// / ASan CI jobs. The contract under test: loading arbitrarily corrupted
// bytes either succeeds (the corruption kept all counts consistent) or
// throws std::runtime_error -- never anything else, never a crash or an
// attempted huge allocation.

#include <cstring>   // std::memset
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "BinaryArchive.hpp"
#include "test_helpers.hpp"

using ls2g_test::FakeContainer;

namespace {

// writes the corrupted bytes and attempts a load; anything other than a
// clean success or a std::runtime_error fails the test at this offset
void check_load(const std::string & path, const std::vector<char> & bytes, std::size_t offset)
{
    ls2g_test::write_file(path, bytes);
    INFO("corruption at byte offset " << offset << " of " << bytes.size());
    try {
        (void) ls2g::load_binary_generic<FakeContainer>(path, "0", "14", "0");
    } catch (const std::runtime_error &) {
        // the expected rejection path
    }
    // any other exception type (std::bad_alloc, std::length_error, ...)
    // propagates out of the loop and fails the surrounding TEST_CASE
    SUCCEED();
}

std::vector<char> reference_bytes(const std::string & path)
{
    ls2g::save_binary_generic(ls2g_test::make_reference_container(), path, "0", "14", "0");
    return ls2g_test::read_file(path);
}

}  // anonymous namespace

TEST_CASE("every single-byte flip loads cleanly or throws runtime_error", "[corruption]")
{
    ls2g_test::TempFile file;
    const std::vector<char> ref = reference_bytes(file.str());
    REQUIRE(ref.size() > 0);
    for (std::size_t i = 0; i < ref.size(); ++i) {
        std::vector<char> corrupted = ref;
        corrupted[i] = static_cast<char>(corrupted[i] ^ 0xFF);
        check_load(file.str(), corrupted, i);
    }
}

TEST_CASE("every 8-byte 0xFF stamp loads cleanly or throws runtime_error", "[corruption]")
{
    // stamps a maximal uint64 over every alignment, the worst case for any
    // count / length field it happens to cover
    ls2g_test::TempFile file;
    const std::vector<char> ref = reference_bytes(file.str());
    REQUIRE(ref.size() >= 8);
    for (std::size_t i = 0; i + 8 <= ref.size(); ++i) {
        std::vector<char> corrupted = ref;
        std::memset(&corrupted[i], 0xFF, 8);
        check_load(file.str(), corrupted, i);
    }
}

TEST_CASE("every truncation loads cleanly or throws runtime_error", "[corruption]")
{
    ls2g_test::TempFile file;
    const std::vector<char> ref = reference_bytes(file.str());
    for (std::size_t len = 0; len < ref.size(); ++len) {
        const std::vector<char> truncated(ref.begin(), ref.begin() + len);
        check_load(file.str(), truncated, len);
    }
}
