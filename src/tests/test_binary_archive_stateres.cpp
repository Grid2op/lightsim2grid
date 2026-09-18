// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// ValueArchiver dispatch + the generic save_binary_generic / load_binary_generic
// entry points, over synthetic StateRes types (test_helpers.hpp) covering every
// field shape the real containers use -- no grid, no Eigen data, no python.

#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "BinaryArchive.hpp"
#include "test_helpers.hpp"

using ls2g_test::FakeContainer;
using ls2g_test::TempFile;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;

namespace {

// number of bytes the header occupies for a given tag, with the "0"/"14"/"0"
// version strings every test here writes: magic (4) + format (uint32) + four
// length-prefixed strings (uint32 length each)
std::size_t header_size(const std::string & tag)
{
    return 4 + 4 + (4 + tag.size()) + (4 + 1) + (4 + 2) + (4 + 1);
}

template<typename T>
void save(const T & obj, const std::string & path, bool atomic = true)
{
    ls2g::save_binary_generic(obj, path, "0", "14", "0", atomic);
}

template<typename T>
T load(const std::string & path)
{
    return ls2g::load_binary_generic<T>(path, "0", "14", "0");
}

}  // anonymous namespace

TEST_CASE("full StateRes round-trip covering every ValueArchiver specialization", "[StateRes]")
{
    const FakeContainer ref = ls2g_test::make_reference_container();
    TempFile file;
    save(ref, file.str());
    const FakeContainer loaded = load<FakeContainer>(file.str());
    // std::tuple operator== compares field by field, recursing into the
    // nested tuple and the nested vectors
    CHECK(loaded.get_state() == ref.get_state());
}

TEST_CASE("an all-default (empty) StateRes round-trips", "[StateRes]")
{
    const FakeContainer ref{};  // zero scalars, empty strings and vectors
    TempFile file;
    save(ref, file.str());
    CHECK(load<FakeContainer>(file.str()).get_state() == ref.get_state());
}

TEST_CASE("non-atomic save round-trips too", "[StateRes]")
{
    const FakeContainer ref = ls2g_test::make_reference_container();
    TempFile file;
    save(ref, file.str(), false);
    CHECK_FALSE(ls2g_test::file_exists(file.str() + ".lsb_tmp"));
    CHECK(load<FakeContainer>(file.str()).get_state() == ref.get_state());
}

TEST_CASE("a corrupted on-disk bool byte is normalized, not undefined behavior", "[StateRes]")
{
    // reading a byte holding eg 255 straight into a bool object would be UB
    // when that bool is later loaded; ValueArchiver<bool> normalizes instead
    const FakeContainer ref = ls2g_test::make_reference_container();
    TempFile file;
    save(ref, file.str());

    // the bool is the third StateRes field: after the int and the real_type
    std::vector<char> bytes = ls2g_test::read_file(file.str());
    const std::size_t bool_offset =
        header_size(FakeContainer::binary_type_tag()) + sizeof(int) + sizeof(ls2g::real_type);
    REQUIRE(bytes.at(bool_offset) == 1);  // written by the reference container as true
    bytes.at(bool_offset) = static_cast<char>(0xFF);
    ls2g_test::write_file(file.str(), bytes);

    const FakeContainer loaded = load<FakeContainer>(file.str());
    CHECK(std::get<2>(loaded.get_state()) == true);

    bytes.at(bool_offset) = 0;
    ls2g_test::write_file(file.str(), bytes);
    CHECK(std::get<2>(load<FakeContainer>(file.str()).get_state()) == false);
}

TEST_CASE("identical StateRes layouts are told apart by the type tag", "[StateRes]")
{
    // FakeLoad and FakeStorage have byte-identical payloads (like the real
    // LoadContainer / StorageContainer): only the header tag differs
    ls2g_test::FakeLoad load_obj;
    load_obj.state = std::make_tuple(3, std::vector<ls2g::real_type>{1., 2.});
    TempFile file;
    save(load_obj, file.str());

    CHECK(load<ls2g_test::FakeLoad>(file.str()).get_state() == load_obj.get_state());
    REQUIRE_THROWS_MATCHES(
        load<ls2g_test::FakeStorage>(file.str()), std::runtime_error,
        MessageMatches(ContainsSubstring("wrong object type") &&
                       ContainsSubstring("'FakeLoad'") &&
                       ContainsSubstring("'FakeStorage'")));
}

TEST_CASE("fixup_binary_state is applied when declared, and only then", "[StateRes]")
{
    SECTION("a type declaring the hook gets it applied after the read") {
        ls2g_test::FakeWithFixup obj;
        obj.state = std::make_tuple(5, std::string("original"));
        TempFile file;
        save(obj, file.str());
        const ls2g_test::FakeWithFixup loaded = load<ls2g_test::FakeWithFixup>(file.str());
        CHECK(std::get<0>(loaded.get_state()) == 5);
        CHECK(std::get<1>(loaded.get_state()) == "fixed up");
    }
    SECTION("a type without the hook is loaded untouched") {
        ls2g_test::FakeLoad obj;
        obj.state = std::make_tuple(5, std::vector<ls2g::real_type>{7.});
        TempFile file;
        save(obj, file.str());
        CHECK(load<ls2g_test::FakeLoad>(file.str()).get_state() == obj.get_state());
    }
}

TEST_CASE("loading rejects a file whose payload is longer than the StateRes", "[StateRes]")
{
    // a FakeContainer payload after a FakeLoad-sized read leaves trailing
    // bytes: check_fully_consumed must reject them (wrong-class protection
    // beyond the tag, eg after a manual rename attack on the header)
    ls2g_test::FakeLoad obj;
    obj.state = std::make_tuple(1, std::vector<ls2g::real_type>{1.});
    TempFile file;
    save(obj, file.str());
    std::vector<char> bytes = ls2g_test::read_file(file.str());
    bytes.push_back('\0');
    ls2g_test::write_file(file.str(), bytes);
    REQUIRE_THROWS_MATCHES(
        load<ls2g_test::FakeLoad>(file.str()), std::runtime_error,
        MessageMatches(ContainsSubstring("unexpected byte(s)")));
}
