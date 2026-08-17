// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Low-level BinaryArchive API: primitive round-trips, header / magic / format
// / tag validation, every require_count path, truncation and trailing-byte
// detection, and the atomic temp-file commit / rollback lifecycle.

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "BinaryArchive.hpp"
#include "test_helpers.hpp"

using ls2g::BinaryArchive;
using ls2g_test::TempFile;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;

namespace {

// most tests start from a freshly written file
void with_written(const std::string & path, void (*writer)(BinaryArchive &))
{
    BinaryArchive ar(path, BinaryArchive::Mode::Write);
    writer(ar);
    ar.commit();
}

}  // anonymous namespace

TEST_CASE("scalars round-trip", "[BinaryArchive]")
{
    TempFile file;
    with_written(file.str(), [](BinaryArchive & ar) {
        ar.write_scalar(std::int32_t(-123));
        ar.write_scalar(std::uint64_t(1) << 60);
        ar.write_scalar(2.75);
        ar.write_scalar(ls2g::cplx_type(1., -1.));
    });

    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
    std::int32_t i = 0;
    std::uint64_t u = 0;
    double d = 0.;
    ls2g::cplx_type c;
    ar.read_scalar(i);
    ar.read_scalar(u);
    ar.read_scalar(d);
    ar.read_scalar(c);
    CHECK(i == -123);
    CHECK(u == (std::uint64_t(1) << 60));
    CHECK(d == 2.75);
    CHECK(c == ls2g::cplx_type(1., -1.));
    CHECK_NOTHROW(ar.check_fully_consumed());
}

TEST_CASE("strings round-trip, including empty and embedded NUL", "[BinaryArchive]")
{
    TempFile file;
    const std::string with_nul("a\0b", 3);
    {
        BinaryArchive ar(file.str(), BinaryArchive::Mode::Write);
        ar.write_string("");
        ar.write_string("plain");
        ar.write_string(with_nul);
        ar.commit();
    }
    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
    std::string a, b, c;
    ar.read_string(a);
    ar.read_string(b);
    ar.read_string(c);
    CHECK(a == "");
    CHECK(b == "plain");
    CHECK(c == with_nul);
    CHECK_NOTHROW(ar.check_fully_consumed());
}

TEST_CASE("bool / string / raw vectors round-trip, including empty", "[BinaryArchive]")
{
    TempFile file;
    const std::vector<bool> bools{true, false, true};
    const std::vector<std::string> strings{"", "x", "yz"};
    const std::vector<double> doubles{-1., 0., 1.5};
    const std::vector<int> ints{7, -8};
    const std::vector<ls2g::cplx_type> cplxs{ls2g::cplx_type(0., 2.)};
    {
        BinaryArchive ar(file.str(), BinaryArchive::Mode::Write);
        ar.write_bool_vector(bools);
        ar.write_bool_vector({});
        ar.write_string_vector(strings);
        ar.write_string_vector({});
        ar.write_vector_raw(doubles);
        ar.write_vector_raw(std::vector<double>{});
        ar.write_vector_raw(ints);
        ar.write_vector_raw(cplxs);
        ar.commit();
    }
    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
    std::vector<bool> rb, rb_empty;
    std::vector<std::string> rs, rs_empty;
    std::vector<double> rd, rd_empty;
    std::vector<int> ri;
    std::vector<ls2g::cplx_type> rc;
    ar.read_bool_vector(rb);
    ar.read_bool_vector(rb_empty);
    ar.read_string_vector(rs);
    ar.read_string_vector(rs_empty);
    ar.read_vector_raw(rd);
    ar.read_vector_raw(rd_empty);
    ar.read_vector_raw(ri);
    ar.read_vector_raw(rc);
    CHECK(rb == bools);
    CHECK(rb_empty.empty());
    CHECK(rs == strings);
    CHECK(rs_empty.empty());
    CHECK(rd == doubles);
    CHECK(rd_empty.empty());
    CHECK(ri == ints);
    CHECK(rc == cplxs);
    CHECK_NOTHROW(ar.check_fully_consumed());
}

TEST_CASE("header round-trips; writer version strings are not compared", "[BinaryArchive]")
{
    TempFile file;
    {
        BinaryArchive ar(file.str(), BinaryArchive::Mode::Write);
        ar.write_header("SomeTag", "0", "14", "0");
        ar.commit();
    }
    {
        BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
        CHECK_NOTHROW(ar.check_header("SomeTag", "0", "14", "0"));
        CHECK_NOTHROW(ar.check_fully_consumed());
    }
    // a build with a different version but the same binary format reads the file
    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
    CHECK_NOTHROW(ar.check_header("SomeTag", "99", "0", "1"));
}

TEST_CASE("garbage magic number is rejected", "[BinaryArchive]")
{
    TempFile file;
    ls2g_test::write_file(file.str(), {'G', 'A', 'R', 'B', 'A', 'G', 'E', '!'});
    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
    REQUIRE_THROWS_MATCHES(
        ar.check_magic(), std::runtime_error,
        MessageMatches(ContainsSubstring("magic number mismatch")));
}

TEST_CASE("unsupported binary format version is rejected with a precise message", "[BinaryArchive]")
{
    TempFile file;
    {
        // hand-crafted header claiming format 999: magic is valid, so the
        // header must still parse and the error must cite both versions
        BinaryArchive ar(file.str(), BinaryArchive::Mode::Write);
        ar.write_magic();
        ar.write_scalar(std::uint32_t(999));
        ar.write_string("SomeTag");
        ar.write_string("12");
        ar.write_string("34");
        ar.write_string("56");
        ar.commit();
    }
    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
    REQUIRE_THROWS_MATCHES(
        ar.check_header("SomeTag", "0", "14", "0"), std::runtime_error,
        MessageMatches(ContainsSubstring("binary format 999") &&
                       ContainsSubstring("12.34.56") &&
                       ContainsSubstring("only supports binary format(s)")));
}

TEST_CASE("wrong type tag is rejected", "[BinaryArchive]")
{
    TempFile file;
    with_written(file.str(), [](BinaryArchive & ar) {
        ar.write_header("LoadContainer", "0", "14", "0");
    });
    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
    REQUIRE_THROWS_MATCHES(
        ar.check_header("StorageContainer", "0", "14", "0"), std::runtime_error,
        MessageMatches(ContainsSubstring("wrong object type") &&
                       ContainsSubstring("'LoadContainer'") &&
                       ContainsSubstring("'StorageContainer'")));
}

TEST_CASE("corrupt strings are escaped and truncated in error messages", "[BinaryArchive]")
{
    SECTION("non-printable bytes become \\x escapes") {
        TempFile file;
        with_written(file.str(), [](BinaryArchive & ar) {
            ar.write_header(std::string("bad\x01tag\xff", 8), "0", "14", "0");
        });
        BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
        REQUIRE_THROWS_MATCHES(
            ar.check_header("SomeTag", "0", "14", "0"), std::runtime_error,
            MessageMatches(ContainsSubstring("\\x01") && ContainsSubstring("\\xff")));
    }
    SECTION("over-long strings are truncated with a length note") {
        TempFile file;
        const std::string long_tag(200, 'A');
        {
            BinaryArchive ar(file.str(), BinaryArchive::Mode::Write);
            ar.write_header(long_tag, "0", "14", "0");
            ar.commit();
        }
        BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
        REQUIRE_THROWS_MATCHES(
            ar.check_header("SomeTag", "0", "14", "0"), std::runtime_error,
            MessageMatches(ContainsSubstring("... (200 chars)") &&
                           !ContainsSubstring(long_tag)));
    }
}

TEST_CASE("require_count", "[BinaryArchive]")
{
    TempFile file;
    // 16 bytes of payload, nothing consumed yet
    ls2g_test::write_file(file.str(), std::vector<char>(16, '\0'));
    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);

    SECTION("counts that fit do not throw") {
        CHECK_NOTHROW(ar.require_count(16, 1));
        CHECK_NOTHROW(ar.require_count(2, 8));
        CHECK_NOTHROW(ar.require_count(0, 8));
    }
    SECTION("one element too many throws before any allocation") {
        REQUIRE_THROWS_MATCHES(
            ar.require_count(17, 1), std::runtime_error,
            MessageMatches(ContainsSubstring("truncated or corrupted")));
        REQUIRE_THROWS_AS(ar.require_count(3, 8), std::runtime_error);
    }
    SECTION("huge counts cannot overflow back into acceptance") {
        // 2^61 * 8 == 2^64 wraps to 0 with a multiplication-based check
        REQUIRE_THROWS_AS(ar.require_count(std::uint64_t(1) << 61, 8), std::runtime_error);
        REQUIRE_THROWS_AS(ar.require_count(~std::uint64_t(0), 8), std::runtime_error);
    }
    SECTION("elem_size 0 never throws") {
        CHECK_NOTHROW(ar.require_count(~std::uint64_t(0), 0));
    }
    SECTION("the check is relative to the current read position") {
        std::uint64_t dummy = 0;
        ar.read_scalar(dummy);  // consume 8 of the 16 bytes
        CHECK_NOTHROW(ar.require_count(8, 1));
        REQUIRE_THROWS_AS(ar.require_count(9, 1), std::runtime_error);
    }
}

TEST_CASE("trailing bytes are treated as corruption", "[BinaryArchive]")
{
    TempFile file;
    with_written(file.str(), [](BinaryArchive & ar) {
        ar.write_scalar(std::int32_t(1));
        ar.write_scalar(std::uint8_t(0));  // the trailing byte
    });
    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
    std::int32_t v = 0;
    ar.read_scalar(v);
    REQUIRE_THROWS_MATCHES(
        ar.check_fully_consumed(), std::runtime_error,
        MessageMatches(ContainsSubstring("1 unexpected byte(s)")));
}

TEST_CASE("reading past the end of a truncated file throws", "[BinaryArchive]")
{
    TempFile file;
    ls2g_test::write_file(file.str(), {'\x01', '\x02', '\x03', '\x04'});
    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
    double d = 0.;  // needs 8 bytes, only 4 available
    REQUIRE_THROWS_MATCHES(
        ar.read_scalar(d), std::runtime_error,
        MessageMatches(ContainsSubstring("unexpected end of file")));
}

TEST_CASE("constructor failures name the file", "[BinaryArchive]")
{
    SECTION("reading a nonexistent file") {
        REQUIRE_THROWS_MATCHES(
            BinaryArchive("no_such_file.lsb", BinaryArchive::Mode::Read),
            std::runtime_error,
            MessageMatches(ContainsSubstring("cannot open file for reading") &&
                           ContainsSubstring("no_such_file.lsb")));
    }
    SECTION("atomic write into a nonexistent directory names the temp file") {
        REQUIRE_THROWS_MATCHES(
            BinaryArchive("no_such_dir_ls2g/x.lsb", BinaryArchive::Mode::Write),
            std::runtime_error,
            MessageMatches(ContainsSubstring("cannot open file for writing") &&
                           ContainsSubstring(".lsb_tmp")));
    }
    SECTION("non-atomic write into a nonexistent directory") {
        REQUIRE_THROWS_MATCHES(
            BinaryArchive("no_such_dir_ls2g/x.lsb", BinaryArchive::Mode::Write, false),
            std::runtime_error,
            MessageMatches(ContainsSubstring("cannot open file for writing")));
    }
}

TEST_CASE("commit() on a read-mode archive throws", "[BinaryArchive]")
{
    TempFile file;
    ls2g_test::write_file(file.str(), {'\0'});
    BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
    REQUIRE_THROWS_MATCHES(
        ar.commit(), std::runtime_error,
        MessageMatches(ContainsSubstring("read-mode")));
}

TEST_CASE("atomic write lifecycle", "[BinaryArchive]")
{
    TempFile file;
    const std::string tmp = file.str() + ".lsb_tmp";

    SECTION("data goes to the temp file, the destination only appears on commit") {
        BinaryArchive ar(file.str(), BinaryArchive::Mode::Write);
        ar.write_scalar(std::int32_t(1));
        CHECK(ls2g_test::file_exists(tmp));
        CHECK_FALSE(ls2g_test::file_exists(file.str()));
        ar.commit();
        CHECK_FALSE(ls2g_test::file_exists(tmp));
        CHECK(ls2g_test::file_exists(file.str()));
    }

    SECTION("an interrupted write leaves a previously good file untouched") {
        // a good file first
        with_written(file.str(), [](BinaryArchive & ar) {
            ar.write_scalar(std::int32_t(42));
        });
        const std::vector<char> before = ls2g_test::read_file(file.str());
        {
            // then a write abandoned before commit() (as after an exception
            // during serialization): rollback in the destructor
            BinaryArchive ar(file.str(), BinaryArchive::Mode::Write);
            ar.write_scalar(std::int32_t(-1));
            ar.write_scalar(std::int32_t(-2));
        }
        CHECK_FALSE(ls2g_test::file_exists(tmp));
        CHECK(ls2g_test::read_file(file.str()) == before);
    }
}

TEST_CASE("non-atomic write goes straight to the destination", "[BinaryArchive]")
{
    TempFile file;
    const std::string tmp = file.str() + ".lsb_tmp";

    SECTION("no temp file is involved") {
        BinaryArchive ar(file.str(), BinaryArchive::Mode::Write, false);
        ar.write_scalar(std::int32_t(1));
        CHECK_FALSE(ls2g_test::file_exists(tmp));
        CHECK(ls2g_test::file_exists(file.str()));
        ar.commit();
        CHECK(ls2g_test::file_exists(file.str()));
    }

    SECTION("an abandoned partial file remains and is rejected as truncated at load time") {
        {
            // declares a 5-element block whose data is never written
            BinaryArchive ar(file.str(), BinaryArchive::Mode::Write, false);
            ar.write_scalar(std::uint64_t(5));
            // no commit(): simulates an interrupted save
        }
        CHECK(ls2g_test::file_exists(file.str()));
        BinaryArchive ar(file.str(), BinaryArchive::Mode::Read);
        std::vector<double> v;
        REQUIRE_THROWS_AS(ar.read_vector_raw(v), std::runtime_error);
    }
}
