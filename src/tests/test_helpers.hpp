// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LS2G_TEST_HELPERS_H
#define LS2G_TEST_HELPERS_H

// Shared helpers for the C++ unit tests (src/tests): a self-cleaning
// temporary file, raw file IO, and the synthetic serializable types used to
// exercise BinaryArchive without a real grid. C++14 only (project policy).

#include <cstdio>    // std::remove
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "Utils.hpp"  // ls2g::real_type, ls2g::cplx_type

namespace ls2g_test {

// A unique file path in the current working directory (ctest / CI run the
// binary from the build directory), removed on destruction together with the
// ".lsb_tmp" sibling a BinaryArchive atomic write may have left behind.
// mkdtemp/std::filesystem are avoided on purpose: C++14, no platform #ifdef.
class TempFile
{
    public:
        explicit TempFile(const std::string & suffix = ".lsb") {
            static int counter = 0;
            std::ostringstream oss;
            oss << "ls2g_unit_test_" << counter++ << suffix;
            path_ = oss.str();
            // in case a previous crashed run left files behind
            std::remove(path_.c_str());
            std::remove((path_ + ".lsb_tmp").c_str());
        }
        ~TempFile() {
            std::remove(path_.c_str());
            std::remove((path_ + ".lsb_tmp").c_str());
        }
        TempFile(const TempFile &) = delete;
        TempFile & operator=(const TempFile &) = delete;

        const std::string & str() const { return path_; }

    private:
        std::string path_;
};

inline bool file_exists(const std::string & path)
{
    std::ifstream f(path, std::ios::binary);
    return f.is_open();
}

inline std::vector<char> read_file(const std::string & path)
{
    std::ifstream f(path, std::ios::binary);
    return std::vector<char>((std::istreambuf_iterator<char>(f)),
                             std::istreambuf_iterator<char>());
}

inline void write_file(const std::string & path, const std::vector<char> & bytes)
{
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    f.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
}

// ---- synthetic serializable types ------------------------------------------
// They implement the same contract as the real containers (StateRes typedef +
// get_state / set_state / binary_type_tag) with none of their dependencies.

enum class FakeEnum : int { kOff = 0, kOn = 1, kAuto = 42 };

// One field per ValueArchiver specialization in BinaryArchive.hpp: arithmetic
// scalars, bool, cplx_type, enum, string, raw vectors (real / int / cplx),
// vector<bool>, vector<string>, vector<vector<T>> and a nested tuple.
struct FakeContainer
{
    using SubState = std::tuple<int, std::vector<ls2g::real_type> >;
    using StateRes = std::tuple<
        int,
        ls2g::real_type,
        bool,
        ls2g::cplx_type,
        FakeEnum,
        std::string,
        std::vector<ls2g::real_type>,
        std::vector<int>,
        std::vector<ls2g::cplx_type>,
        std::vector<bool>,
        std::vector<std::string>,
        std::vector<std::vector<ls2g::real_type> >,
        SubState
    >;

    StateRes state{};

    StateRes get_state() const { return state; }
    void set_state(StateRes & s) { state = s; }
    static const char * binary_type_tag() { return "FakeContainer"; }
};

inline FakeContainer make_reference_container()
{
    FakeContainer res;
    res.state = FakeContainer::StateRes(
        -7,
        3.141592653589793,
        true,
        ls2g::cplx_type(1.5, -2.5),
        FakeEnum::kAuto,
        "hello archive",
        {0., -1.5, 2.25e3},
        {1, -2, 3},
        {ls2g::cplx_type(0., 1.), ls2g::cplx_type(-1., 0.)},
        {true, false, true, true},
        {"", "one", "two words"},
        {{1., 2.}, {}, {3.}},
        FakeContainer::SubState(99, {4., 5., 6.})
    );
    return res;
}

// Two types with identical StateRes layouts but different tags: the
// LoadContainer-vs-StorageContainer confusion the type tag exists to reject.
struct FakeLoad
{
    using StateRes = std::tuple<int, std::vector<ls2g::real_type> >;
    StateRes state{};
    StateRes get_state() const { return state; }
    void set_state(StateRes & s) { state = s; }
    static const char * binary_type_tag() { return "FakeLoad"; }
};

struct FakeStorage
{
    using StateRes = std::tuple<int, std::vector<ls2g::real_type> >;
    StateRes state{};
    StateRes get_state() const { return state; }
    void set_state(StateRes & s) { state = s; }
    static const char * binary_type_tag() { return "FakeStorage"; }
};

// Opts into the optional post-read hook (like LSGrid does): detected by
// BinaryLoadFixup via the C++14 detection idiom and applied after the state
// is read, before set_state().
struct FakeWithFixup
{
    using StateRes = std::tuple<int, std::string>;
    StateRes state{};
    StateRes get_state() const { return state; }
    void set_state(StateRes & s) { state = s; }
    static const char * binary_type_tag() { return "FakeWithFixup"; }
    static void fixup_binary_state(StateRes & s) { std::get<1>(s) = "fixed up"; }
};

}  // namespace ls2g_test

#endif  // LS2G_TEST_HELPERS_H
