// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BinaryArchive.hpp"

#include <cstring>
#include <sstream>
#include <stdexcept>

namespace ls2g {

namespace {
const char MAGIC[4] = {'L', 'S', 'B', '1'};
}  // anonymous namespace

BinaryArchive::BinaryArchive(const std::string & path, Mode mode):
    mode_(mode),
    path_(path)
{
    if (mode_ == Mode::Write) {
        ofs_.open(path_, std::ios::binary | std::ios::out | std::ios::trunc);
        if (!ofs_.is_open()) {
            throw std::runtime_error("BinaryArchive: cannot open file for writing: '" + path_ + "'");
        }
    } else {
        ifs_.open(path_, std::ios::binary | std::ios::in);
        if (!ifs_.is_open()) {
            throw std::runtime_error("BinaryArchive: cannot open file for reading: '" + path_ + "'");
        }
    }
}

BinaryArchive::~BinaryArchive() = default;

void BinaryArchive::write_raw(const void * data, std::size_t nbytes)
{
    if (nbytes == 0) return;
    ofs_.write(static_cast<const char *>(data), static_cast<std::streamsize>(nbytes));
    if (!ofs_) {
        throw std::runtime_error("BinaryArchive: failed to write to file '" + path_ + "'");
    }
}

void BinaryArchive::read_raw(void * data, std::size_t nbytes)
{
    if (nbytes == 0) return;
    ifs_.read(static_cast<char *>(data), static_cast<std::streamsize>(nbytes));
    if (!ifs_) {
        throw std::runtime_error(
            "BinaryArchive: unexpected end of file (or read error) while reading '" + path_ +
            "'. The file is likely truncated or corrupted.");
    }
}

void BinaryArchive::write_magic()
{
    write_raw(MAGIC, sizeof(MAGIC));
}

void BinaryArchive::check_magic()
{
    char buf[sizeof(MAGIC)];
    read_raw(buf, sizeof(buf));
    if (std::memcmp(buf, MAGIC, sizeof(MAGIC)) != 0) {
        throw std::runtime_error(
            "BinaryArchive: invalid file '" + path_ +
            "': magic number mismatch. This is not a lightsim2grid binary file, or it is corrupted.");
    }
}

void BinaryArchive::write_header(const std::string & v_major, const std::string & v_medium, const std::string & v_minor)
{
    write_magic();
    write_string(v_major);
    write_string(v_medium);
    write_string(v_minor);
}

void BinaryArchive::check_header(const std::string & v_major, const std::string & v_medium, const std::string & v_minor)
{
    check_magic();
    std::string file_major, file_medium, file_minor;
    read_string(file_major);
    read_string(file_medium);
    read_string(file_minor);
    if (file_major != v_major || file_medium != v_medium || file_minor != v_minor) {
        std::ostringstream oss;
        oss << "BinaryArchive: version mismatch when loading '" << path_ << "'. "
            << "This file was saved with lightsim2grid version " << file_major << "." << file_medium << "." << file_minor
            << ", but the currently installed version is " << v_major << "." << v_medium << "." << v_minor
            << ". Binary files can only be reloaded with the exact same lightsim2grid version they were saved with "
            << "(unlike pickle, this format does not attempt any cross-version migration).";
        throw std::runtime_error(oss.str());
    }
}

void BinaryArchive::write_string(const std::string & s)
{
    std::uint32_t len = static_cast<std::uint32_t>(s.size());
    write_scalar(len);
    if (len) write_raw(s.data(), len);
}

void BinaryArchive::read_string(std::string & s)
{
    std::uint32_t len = 0;
    read_scalar(len);
    s.resize(len);
    if (len) read_raw(&s[0], len);
}

void BinaryArchive::write_bool_vector(const std::vector<bool> & v)
{
    std::uint64_t n = static_cast<std::uint64_t>(v.size());
    write_scalar(n);
    std::vector<std::uint8_t> buf(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) buf[i] = v[i] ? 1 : 0;
    if (!buf.empty()) write_raw(buf.data(), buf.size());
}

void BinaryArchive::read_bool_vector(std::vector<bool> & v)
{
    std::uint64_t n = 0;
    read_scalar(n);
    std::vector<std::uint8_t> buf(static_cast<std::size_t>(n));
    if (!buf.empty()) read_raw(buf.data(), buf.size());
    v.resize(static_cast<std::size_t>(n));
    for (std::size_t i = 0; i < buf.size(); ++i) v[i] = (buf[i] != 0);
}

void BinaryArchive::write_string_vector(const std::vector<std::string> & v)
{
    std::uint64_t n = static_cast<std::uint64_t>(v.size());
    write_scalar(n);
    for (const auto & s : v) write_string(s);
}

void BinaryArchive::read_string_vector(std::vector<std::string> & v)
{
    std::uint64_t n = 0;
    read_scalar(n);
    v.clear();
    v.reserve(static_cast<std::size_t>(n));
    for (std::uint64_t i = 0; i < n; ++i) {
        std::string s;
        read_string(s);
        v.push_back(std::move(s));
    }
}

}  // namespace ls2g
