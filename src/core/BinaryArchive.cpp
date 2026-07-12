// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BinaryArchive.hpp"

#include <cstdio>   // std::remove, std::rename
#include <cstring>
#include <sstream>
#include <stdexcept>

namespace ls2g {

namespace {

const char MAGIC[4] = {'L', 'S', 'B', '2'};

// Binary format versions this build can load. Almost always just the
// current one: reading an older format requires dedicated migration code
// (the StateRes tuples are read structurally, with the current layout).
// If such code is ever written, add the older format number here.
const std::uint32_t SUPPORTED_BINARY_FORMATS[] = {BINARY_FORMAT_VERSION};

bool is_format_supported(std::uint32_t format)
{
    for (std::uint32_t f : SUPPORTED_BINARY_FORMATS) {
        if (f == format) return true;
    }
    return false;
}

std::string supported_formats_str()
{
    std::ostringstream oss;
    bool first = true;
    for (std::uint32_t f : SUPPORTED_BINARY_FORMATS) {
        if (!first) oss << ", ";
        oss << f;
        first = false;
    }
    return oss.str();
}

}  // anonymous namespace

BinaryArchive::BinaryArchive(const std::string & path, Mode mode, bool atomic_write):
    mode_(mode),
    path_(path),
    atomic_write_(atomic_write),
    tmp_path_(),
    committed_(false),
    file_size_(0)
{
    if (mode_ == Mode::Write) {
        if (atomic_write_) {
            // write to a temporary file: the destination is only replaced in
            // commit(), after the whole state has been written successfully.
            tmp_path_ = path_ + ".lsb_tmp";
            ofs_.open(tmp_path_, std::ios::binary | std::ios::out | std::ios::trunc);
            if (!ofs_.is_open()) {
                throw std::runtime_error(
                    "BinaryArchive: cannot open file for writing: '" + path_ +
                    "' (temporary file '" + tmp_path_ + "' could not be created)");
            }
        } else {
            // fast path (atomic=false): truncate + write the destination in
            // place, no rename; an interrupted save leaves a truncated file
            ofs_.open(path_, std::ios::binary | std::ios::out | std::ios::trunc);
            if (!ofs_.is_open()) {
                throw std::runtime_error("BinaryArchive: cannot open file for writing: '" + path_ + "'");
            }
        }
    } else {
        ifs_.open(path_, std::ios::binary | std::ios::in);
        if (!ifs_.is_open()) {
            throw std::runtime_error("BinaryArchive: cannot open file for reading: '" + path_ + "'");
        }
        // capture the file size once: every count / length field read from
        // the file is validated against it before any allocation.
        ifs_.seekg(0, std::ios::end);
        file_size_ = static_cast<std::uint64_t>(ifs_.tellg());
        ifs_.seekg(0, std::ios::beg);
    }
}

BinaryArchive::~BinaryArchive()
{
    if (mode_ == Mode::Write && atomic_write_ && !committed_) {
        // interrupted write (exception during serialization, or commit()
        // failed): drop the temporary file, the destination is untouched.
        if (ofs_.is_open()) ofs_.close();
        std::remove(tmp_path_.c_str());
    }
}

void BinaryArchive::commit()
{
    if (mode_ != Mode::Write) {
        throw std::runtime_error("BinaryArchive: commit() called on a read-mode archive for '" + path_ + "'");
    }
    ofs_.flush();
    if (!ofs_) {
        throw std::runtime_error("BinaryArchive: failed to write to file '" + path_ + "'");
    }
    ofs_.close();
    if (atomic_write_) {
        // atomically replace the destination. std::rename overwrites on POSIX
        // but fails on Windows when the destination exists, so remove it first
        // (failure to remove is ignored: rename below reports the real error).
        std::remove(path_.c_str());
        if (std::rename(tmp_path_.c_str(), path_.c_str()) != 0) {
            throw std::runtime_error(
                "BinaryArchive: could not move the temporary file '" + tmp_path_ +
                "' to its final destination '" + path_ + "'");
        }
    }
    committed_ = true;
}

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

std::uint64_t BinaryArchive::remaining_bytes()
{
    const std::uint64_t pos = static_cast<std::uint64_t>(ifs_.tellg());
    return (pos <= file_size_) ? (file_size_ - pos) : 0;
}

void BinaryArchive::require_count(std::uint64_t count, std::uint64_t elem_size)
{
    // division (not multiplication) so a huge corrupted count cannot
    // overflow back into an "acceptable" byte total
    if (elem_size == 0) return;
    if (count > remaining_bytes() / elem_size) {
        std::ostringstream oss;
        oss << "BinaryArchive: invalid file '" << path_ << "': it declares a block of "
            << count << " element(s) of " << elem_size << " byte(s), but only "
            << remaining_bytes() << " byte(s) remain in the file. "
            << "The file is truncated or corrupted.";
        throw std::runtime_error(oss.str());
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

void BinaryArchive::write_header(const std::string & type_tag,
                                 const std::string & v_major, const std::string & v_medium, const std::string & v_minor)
{
    write_magic();
    write_scalar(BINARY_FORMAT_VERSION);
    write_string(type_tag);
    write_string(v_major);
    write_string(v_medium);
    write_string(v_minor);
}

void BinaryArchive::check_header(const std::string & type_tag,
                                 const std::string & v_major, const std::string & v_medium, const std::string & v_minor)
{
    check_magic();
    // the header layout (magic, format, tag, writer version) is a stable
    // contract across binary format versions: it can be parsed even when
    // the format number does not match, to produce a precise error below.
    std::uint32_t file_format = 0;
    read_scalar(file_format);
    std::string file_tag, file_major, file_medium, file_minor;
    read_string(file_tag);
    read_string(file_major);
    read_string(file_medium);
    read_string(file_minor);
    if (!is_format_supported(file_format)) {
        std::ostringstream oss;
        oss << "BinaryArchive: incompatible file '" << path_ << "': it was written by lightsim2grid version "
            << file_major << "." << file_medium << "." << file_minor
            << " using binary format " << file_format
            << ", but the currently installed lightsim2grid (version "
            << v_major << "." << v_medium << "." << v_minor
            << ") only supports binary format(s) {" << supported_formats_str() << "}. "
            << "Re-save the file with a compatible lightsim2grid version (the binary format only "
            << "changes when the serialized layout changes, so versions sharing the same format "
            << "number can read each other's files).";
        throw std::runtime_error(oss.str());
    }
    if (file_tag != type_tag) {
        std::ostringstream oss;
        oss << "BinaryArchive: wrong object type in file '" << path_ << "': it contains a '"
            << file_tag << "', not a '" << type_tag << "'.";
        throw std::runtime_error(oss.str());
    }
}

void BinaryArchive::check_fully_consumed()
{
    const std::uint64_t left = remaining_bytes();
    if (left != 0) {
        std::ostringstream oss;
        oss << "BinaryArchive: invalid file '" << path_ << "': " << left
            << " unexpected byte(s) remain after the end of the data. "
            << "The file is corrupted, or was written by an incompatible build.";
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
    require_count(len, 1);
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
    require_count(n, 1);
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
    // each string takes at least its own 4-byte length prefix
    require_count(n, sizeof(std::uint32_t));
    v.clear();
    v.reserve(static_cast<std::size_t>(n));
    for (std::uint64_t i = 0; i < n; ++i) {
        std::string s;
        read_string(s);
        v.push_back(std::move(s));
    }
}

}  // namespace ls2g
