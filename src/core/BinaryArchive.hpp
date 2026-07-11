// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BINARYARCHIVE_H
#define BINARYARCHIVE_H

// Fast, additive binary serialization for any class exposing a `StateRes`
// tuple + get_state()/set_state() (the same contract already used for
// python pickle support, see pickle_helpers.hpp). This is NOT meant to
// replace pickle: it trades portability for speed (no python-level
// marshalling, raw contiguous writes for vector data).
//
// Robustness guarantees (all failures are std::runtime_error naming the
// file, never a crash, a std::bad_alloc or an out-of-memory situation):
//  - a magic number rejects files that are not lightsim2grid binary files;
//  - a *binary format version* (see BINARY_FORMAT_VERSION below) decides
//    compatibility across lightsim2grid releases;
//  - a type tag rejects loading a file into the wrong class (eg a
//    LoadContainer file into StorageContainer.load_binary, whose StateRes
//    layouts are otherwise identical);
//  - every count/length field read from the file is checked against the
//    actual file size *before* any allocation, so a corrupted count cannot
//    trigger a multi-gigabyte allocation;
//  - after the state is fully read, the file must be fully consumed
//    (trailing bytes are treated as corruption);
//  - writes go to a temporary file that is atomically renamed over the
//    destination only on success, so a crash / full disk mid-save never
//    destroys a previously good file.
// NOT covered (possible future work): a payload checksum (bit flips inside
// numeric data that keep all counts consistent load silently), and an
// ABI descriptor (endianness / sizeof(real_type)): files are assumed to be
// written and read by builds with the same data layout.
//
// C++14 only (see project policy, eg LSGrid.hpp): no `if constexpr`, no
// `std::apply`, no fold expressions. The recursive walk over an arbitrary
// StateRes tuple is done via partial-specialization of the `ValueArchiver`
// class template (not free-function SFINAE overloads): class template
// partial specialization resolution does not depend on declaration order
// or two-phase-lookup/ADL subtleties the way a set of mutually-recursive
// free function template overloads would, which matters here since a
// StateRes tuple can nest other StateRes tuples several levels deep.

#include <cstdint>
#include <fstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "ls2g_api.hpp"
#include "Utils.hpp"  // real_type, cplx_type

namespace ls2g {

// Version of the binary *format* itself, decoupled from the lightsim2grid
// package version: all lightsim2grid releases sharing the same format number
// can read each other's binary files (most releases do not touch the
// serialized layout at all). The lightsim2grid version that wrote a file is
// still stored in its header, but only for diagnostics / error messages.
//
// >>> BUMP THIS NUMBER whenever anything about the serialized layout
// >>> changes: any StateRes tuple (fields added / removed / reordered /
// >>> retyped, in LSGrid or any container), or the encoding in this file.
// A bumped format makes older files fail to load with a clear error instead
// of being silently mis-read. If a future version wants to keep reading an
// older format, add that format number to SUPPORTED_BINARY_FORMATS (in
// BinaryArchive.cpp) together with the required migration code.
constexpr std::uint32_t BINARY_FORMAT_VERSION = 1;

class LS2G_API BinaryArchive
{
    public:
        enum class Mode { Write, Read };

        // Write mode writes to a temporary file (path + ".lsb_tmp"); the
        // data only replaces `path` when commit() is called after a fully
        // successful write. Read mode opens `path` directly and captures
        // its size for the bounds checks below.
        BinaryArchive(const std::string & path, Mode mode);
        ~BinaryArchive();
        BinaryArchive(const BinaryArchive &) = delete;
        BinaryArchive & operator=(const BinaryArchive &) = delete;

        // Write mode only: flush + close the temporary file and atomically
        // rename it over the destination path. Throws std::runtime_error on
        // any failure (the destination is left untouched in that case). If
        // commit() is never called (eg an exception during the write), the
        // destructor removes the temporary file and the destination keeps
        // its previous content.
        void commit();

        // low level: the only methods that touch the underlying stream.
        // Both throw std::runtime_error (including the file path) on any
        // stream failure, which uniformly covers "file does not exist",
        // "corrupted file" and "truncated file" (unexpected EOF) cases.
        void write_raw(const void * data, std::size_t nbytes);
        void read_raw(void * data, std::size_t nbytes);

        // magic number ("LSB2"), to catch garbage files early
        void write_magic();
        void check_magic();  // throws std::runtime_error on mismatch

        // Full header: magic, BINARY_FORMAT_VERSION (uint32), type tag
        // (length-prefixed string, eg "LoadContainer"), then 3
        // length-prefixed strings with the lightsim2grid version that wrote
        // the file. This header layout is a stable contract across all
        // binary format versions (so any version can at least *parse* the
        // header of any file and produce a precise error message); only the
        // state payload after the header is governed by the format number.
        //
        // check_header() throws std::runtime_error when the file's format
        // version is not supported by this build, or when the type tag does
        // not match. The writer version strings are NOT compared: two
        // different lightsim2grid versions sharing the same binary format
        // read each other's files.
        //
        // Version strings are passed in explicitly by the caller (never read
        // from the VERSION_MAJOR/MEDIUM/MINOR macros in this shared header):
        // those macros are only defined via target_compile_definitions on
        // the python bindings target, not on lightsim2grid_core, so
        // referencing them here would silently embed different values (even
        // risk an ODR violation) depending on which translation unit
        // instantiates this header first. Callers (each container's own
        // .cpp, or binary_helpers.hpp on the bindings side) supply the
        // macros themselves, each from a single, consistent translation unit.
        void write_header(const std::string & type_tag,
                          const std::string & v_major, const std::string & v_medium, const std::string & v_minor);
        void check_header(const std::string & type_tag,
                          const std::string & v_major, const std::string & v_medium, const std::string & v_minor);

        // Read mode only: throws std::runtime_error if the whole file has
        // not been consumed (trailing bytes = corrupted / wrong file).
        void check_fully_consumed();

        void write_string(const std::string & s);
        void read_string(std::string & s);

        // std::vector<bool> is bit-packed (not contiguous): converts to/from
        // std::vector<uint8_t> (1 byte per bool) rather than touching raw
        // storage.
        void write_bool_vector(const std::vector<bool> & v);
        void read_bool_vector(std::vector<bool> & v);

        // count + length-prefixed strings
        void write_string_vector(const std::vector<std::string> & v);
        void read_string_vector(std::vector<std::string> & v);

        // any trivially-copyable scalar (arithmetic incl. bool, or cplx_type)
        template<typename T>
        void write_scalar(const T & v) { write_raw(&v, sizeof(T)); }
        template<typename T>
        void read_scalar(T & v) { read_raw(&v, sizeof(T)); }

        // std::vector<T> of a trivially-copyable, contiguous element type
        // (real_type, int, cplx_type, ...): count + raw bytes straight from
        // / into .data(), no per-element conversion.
        template<typename T>
        void write_vector_raw(const std::vector<T> & v) {
            std::uint64_t n = static_cast<std::uint64_t>(v.size());
            write_scalar(n);
            if (n) write_raw(v.data(), static_cast<std::size_t>(n) * sizeof(T));
        }
        template<typename T>
        void read_vector_raw(std::vector<T> & v) {
            std::uint64_t n = 0;
            read_scalar(n);
            // reject a corrupted count before allocating anything
            require_count(n, sizeof(T));
            v.resize(static_cast<std::size_t>(n));
            if (n) read_raw(v.data(), static_cast<std::size_t>(n) * sizeof(T));
        }

        // Read mode: throws std::runtime_error if the file does not contain
        // at least `count * elem_size` more bytes after the current read
        // position (overflow-safe: does the check by division). Called
        // before every count-driven allocation so ill-formed counts turn
        // into clean errors instead of std::bad_alloc / OOM.
        void require_count(std::uint64_t count, std::uint64_t elem_size);

    private:
        std::uint64_t remaining_bytes();  // read mode: bytes after the current position

        std::ofstream ofs_;
        std::ifstream ifs_;
        Mode mode_;
        std::string path_;
        std::string tmp_path_;         // write mode: temporary file actually written
        bool committed_;               // write mode: commit() completed
        std::uint64_t file_size_;      // read mode: total size of the file
};

namespace archive_detail {

// C++14 substitute for std::void_t (C++17)
template<typename... Ts> struct make_void { using type = void; };
template<typename... Ts> using void_t = typename make_void<Ts...>::type;

// detects any std::tuple<...> (or more generally any type with a
// std::tuple_size specialization) via the standard pre-C++17
// detection-idiom pattern.
template<typename T, typename = void>
struct is_tuple_like : std::false_type {};
template<typename T>
struct is_tuple_like<T, void_t<decltype(std::tuple_size<T>::value)> > : std::true_type {};

// element types safe to write/read as raw contiguous bytes inside a
// std::vector<T> (excludes bool: std::vector<bool> is bit-packed and has
// no .data(), and is handled by BinaryArchive::write_bool_vector instead).
template<typename T>
struct is_raw_vector_elem : std::integral_constant<bool,
    (std::is_arithmetic<T>::value && !std::is_same<T, bool>::value) || std::is_same<T, cplx_type>::value> {};

// ---- ValueArchiver: recursive dispatcher, one partial specialization per
// StateRes field shape actually found in the codebase. Adding a field type
// not covered here leaves the primary (undefined) template selected, which
// fails to compile loudly rather than silently mis-serializing.
template<typename T, typename Enable = void>
struct ValueArchiver;  // intentionally incomplete: every real field type below has a specialization

// arithmetic scalars (real_type, int, bool, ...)
template<typename T>
struct ValueArchiver<T, typename std::enable_if<std::is_arithmetic<T>::value>::type> {
    static void write(BinaryArchive & ar, const T & v) { ar.write_scalar(v); }
    static void read(BinaryArchive & ar, T & v) { ar.read_scalar(v); }
};

// cplx_type (std::complex<real_type>), standard-guaranteed layout
// compatible with real_type[2] -- safe as a raw scalar for a same-layout
// round trip (no cross-platform guarantee needed: the binary format is
// only supported on builds sharing the same data layout, see the
// robustness notes at the top of this file).
template<>
struct ValueArchiver<cplx_type> {
    static void write(BinaryArchive & ar, const cplx_type & v) { ar.write_scalar(v); }
    static void read(BinaryArchive & ar, cplx_type & v) { ar.read_scalar(v); }
};

// enum types (eg AlgorithmType): stored as their underlying integer type
template<typename T>
struct ValueArchiver<T, typename std::enable_if<std::is_enum<T>::value>::type> {
    using Under = typename std::underlying_type<T>::type;
    static void write(BinaryArchive & ar, const T & v) { ar.write_scalar(static_cast<Under>(v)); }
    static void read(BinaryArchive & ar, T & v) { Under u{}; ar.read_scalar(u); v = static_cast<T>(u); }
};

template<>
struct ValueArchiver<std::string> {
    static void write(BinaryArchive & ar, const std::string & v) { ar.write_string(v); }
    static void read(BinaryArchive & ar, std::string & v) { ar.read_string(v); }
};

template<>
struct ValueArchiver<std::vector<bool> > {
    static void write(BinaryArchive & ar, const std::vector<bool> & v) { ar.write_bool_vector(v); }
    static void read(BinaryArchive & ar, std::vector<bool> & v) { ar.read_bool_vector(v); }
};

template<>
struct ValueArchiver<std::vector<std::string> > {
    static void write(BinaryArchive & ar, const std::vector<std::string> & v) { ar.write_string_vector(v); }
    static void read(BinaryArchive & ar, std::vector<std::string> & v) { ar.read_string_vector(v); }
};

// std::vector<std::vector<T>> (only TrafoContainer's rx_corr_alpha/pct):
// outer count, then each inner vector via the raw-vector path.
template<typename T>
struct ValueArchiver<std::vector<std::vector<T> > > {
    static void write(BinaryArchive & ar, const std::vector<std::vector<T> > & v) {
        std::uint64_t n = static_cast<std::uint64_t>(v.size());
        ar.write_scalar(n);
        for (const auto & inner : v) ValueArchiver<std::vector<T> >::write(ar, inner);
    }
    static void read(BinaryArchive & ar, std::vector<std::vector<T> > & v) {
        std::uint64_t n = 0;
        ar.read_scalar(n);
        // each inner vector takes at least its own 8-byte count
        ar.require_count(n, sizeof(std::uint64_t));
        v.resize(static_cast<std::size_t>(n));
        for (auto & inner : v) ValueArchiver<std::vector<T> >::read(ar, inner);
    }
};

// std::vector<T> of raw-copyable T (real_type, int, cplx_type, ...)
template<typename T>
struct ValueArchiver<std::vector<T>, typename std::enable_if<is_raw_vector_elem<T>::value>::type> {
    static void write(BinaryArchive & ar, const std::vector<T> & v) { ar.write_vector_raw(v); }
    static void read(BinaryArchive & ar, std::vector<T> & v) { ar.read_vector_raw(v); }
};

// any std::tuple<...> (a StateRes, possibly nesting other StateRes tuples):
// recurse field-by-field via the index_sequence/dummy-array fold idiom
// already used in NRSystem.hpp (C++14-safe, no std::apply).
template<typename Tuple>
struct ValueArchiver<Tuple, typename std::enable_if<is_tuple_like<Tuple>::value>::type> {
    template<std::size_t... Is>
    static void write_impl(BinaryArchive & ar, const Tuple & t, std::index_sequence<Is...>) {
        int dummy[] = { 0, (ValueArchiver<typename std::tuple_element<Is, Tuple>::type>::write(ar, std::get<Is>(t)), 0)... };
        (void)dummy;
    }
    template<std::size_t... Is>
    static void read_impl(BinaryArchive & ar, Tuple & t, std::index_sequence<Is...>) {
        int dummy[] = { 0, (ValueArchiver<typename std::tuple_element<Is, Tuple>::type>::read(ar, std::get<Is>(t)), 0)... };
        (void)dummy;
    }
    static void write(BinaryArchive & ar, const Tuple & t) {
        write_impl(ar, t, std::make_index_sequence<std::tuple_size<Tuple>::value>{});
    }
    static void read(BinaryArchive & ar, Tuple & t) {
        read_impl(ar, t, std::make_index_sequence<std::tuple_size<Tuple>::value>{});
    }
};

// ---- optional per-class fixup applied after a state has been read from a
// binary file, before set_state(). Detected via the standard C++14
// detection idiom: a class opting in declares
//     static void fixup_binary_state(StateRes & state);
// Today only LSGrid uses it (its StateRes embeds the writing build's
// version strings, which its pickle-shared set_state() checks against the
// *current* build -- fine for pickle, but it would wrongly reject a
// cross-version binary load that the format number allows, so the fixup
// rewrites those fields to the current build's values).
template<typename T, typename = void>
struct BinaryLoadFixup {
    static void apply(typename T::StateRes &) {}
};
template<typename T>
struct BinaryLoadFixup<T, void_t<decltype(T::fixup_binary_state(std::declval<typename T::StateRes &>()))> > {
    static void apply(typename T::StateRes & state) { T::fixup_binary_state(state); }
};

}  // namespace archive_detail

// Public entry points into the dispatcher above.
template<typename T>
void archive_write_value(BinaryArchive & ar, const T & v) {
    archive_detail::ValueArchiver<typename std::decay<T>::type>::write(ar, v);
}
template<typename T>
void archive_read_value(BinaryArchive & ar, T & v) {
    archive_detail::ValueArchiver<typename std::decay<T>::type>::read(ar, v);
}

// Top-level helpers reused by every serializable class's save_binary()/
// load_binary() (LSGrid + every element container with its own StateRes) --
// the single shared implementation the per-class methods delegate to, so
// none of them duplicate any serialization logic. T must additionally
// expose `static const char * binary_type_tag()` (its class name), written
// into / checked against the file header.
template<typename T>
void save_binary_generic(const T & obj, const std::string & path,
                          const std::string & v_major, const std::string & v_medium, const std::string & v_minor) {
    BinaryArchive ar(path, BinaryArchive::Mode::Write);
    ar.write_header(T::binary_type_tag(), v_major, v_medium, v_minor);
    typename T::StateRes state = obj.get_state();
    archive_write_value(ar, state);
    ar.commit();
}

template<typename T>
T load_binary_generic(const std::string & path,
                       const std::string & v_major, const std::string & v_medium, const std::string & v_minor) {
    BinaryArchive ar(path, BinaryArchive::Mode::Read);
    ar.check_header(T::binary_type_tag(), v_major, v_medium, v_minor);
    typename T::StateRes state{};
    archive_read_value(ar, state);
    ar.check_fully_consumed();
    archive_detail::BinaryLoadFixup<T>::apply(state);
    T res{};
    res.set_state(state);
    return res;
}

}  // namespace ls2g

#endif  // BINARYARCHIVE_H
