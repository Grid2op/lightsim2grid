// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ALGORITHM_REGISTRY_H
#define ALGORITHM_REGISTRY_H

#include <cstddef>
#include <exception>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "powerflow_algorithm/BaseAlgo.hpp"
#include "Ls2gAbiTag.hpp"

namespace ls2g {

// Longest accepted solver name (see is_valid_solver_name). Matches the length at
// which printable() truncates, so a name always fits whole in an error message.
constexpr std::size_t SOLVER_NAME_MAX_LEN = 64;

/**
 * Is `name` an acceptable solver name?
 *
 *   first character  : A-Z a-z _
 *   other characters : A-Z a-z 0-9 _ .
 *   length           : 1 to SOLVER_NAME_MAX_LEN characters
 *
 * A solver name is an identity, not just a label: it is written into every
 * serialized grid (see LSGrid::StateRes) and is what selects the solver again on
 * load. Restricting it to this ASCII subset is what makes that identity
 * trustworthy and safe to display:
 *
 *  - **no non-ASCII**: otherwise a plugin could register a name that is visually
 *    identical to a built-in one (say "NR_SparseLU" with a Cyrillic 'a') and
 *    masquerade as it in every error message, doc and available_algorithm_names()
 *    listing, while being a different registry key;
 *  - **no control characters**: names end up in exception messages and logs, so a
 *    newline would allow log injection and an ESC byte terminal-escape injection;
 *  - **bounded length**: bounds what goes into a file, into an error message and
 *    into a plugin's diagnostic buffer;
 *  - **'.' allowed (but not first)**: lets plugin authors namespace their solvers
 *    ("acme.NR_fast"), which is the recommended way to avoid name collisions --
 *    see docs/solver_plugin.rst;
 *  - **digits allowed (but not first)**: "solver_v2" is natural, while a leading
 *    digit would let a name look like a number.
 *
 * The set is deliberately minimal: it can be widened later without breaking
 * anyone, whereas narrowing it could not.
 */
inline bool is_valid_solver_name(const std::string& name) noexcept {
    if (name.empty() || name.size() > SOLVER_NAME_MAX_LEN) return false;
    // NB explicit ranges, not std::isalpha & co: those are locale-dependent and
    // are UB for negative char values.
    const unsigned char first = static_cast<unsigned char>(name[0]);
    const bool first_ok = (first >= 'A' && first <= 'Z') ||
                          (first >= 'a' && first <= 'z') ||
                          (first == '_');
    if (!first_ok) return false;
    for (std::size_t i = 1; i < name.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(name[i]);
        const bool ok = (c >= 'A' && c <= 'Z') ||
                        (c >= 'a' && c <= 'z') ||
                        (c >= '0' && c <= '9') ||
                        (c == '_') || (c == '.');
        if (!ok) return false;
    }
    return true;
}

// Human-readable statement of the rule above, for error messages.
inline std::string solver_name_rule() {
    return "a solver name must be 1 to " + std::to_string(SOLVER_NAME_MAX_LEN) +
           " characters, start with an ASCII letter or '_', and contain only ASCII "
           "letters, digits, '_' or '.'";
}

class LS2G_API AlgorithmRegistry {
public:
    using Factory = std::function<std::unique_ptr<BaseAlgo>()>;
    using FactoryMap = std::unordered_map<std::string, Factory>;

    static AlgorithmRegistry& instance();

    // `caller_tag` defaults to `ls2g_current_abi_tag()` evaluated at the call
    // site: since it is compiled fresh into whichever translation unit calls
    // this, the default captures that TU's own Eigen SIMD/alignment settings.
    // Throws if it disagrees with lightsim2grid_core's own tag -- see
    // Ls2gAbiTag.hpp. Used for the built-in solvers (BuiltinSolversRegistration).
    void register_solver(const std::string& name, Factory f,
                          Ls2gAbiTag caller_tag = ls2g_current_abi_tag());

    // Register a whole batch of solvers atomically: either every name in
    // `batch` is added, or none is and the registry is left exactly as it was.
    // This is what a plugin registers through (see LS2G_PLUGIN_ENTRY below) so
    // that a plugin exposing several solvers can never half-register when a
    // later name collides with an already-registered one. Throws (leaving the
    // registry unchanged) on an ABI-tag mismatch or a name already present.
    void register_all(FactoryMap batch, Ls2gAbiTag caller_tag = ls2g_current_abi_tag());

    std::unique_ptr<BaseAlgo> make(const std::string& name) const;
    bool is_registered(const std::string& name) const;
    std::vector<std::string> available_algorithm_names() const;

private:
    AlgorithmRegistry() = default;
    FactoryMap _factories;
};

// Staging area a plugin fills in with its solvers before they are committed to
// the real registry. Handed to the plugin's registration function; `add()`
// only stages (and rejects a name the same plugin already staged) -- the actual
// commit happens atomically once the whole function has run without throwing.
class LS2G_API PluginRegistrar {
public:
    void add(const std::string& name, AlgorithmRegistry::Factory f) {
        if (_staged.count(name)) {
            throw std::invalid_argument(
                "PluginRegistrar: this plugin registers the solver name '" + name + "' twice.");
        }
        _staged.emplace(name, std::move(f));
    }
    // Hand the staged solvers over to the caller (the entry-point helper),
    // leaving this registrar empty.
    AlgorithmRegistry::FactoryMap take() { return std::move(_staged); }

private:
    AlgorithmRegistry::FactoryMap _staged;
};

// Function pointer type every plugin's exported entry point matches. `registry`
// is an `AlgorithmRegistry*` passed as void* to keep the C boundary free of C++
// types; `errbuf`/`errbuf_len` receive a NUL-terminated diagnostic on failure.
// Returns 0 on success, non-zero on failure (never throws across this boundary).
extern "C" {
    typedef int (*ls2g_plugin_entry_fn)(void* registry, char* errbuf, std::size_t errbuf_len);
}

namespace detail {

// Body of the exported ls2g_register_plugin entry point (see LS2G_PLUGIN_ENTRY).
// `noexcept` and catching everything is the whole point: registration runs
// C++ code that can throw (ABI mismatch, duplicate name, bad_alloc, a bug in
// the plugin's own factory setup), and none of it may escape across the
// extern "C" boundary back into the C loader frames -- there it would be an
// uncatchable std::terminate(). Instead every failure becomes a return code
// plus a message the host turns into a normal Python exception.
inline int ls2g_run_plugin_entry(void* registry_ptr,
                                  char* errbuf,
                                  std::size_t errbuf_len,
                                  void (*register_fn)(PluginRegistrar&),
                                  Ls2gAbiTag plugin_tag) noexcept {
    // Bounded, overflow-proof copy of `msg` into `errbuf`: writes at most
    // errbuf_len-1 bytes and always NUL-terminates; a null or zero-length
    // buffer is simply left untouched. Written out by hand (rather than
    // strncpy) so the bound is obvious to a reader and to valgrind.
    auto write_err = [errbuf, errbuf_len](const char* msg) noexcept {
        if (errbuf == nullptr || errbuf_len == 0) return;
        std::size_t i = 0;
        const std::size_t last = errbuf_len - 1;
        if (msg != nullptr) {
            for (; i < last && msg[i] != '\0'; ++i) errbuf[i] = msg[i];
        }
        errbuf[i] = '\0';
    };
    try {
        if (registry_ptr == nullptr) {
            write_err("lightsim2grid plugin loader passed a null registry pointer.");
            return 2;
        }
        AlgorithmRegistry& reg = *static_cast<AlgorithmRegistry*>(registry_ptr);
        PluginRegistrar staging;
        register_fn(staging);                       // plugin stages its solvers
        reg.register_all(staging.take(), plugin_tag);  // atomic commit + ABI check
        return 0;
    } catch (const std::exception& e) {
        write_err(e.what());
        return 1;
    } catch (...) {
        write_err("unknown C++ exception during lightsim2grid plugin registration.");
        return 1;
    }
}

} // namespace detail

} // namespace ls2g

// Declare the exported entry point of a solver plugin. `register_fn` is a
// `void(ls2g::PluginRegistrar&)` the plugin author writes, calling reg.add(...)
// once per solver it provides. The generated ls2g_register_plugin function is
// what lightsim2grid.load_algorithm_plugin() looks up (via dlsym/GetProcAddress)
// and calls *after* the library is loaded -- so registration no longer happens
// inside a static constructor during dlopen(), and any failure comes back as a
// catchable Python exception instead of aborting the interpreter.
//
//   static void register_plugin(ls2g::PluginRegistrar& reg) {
//       reg.add("MyAlgo", []{ return std::make_unique<MyAlgo>(); });
//   }
//   LS2G_PLUGIN_ENTRY(register_plugin)
//
// ls2g_current_abi_tag() is evaluated here, in the plugin's own translation
// unit, so it captures the plugin's Eigen SIMD/alignment flags automatically.
#define LS2G_PLUGIN_ENTRY(register_fn)                                            \
    extern "C" LS2G_PLUGIN_EXPORT int                                             \
    ls2g_register_plugin(void* registry_ptr, char* errbuf, std::size_t errbuf_len) { \
        return ::ls2g::detail::ls2g_run_plugin_entry(                             \
            registry_ptr, errbuf, errbuf_len, &register_fn,                       \
            ::ls2g::ls2g_current_abi_tag());                                      \
    }

#endif // ALGORITHM_REGISTRY_H
