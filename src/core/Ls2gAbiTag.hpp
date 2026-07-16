// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LS2G_ABI_TAG_H
#define LS2G_ABI_TAG_H

#include <sstream>
#include <string>

#include "ls2g_api.hpp"
#include "Eigen/Core"

namespace ls2g {

// Fingerprint of the Eigen SIMD/alignment configuration a given translation
// unit was compiled with. -march=native (or any -m<isa> flag) changes which
// of EIGEN_VECTORIZE_* get defined, which changes EIGEN_MAX_ALIGN_BYTES,
// which changes how Eigen aligns/allocates/frees its dynamic-size matrices.
// If two separately-compiled binaries that share Eigen objects (core vs. the
// Python bindings, or core vs. a solver plugin) disagree on this, an object
// allocated under one alignment assumption and freed under another silently
// corrupts the heap -- this is a real, previously-hit bug, not a
// hypothetical one (see CHANGELOG.rst / docs/solver_plugin.rst).
//
// `ls2g_current_abi_tag()` below is `inline` on purpose: it must be
// recompiled in each translation unit that calls it (core, bindings, a
// plugin) so that it reflects *that* TU's own compiler flags, not a value
// baked once somewhere else. Comparing two independently-evaluated tags is
// the only way to observe this kind of drift from outside a single
// compilation -- a static_assert cannot do it, because no single TU has
// visibility into what flags another, separately-compiled TU used.
struct Ls2gAbiTag {
    // Bump only if this struct's meaning changes (e.g. a field is added or
    // reinterpreted), not for ordinary Eigen/compiler upgrades.
    unsigned tag_version = 1;

    unsigned eigen_max_align_bytes = EIGEN_MAX_ALIGN_BYTES;
    unsigned eigen_world_version = EIGEN_WORLD_VERSION;
    unsigned eigen_major_version = EIGEN_MAJOR_VERSION;
    unsigned eigen_minor_version = EIGEN_MINOR_VERSION;

#ifdef EIGEN_VECTORIZE_SSE
    bool vectorize_sse = true;
#else
    bool vectorize_sse = false;
#endif
#ifdef EIGEN_VECTORIZE_AVX
    bool vectorize_avx = true;
#else
    bool vectorize_avx = false;
#endif
#ifdef EIGEN_VECTORIZE_AVX2
    bool vectorize_avx2 = true;
#else
    bool vectorize_avx2 = false;
#endif
#ifdef EIGEN_VECTORIZE_AVX512
    bool vectorize_avx512 = true;
#else
    bool vectorize_avx512 = false;
#endif
#ifdef EIGEN_VECTORIZE_FMA
    bool vectorize_fma = true;
#else
    bool vectorize_fma = false;
#endif
#ifdef EIGEN_VECTORIZE_NEON
    bool vectorize_neon = true;
#else
    bool vectorize_neon = false;
#endif

    bool operator==(const Ls2gAbiTag& o) const {
        return tag_version == o.tag_version
            && eigen_max_align_bytes == o.eigen_max_align_bytes
            && eigen_world_version == o.eigen_world_version
            && eigen_major_version == o.eigen_major_version
            && eigen_minor_version == o.eigen_minor_version
            && vectorize_sse == o.vectorize_sse
            && vectorize_avx == o.vectorize_avx
            && vectorize_avx2 == o.vectorize_avx2
            && vectorize_avx512 == o.vectorize_avx512
            && vectorize_fma == o.vectorize_fma
            && vectorize_neon == o.vectorize_neon;
    }
    bool operator!=(const Ls2gAbiTag& o) const { return !(*this == o); }

    // Human-readable summary for error messages, e.g. "align=32 bytes,
    // Eigen 5.0.1, SSE+AVX+AVX2+FMA".
    std::string describe() const {
        std::ostringstream oss;
        oss << "align=" << eigen_max_align_bytes << " bytes, Eigen "
            << eigen_world_version << "." << eigen_major_version << "." << eigen_minor_version << ", ";
        bool first = true;
        auto add = [&](bool present, const char* name) {
            if (!present) return;
            if (!first) oss << "+";
            oss << name;
            first = false;
        };
        add(vectorize_sse, "SSE");
        add(vectorize_avx, "AVX");
        add(vectorize_avx2, "AVX2");
        add(vectorize_avx512, "AVX512");
        add(vectorize_fma, "FMA");
        add(vectorize_neon, "NEON");
        if (first) oss << "no vectorization";
        return oss.str();
    }
};

// Evaluated in the caller's own translation unit -- see the class comment.
inline Ls2gAbiTag ls2g_current_abi_tag() { return Ls2gAbiTag{}; }

// Non-inline: always executes as compiled into lightsim2grid_core, so it
// reports core's real build flags regardless of who calls it (bindings, a
// plugin, ...). Defined in Ls2gAbiTag.cpp.
LS2G_API Ls2gAbiTag core_abi_tag();

} // namespace ls2g

#endif // LS2G_ABI_TAG_H
