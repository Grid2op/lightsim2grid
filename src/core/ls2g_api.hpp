// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LS2G_API_H
#define LS2G_API_H

// Symbols that must be visible across .so/.dll boundaries carry LS2G_API.
// On Linux/macOS: visibility("default") overrides -fvisibility=hidden.
// On Windows/MSVC: dllexport when building the core library (LS2G_BUILDING_CORE
//   defined by CMake), dllimport when included by a plugin or the Python wrapper.
#ifndef LS2G_API
#  if defined(_MSC_VER)
#    ifdef LS2G_BUILDING_CORE
#      define LS2G_API __declspec(dllexport)
#    else
#      define LS2G_API __declspec(dllimport)
#    endif
#  elif defined(__GNUC__) || defined(__clang__)
#    define LS2G_API __attribute__((visibility("default")))
#  else
#    define LS2G_API
#  endif
#endif

// A solver plugin (a separate .so/.dll loaded at runtime via
// lightsim2grid.load_algorithm_plugin()) must *export* its registration entry
// point, regardless of which side it is built on. LS2G_API resolves to
// dllimport inside a plugin (LS2G_BUILDING_CORE is not defined there), so it is
// the wrong macro for a symbol the plugin itself provides -- hence this
// dedicated always-export macro. See LS2G_PLUGIN_ENTRY in AlgorithmRegistry.hpp.
#ifndef LS2G_PLUGIN_EXPORT
#  if defined(_MSC_VER)
#    define LS2G_PLUGIN_EXPORT __declspec(dllexport)
#  elif defined(__GNUC__) || defined(__clang__)
#    define LS2G_PLUGIN_EXPORT __attribute__((visibility("default")))
#  else
#    define LS2G_PLUGIN_EXPORT
#  endif
#endif

#endif // LS2G_API_H
