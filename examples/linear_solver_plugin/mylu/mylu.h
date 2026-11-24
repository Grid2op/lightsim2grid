// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

/* mylu: a stand-in for a third-party sparse direct solver with a C API, so that this example
   builds and runs on its own. Underneath it is a dense LU with partial pivoting: fine for the
   small grids of test_plugin.py, far too slow for anything else. A real library replaces this
   folder; the wrapper (MyLULinearSolver.hpp) is written against this header only. */
#ifndef MYLU_H
#define MYLU_H
#ifdef __cplusplus
extern "C" {
#endif

typedef struct mylu_handle mylu_handle;

/* every function returns 0 on success, a negative value on failure */
int  mylu_create(mylu_handle ** h);
void mylu_destroy(mylu_handle * h);
/* the sparsity pattern: compressed sparse column, 0-based, n + 1 column pointers */
int  mylu_analyze(mylu_handle * h, int n, const int * colptr, const int * rowind);
/* the values, in the order of rowind: chooses the pivots */
int  mylu_factor(mylu_handle * h, const double * values);
/* same pattern, new values: reuses the pivots of the last mylu_factor; fails on a zero pivot */
int  mylu_refactor(mylu_handle * h, const double * values);
/* solves A x = b in place */
int  mylu_solve(mylu_handle * h, double * b);

#ifdef __cplusplus
}
#endif
#endif
