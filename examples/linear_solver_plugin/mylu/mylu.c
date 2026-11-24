// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

/* see mylu.h: a dense LU with partial pivoting behind a sparse-solver-like C API */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mylu.h"

struct mylu_handle {
    int n;
    int * colptr;
    int * rowind;
    int * perm;     /* row i of the factors is row perm[i] of the matrix */
    double * lu;    /* n x n, row major, L (unit diagonal) and U in place */
};

int mylu_create(mylu_handle ** h){
    *h = (mylu_handle *) calloc(1, sizeof(mylu_handle));
    return *h != NULL ? 0 : -1;
}

void mylu_destroy(mylu_handle * h){
    if(h == NULL) return;
    free(h->colptr); free(h->rowind); free(h->perm); free(h->lu);
    free(h);
}

int mylu_analyze(mylu_handle * h, int n, const int * colptr, const int * rowind){
    const int nnz = colptr[n];
    free(h->colptr); free(h->rowind); free(h->perm); free(h->lu);
    h->n = n;
    h->colptr = (int *) malloc(sizeof(int) * (size_t)(n + 1));
    h->rowind = (int *) malloc(sizeof(int) * (size_t)(nnz > 0 ? nnz : 1));
    h->perm = (int *) malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    h->lu = (double *) malloc(sizeof(double) * (size_t)(n > 0 ? n * n : 1));
    if(!h->colptr || !h->rowind || !h->perm || !h->lu) return -1;
    memcpy(h->colptr, colptr, sizeof(int) * (size_t)(n + 1));
    memcpy(h->rowind, rowind, sizeof(int) * (size_t)nnz);
    for(int i = 0; i < n; ++i) h->perm[i] = i;
    return 0;
}

/* writes the matrix into lu, its rows in the order of perm */
static void scatter(mylu_handle * h, const double * values){
    const int n = h->n;
    int * where = (int *) malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    for(int i = 0; i < n; ++i) where[h->perm[i]] = i;
    memset(h->lu, 0, sizeof(double) * (size_t)(n * n));
    for(int j = 0; j < n; ++j){
        for(int k = h->colptr[j]; k < h->colptr[j + 1]; ++k){
            h->lu[where[h->rowind[k]] * n + j] = values[k];
        }
    }
    free(where);
}

static int eliminate(mylu_handle * h, int choose_pivots){
    const int n = h->n;
    double * a = h->lu;
    for(int k = 0; k < n; ++k){
        if(choose_pivots){
            int p = k;
            for(int i = k + 1; i < n; ++i) if(fabs(a[i * n + k]) > fabs(a[p * n + k])) p = i;
            if(p != k){
                for(int j = 0; j < n; ++j){ double t = a[k * n + j]; a[k * n + j] = a[p * n + j]; a[p * n + j] = t; }
                int t = h->perm[k]; h->perm[k] = h->perm[p]; h->perm[p] = t;
            }
        }
        if(a[k * n + k] == 0.0 || !isfinite(a[k * n + k])) return -2;
        for(int i = k + 1; i < n; ++i){
            const double f = a[i * n + k] /= a[k * n + k];
            if(f != 0.0) for(int j = k + 1; j < n; ++j) a[i * n + j] -= f * a[k * n + j];
        }
    }
    return 0;
}

int mylu_factor(mylu_handle * h, const double * values){
    if(h->lu == NULL) return -1;
    for(int i = 0; i < h->n; ++i) h->perm[i] = i;
    scatter(h, values);
    return eliminate(h, 1);
}

int mylu_refactor(mylu_handle * h, const double * values){
    if(h->lu == NULL) return -1;
    scatter(h, values);
    return eliminate(h, 0);
}

int mylu_solve(mylu_handle * h, double * b){
    const int n = h->n;
    const double * a = h->lu;
    double * y = (double *) malloc(sizeof(double) * (size_t)(n > 0 ? n : 1));
    if(y == NULL) return -1;
    for(int i = 0; i < n; ++i){
        y[i] = b[h->perm[i]];
        for(int j = 0; j < i; ++j) y[i] -= a[i * n + j] * y[j];
    }
    for(int i = n - 1; i >= 0; --i){
        for(int j = i + 1; j < n; ++j) y[i] -= a[i * n + j] * y[j];
        y[i] /= a[i * n + i];
    }
    memcpy(b, y, sizeof(double) * (size_t)n);
    free(y);
    return 0;
}
