// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef MIN_FILL_ORDERING_H
#define MIN_FILL_ORDERING_H

#include <vector>

#include "ls2g_api.hpp"

namespace ls2g {

/**
A fill-reducing symmetric ordering by minimum fill (minimum deficiency): the index eliminated
next is the one whose elimination creates the fewest new nonzeros, i.e. the one with the
fewest missing edges among its current neighbours in the elimination graph. Ties go to the
smallest current degree, then to the smallest index, so that the order is deterministic.

This is the third of the three ordering schemes of

- W. F. Tinney, J. W. Walker, "Direct solutions of sparse network equations by optimally
  ordered triangular factorization", Proceedings of the IEEE 55(11):1801-1809, 1967.

(scheme 1 orders by the degree in the original graph, scheme 2 is minimum degree). The
elimination graph is the one of

- D. J. Rose, "A graph-theoretic study of the numerical solution of sparse positive definite
  systems of linear equations", in Graph Theory and Computing, R. C. Read (ed.), Academic
  Press, 1972, pp. 183-217.

eliminating a vertex joins its neighbours into a clique. The deficiency of a vertex only
changes when an edge appears between two of its neighbours, so after an elimination only the
neighbours of the eliminated vertex and their own neighbours are re-scored; every score is
exact, none is an approximation.

The ordering is computed on the pattern of A + A^T without its diagonal: (Ap, Ai) is a
compressed-column pattern, which need not be symmetric. Returns false on an invalid pattern
(an index out of range); `perm` (size n) is the elimination order, perm[k] being the index
eliminated k-th.

It costs more than AMD (exact scores, explicit elimination graph): meant for an analysis done
once and amortised over many refactorizations.
**/
LS2G_API bool min_fill_order(int n, const int * Ap, const int * Ai, std::vector<int> & perm);

} // namespace ls2g

#endif // MIN_FILL_ORDERING_H
