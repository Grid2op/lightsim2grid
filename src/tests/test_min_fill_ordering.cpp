// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// min_fill_order: a permutation, deterministic, built on the pattern of A + A^T, that picks a
// vertex creating no fill whenever one exists -- so that a graph with a perfect elimination
// order (a chordal one) is ordered with no fill at all. C++14 only.

#include <algorithm>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "linear_solvers/MinFillOrdering.hpp"

namespace {

// a compressed-column pattern, from (row, col) entries, the diagonal always present
struct Pattern {
    int n;
    std::vector<int> Ap, Ai;
};

Pattern make_pattern(int n, const std::vector<std::pair<int, int> > & entries){
    std::vector<std::set<int> > cols(n);
    for(int j = 0; j < n; ++j) cols[j].insert(j);
    for(const auto & e : entries) cols[e.second].insert(e.first);
    Pattern p{n, {0}, {}};
    for(int j = 0; j < n; ++j){
        p.Ai.insert(p.Ai.end(), cols[j].begin(), cols[j].end());
        p.Ap.push_back(static_cast<int>(p.Ai.size()));
    }
    return p;
}

// both directions of every edge
std::vector<std::pair<int, int> > symmetric(const std::vector<std::pair<int, int> > & edges){
    std::vector<std::pair<int, int> > out;
    for(const auto & e : edges){
        out.push_back(e);
        out.push_back(std::make_pair(e.second, e.first));
    }
    return out;
}

// the fill (new edges) of eliminating the vertices in `perm` order, on the graph of A + A^T
long fill_of(const Pattern & p, const std::vector<int> & perm){
    std::vector<std::set<int> > adj(p.n);
    for(int j = 0; j < p.n; ++j){
        for(int k = p.Ap[j]; k < p.Ap[j + 1]; ++k){
            const int i = p.Ai[k];
            if(i != j){ adj[i].insert(j); adj[j].insert(i); }
        }
    }
    long fill = 0;
    for(int v : perm){
        const std::set<int> nv = adj[v];
        for(int u : nv){
            adj[u].erase(v);
            for(int w : nv){
                if(w != u && adj[u].insert(w).second && u < w) ++fill;
            }
        }
        adj[v].clear();
    }
    return fill;
}

bool is_permutation_of_n(const std::vector<int> & perm, int n){
    if(static_cast<int>(perm.size()) != n) return false;
    std::vector<int> sorted = perm;
    std::sort(sorted.begin(), sorted.end());
    for(int i = 0; i < n; ++i){
        if(sorted[i] != i) return false;
    }
    return true;
}

// a random 3-tree: each new vertex joined to the three vertices of an existing triangle,
// a chordal graph, then the vertices relabelled at random
std::vector<std::pair<int, int> > random_three_tree(int n, unsigned seed){
    std::mt19937 gen(seed);
    std::vector<std::vector<int> > triangles{{0, 1, 2}};
    std::vector<std::pair<int, int> > edges{{0, 1}, {0, 2}, {1, 2}};
    for(int v = 3; v < n; ++v){
        std::uniform_int_distribution<std::size_t> pick(0, triangles.size() - 1);
        const std::vector<int> t = triangles[pick(gen)];
        for(int u : t) edges.push_back(std::make_pair(u, v));
        triangles.push_back({t[0], t[1], v});
        triangles.push_back({t[0], t[2], v});
        triangles.push_back({t[1], t[2], v});
    }
    std::vector<int> label(n);
    for(int i = 0; i < n; ++i) label[i] = i;
    std::shuffle(label.begin(), label.end(), gen);
    for(auto & e : edges) e = std::make_pair(label[e.first], label[e.second]);
    return edges;
}

}  // namespace

TEST_CASE("min-fill orders an empty pattern and rejects an index out of range", "[linear_solver][ordering][min_fill]")
{
    std::vector<int> perm{1, 2, 3};
    REQUIRE(ls2g::min_fill_order(0, nullptr, nullptr, perm));
    CHECK(perm.empty());

    const int Ap[] = {0, 1, 3};
    const int Ai[] = {0, 1, 2};  // row 2 in a 2-by-2 pattern
    CHECK_FALSE(ls2g::min_fill_order(2, Ap, Ai, perm));
}

TEST_CASE("min-fill gives a deterministic permutation", "[linear_solver][ordering][min_fill]")
{
    const int n = 200;
    const Pattern p = make_pattern(n, symmetric(random_three_tree(n, 3)));
    std::vector<int> perm, again;
    REQUIRE(ls2g::min_fill_order(n, p.Ap.data(), p.Ai.data(), perm));
    REQUIRE(ls2g::min_fill_order(n, p.Ap.data(), p.Ai.data(), again));
    CHECK(is_permutation_of_n(perm, n));
    CHECK(perm == again);
}

TEST_CASE("min-fill orders a chordal graph with no fill", "[linear_solver][ordering][min_fill]")
{
    // a chordal graph always has a vertex whose neighbours form a clique, and eliminating it
    // leaves a chordal graph: picking a zero-fill vertex whenever one exists never fills in
    for(unsigned seed : {1u, 2u, 7u}){
        const int n = 300;
        const Pattern p = make_pattern(n, symmetric(random_three_tree(n, seed)));
        std::vector<int> perm;
        REQUIRE(ls2g::min_fill_order(n, p.Ap.data(), p.Ai.data(), perm));
        REQUIRE(is_permutation_of_n(perm, n));
        CHECK(fill_of(p, perm) == 0);
    }
}

TEST_CASE("min-fill eliminates a vertex creating no fill before a smaller-degree one that would", "[linear_solver][ordering][min_fill]")
{
    // vertex 0 has degree 2 but its neighbours 1 and 2 are not joined (one fill edge);
    // 3, 4, 5, 6 form a clique hanging off vertex 2 through 3: vertices 4, 5 and 6 have
    // degree 3 and no fill. Minimum degree would take 0 (or 1) first; minimum fill takes 4.
    const std::vector<std::pair<int, int> > edges{
        {0, 1}, {0, 2}, {1, 7}, {2, 7}, {2, 3},
        {3, 4}, {3, 5}, {3, 6}, {4, 5}, {4, 6}, {5, 6}};
    const Pattern p = make_pattern(8, symmetric(edges));
    std::vector<int> perm;
    REQUIRE(ls2g::min_fill_order(8, p.Ap.data(), p.Ai.data(), perm));
    REQUIRE(is_permutation_of_n(perm, 8));
    CHECK(perm[0] == 4);  // smallest index among the zero-fill vertices of smallest degree
    CHECK(fill_of(p, perm) == 1);  // the cycle 0-1-7-2 needs one chord, nothing else fills
}

TEST_CASE("min-fill orders the pattern of A + A^T", "[linear_solver][ordering][min_fill]")
{
    // every edge of a 3-tree stored once, in one triangle or the other at random: the order is
    // the one of the symmetric pattern
    const int n = 150;
    const std::vector<std::pair<int, int> > edges = random_three_tree(n, 11);
    std::mt19937 gen(5);
    std::vector<std::pair<int, int> > one_sided;
    for(const auto & e : edges){
        one_sided.push_back((gen() & 1u) ? e : std::make_pair(e.second, e.first));
    }
    const Pattern sym = make_pattern(n, symmetric(edges));
    const Pattern half = make_pattern(n, one_sided);
    std::vector<int> perm_sym, perm_half;
    REQUIRE(ls2g::min_fill_order(n, sym.Ap.data(), sym.Ai.data(), perm_sym));
    REQUIRE(ls2g::min_fill_order(n, half.Ap.data(), half.Ai.data(), perm_half));
    CHECK(perm_sym == perm_half);
}
