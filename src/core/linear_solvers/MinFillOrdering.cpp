// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "MinFillOrdering.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <queue>
#include <tuple>

namespace ls2g {

namespace {

// a vertex's place in the queue: (deficiency, degree, index), smallest first; `version` tells a
// stale entry (the vertex was re-scored since) from the current one
struct Entry {
    long deficiency;
    int degree;
    int vertex;
    unsigned version;
    bool operator>(const Entry & other) const {
        return std::tie(deficiency, degree, vertex) > std::tie(other.deficiency, other.degree, other.vertex);
    }
};

}  // namespace

bool min_fill_order(int n, const int * Ap, const int * Ai, std::vector<int> & perm){
    perm.clear();
    if(n < 0) return false;
    perm.reserve(n);
    if(n == 0) return true;

    // the graph of A + A^T without the diagonal, each adjacency sorted
    std::vector<std::vector<int> > adj(n);
    for(int j = 0; j < n; ++j){
        for(int k = Ap[j]; k < Ap[j + 1]; ++k){
            const int i = Ai[k];
            if(i < 0 || i >= n) return false;
            if(i == j) continue;
            adj[i].push_back(j);
            adj[j].push_back(i);
        }
    }
    for(std::vector<int> & a : adj){
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
    }

    // `mark[x] == stamp` flags a set of vertices without clearing anything between two uses
    std::vector<unsigned> mark(n, 0);
    unsigned stamp = 0;
    auto next_stamp = [&](){
        if(++stamp == 0){  // wrapped around: start again from a clean array
            std::fill(mark.begin(), mark.end(), 0);
            stamp = 1;
        }
        return stamp;
    };

    // missing edges among the neighbours of v: d (d - 1) / 2 pairs, minus the edges present,
    // each counted from both of its ends
    auto deficiency = [&](int v){
        const std::vector<int> & nv = adj[v];
        const unsigned s = next_stamp();
        for(int x : nv) mark[x] = s;
        long present = 0;
        for(int x : nv){
            for(int y : adj[x]) present += (mark[y] == s);
        }
        const long d = static_cast<long>(nv.size());
        return (d * (d - 1) - present) / 2;
    };

    std::vector<unsigned> version(n, 0);
    std::vector<char> alive(n, 1);
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry> > queue;
    auto push = [&](int v){
        ++version[v];
        queue.push(Entry{deficiency(v), static_cast<int>(adj[v].size()), v, version[v]});
    };
    for(int v = 0; v < n; ++v) push(v);

    std::vector<int> merged;
    std::vector<int> rescore;
    while(!queue.empty()){
        const Entry e = queue.top();
        queue.pop();
        const int v = e.vertex;
        if(!alive[v] || e.version != version[v]) continue;
        alive[v] = 0;
        perm.push_back(v);

        // eliminate v: its neighbours become a clique, v leaves their adjacency
        const std::vector<int> nv = std::move(adj[v]);
        adj[v].clear();
        for(int u : nv){
            std::vector<int> & au = adj[u];
            merged.clear();
            merged.reserve(au.size() + nv.size());
            std::set_union(au.begin(), au.end(), nv.begin(), nv.end(), std::back_inserter(merged));
            merged.erase(std::remove_if(merged.begin(), merged.end(),
                                        [v, u](int x){ return x == v || x == u; }),
                         merged.end());
            au.swap(merged);
        }

        // only the neighbours of v and their own neighbours can have changed scores
        const unsigned s = next_stamp();
        rescore.clear();
        for(int u : nv){
            if(mark[u] != s){ mark[u] = s; rescore.push_back(u); }
            for(int w : adj[u]){
                if(mark[w] != s){ mark[w] = s; rescore.push_back(w); }
            }
        }
        // the deficiency computation reuses `mark`: the set is already copied out
        for(int u : rescore){
            if(alive[u]) push(u);
        }
    }
    return static_cast<int>(perm.size()) == n;
}

} // namespace ls2g
