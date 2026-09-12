// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// BusGraph answers "which buses does removing these edges strand" from one DFS tree
// of the base graph, where the batch algorithms used to run a breadth-first search
// per contingency. The reference here IS that search (the labelling of
// BaseBatchSweep::_disconnected_buses, copied), run on the edited matrix: on every
// graph below, for every single edge and for random pairs of edges, the tree's
// answer must be the search's answer bus for bus -- or "unknown", which is
// allowed only where the tree genuinely cannot tell (several tree edges gone).

#include <algorithm>
#include <complex>
#include <queue>
#include <random>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "BusGraph.hpp"

using ls2g::BusGraph;
using ls2g::Coeff;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using Edge = std::pair<int, int>;
using CplxSp = Eigen::SparseMatrix<cplx_type>;

// an admittance-shaped matrix: for every branch (i, j) of admittance y, -y off the
// diagonal on both sides and +y on both diagonal entries. Parallel branches add up.
CplxSp make_matrix(int n, const std::vector<Edge> & branches, const std::vector<cplx_type> & y)
{
    std::vector<Eigen::Triplet<cplx_type> > trip;
    for(size_t k = 0; k < branches.size(); ++k){
        const int i = branches[k].first, j = branches[k].second;
        trip.emplace_back(i, j, -y[k]);
        trip.emplace_back(j, i, -y[k]);
        trip.emplace_back(i, i, y[k]);
        trip.emplace_back(j, j, y[k]);
    }
    CplxSp res(n, n);
    res.setFromTriplets(trip.begin(), trip.end());
    res.makeCompressed();
    return res;
}

// the four coefficients ContingencyAnalysis subtracts to disconnect branch k -- the
// same shape as YbusPolicy::Contingency::_coeffs_for_branch_ids builds
std::vector<Coeff> disconnect(const std::vector<Edge> & branches, const std::vector<cplx_type> & y, int k)
{
    const int i = branches[static_cast<size_t>(k)].first, j = branches[static_cast<size_t>(k)].second;
    const cplx_type yk = y[static_cast<size_t>(k)];
    return {Coeff{i, i, yk}, Coeff{j, j, yk}, Coeff{i, j, -yk}, Coeff{j, i, -yk}};
}

// BaseBatchSweep::_disconnected_buses, verbatim: the solver bus ids NOT in the
// largest connected component of the matrix's pattern (empty if connected)
std::vector<int> reference_stranded(const CplxSp & mat)
{
    const int n = static_cast<int>(mat.cols());
    std::vector<int> comp_of_bus(n, -1);
    int nb_comp = 0;
    std::queue<int> neighborhood;
    for(int start = 0; start < n; ++start){
        if(comp_of_bus[start] != -1) continue;
        comp_of_bus[start] = nb_comp;
        neighborhood.push(start);
        while(!neighborhood.empty()){
            const int col_id = neighborhood.front();
            neighborhood.pop();
            for(CplxSp::InnerIterator it(mat, col_id); it; ++it){
                const int row = static_cast<int>(it.row());
                if(comp_of_bus[row] == -1 && std::abs(it.value()) > 1e-8){
                    comp_of_bus[row] = nb_comp;
                    neighborhood.push(row);
                }
            }
        }
        ++nb_comp;
    }
    if(nb_comp <= 1) return std::vector<int>();
    std::vector<int> nb_bus_per_comp(nb_comp, 0);
    for(int bus_id = 0; bus_id < n; ++bus_id) nb_bus_per_comp[comp_of_bus[bus_id]] += 1;
    const int main_comp = static_cast<int>(std::distance(
        nb_bus_per_comp.begin(), std::max_element(nb_bus_per_comp.begin(), nb_bus_per_comp.end())));
    std::vector<int> masked;
    for(int bus_id = 0; bus_id < n; ++bus_id) if(comp_of_bus[bus_id] != main_comp) masked.push_back(bus_id);
    return masked;
}

// what the batch does with a contingency: edit, search, restore
std::vector<int> reference_after(CplxSp & mat, const std::vector<Coeff> & edits)
{
    for(const Coeff & c : edits) mat.coeffRef(c.row_id, c.col_id) -= c.value;
    std::vector<int> res = reference_stranded(mat);
    for(const Coeff & c : edits) mat.coeffRef(c.row_id, c.col_id) += c.value;
    return res;
}

// what BusGraph does with it: Split -> the stranded list, Connected -> empty,
// Unknown -> nullptr (the caller falls back to the search)
bool tree_after(const BusGraph & graph, const CplxSp & base, const std::vector<Coeff> & edits,
                std::vector<int> & out)
{
    std::vector<Edge> removed;
    if(!BusGraph::removed_edges(base, edits, BusGraph::default_threshold(), removed)) return false;
    const BusGraph::Cut cut = graph.cut(removed);
    out.clear();
    if(cut.verdict == BusGraph::Verdict::Unknown) return false;
    if(cut.verdict == BusGraph::Verdict::Split) graph.stranded_buses(cut, out);
    return true;
}

// every single-branch contingency of `branches` must be settled by the tree and
// agree with the search; the N-2 ones must agree whenever the tree settles them
void check_all(int n, const std::vector<Edge> & branches, const std::vector<cplx_type> & y,
               int nb_pairs, unsigned seed)
{
    CplxSp mat = make_matrix(n, branches, y);
    BusGraph graph;
    graph.build(mat);
    REQUIRE(graph.base_connected());

    std::vector<int> from_tree;
    int nb_split = 0;
    for(int k = 0; k < static_cast<int>(branches.size()); ++k){
        const std::vector<Coeff> edits = disconnect(branches, y, k);
        const std::vector<int> expected = reference_after(mat, edits);
        INFO("branch " << k << " (" << branches[static_cast<size_t>(k)].first << ", "
             << branches[static_cast<size_t>(k)].second << ")");
        REQUIRE(tree_after(graph, mat, edits, from_tree));
        CHECK(from_tree == expected);
        if(!expected.empty()) ++nb_split;
        // is_stranded must say the same thing bus by bus
        if(!expected.empty()){
            std::vector<Edge> removed;
            BusGraph::removed_edges(mat, edits, BusGraph::default_threshold(), removed);
            const BusGraph::Cut cut = graph.cut(removed);
            for(int bus = 0; bus < n; ++bus){
                const bool ref = std::binary_search(expected.begin(), expected.end(), bus);
                CHECK(graph.is_stranded(cut, bus) == ref);
            }
        }
    }
    INFO("the graph must have some bridges for the test to mean anything");
    CHECK(nb_split > 0);

    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> pick(0, static_cast<int>(branches.size()) - 1);
    int nb_settled = 0;
    for(int p = 0; p < nb_pairs; ++p){
        const int a = pick(rng), b = pick(rng);
        if(a == b) continue;
        std::vector<Coeff> edits = disconnect(branches, y, a);
        const std::vector<Coeff> more = disconnect(branches, y, b);
        edits.insert(edits.end(), more.begin(), more.end());
        const std::vector<int> expected = reference_after(mat, edits);
        if(tree_after(graph, mat, edits, from_tree)){
            INFO("branches " << a << " and " << b);
            CHECK(from_tree == expected);
            ++nb_settled;
        }
    }
    INFO("an N-2 removing no tree edge is settled without a search");
    CHECK(nb_settled > 0);
}

}  // namespace


TEST_CASE("BusGraph: a small grid with rings, feeders and a parallel pair", "[bus_graph]")
{
    // a ring 0-1-2-3-4-5-0, a feeder 3-6-7-8 (three bridges), a second ring 8-9-10-8
    // at its end, a parallel pair between 1 and 4 (two branches, one entry), and a
    // pendant bus 11 on 10
    const std::vector<Edge> branches = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0},
        {3, 6}, {6, 7}, {7, 8},
        {8, 9}, {9, 10}, {10, 8},
        {1, 4}, {1, 4},
        {10, 11}};
    std::vector<cplx_type> y;
    for(size_t k = 0; k < branches.size(); ++k) y.emplace_back(1. + 0.1 * static_cast<real_type>(k), -3.);
    check_all(12, branches, y, 200, 1u);
}

TEST_CASE("BusGraph: removing one of two parallel branches keeps the edge", "[bus_graph]")
{
    // 0-1 by two branches, 1-2 by one: only (1, 2) is a bridge
    const std::vector<Edge> branches = {{0, 1}, {0, 1}, {1, 2}};
    const std::vector<cplx_type> y = {{1., -2.}, {0.5, -1.}, {1., -4.}};
    CplxSp mat = make_matrix(3, branches, y);
    BusGraph graph;
    graph.build(mat);

    std::vector<Edge> removed;
    REQUIRE(BusGraph::removed_edges(mat, disconnect(branches, y, 0), BusGraph::default_threshold(), removed));
    CHECK(removed.empty());
    CHECK(graph.cut(removed).verdict == BusGraph::Verdict::Connected);

    REQUIRE(BusGraph::removed_edges(mat, disconnect(branches, y, 2), BusGraph::default_threshold(), removed));
    REQUIRE(removed.size() == 1);
    const BusGraph::Cut cut = graph.cut(removed);
    REQUIRE(cut.verdict == BusGraph::Verdict::Split);
    std::vector<int> stranded;
    graph.stranded_buses(cut, stranded);
    CHECK(stranded == std::vector<int>{2});

    // both parallel branches at once: the entry does go to zero
    std::vector<Coeff> both = disconnect(branches, y, 0);
    const std::vector<Coeff> second = disconnect(branches, y, 1);
    both.insert(both.end(), second.begin(), second.end());
    REQUIRE(BusGraph::removed_edges(mat, both, BusGraph::default_threshold(), removed));
    REQUIRE(removed.size() == 1);
    CHECK(removed[0] == Edge(0, 1));
    CHECK(graph.cut(removed).verdict == BusGraph::Verdict::Split);
}

TEST_CASE("BusGraph: ties keep the side of bus 0, and the stranded side may hold the root", "[bus_graph]")
{
    SECTION("a path split in two equal halves strands the half without bus 0"){
        const std::vector<Edge> branches = {{0, 1}, {1, 2}, {2, 3}};
        const std::vector<cplx_type> y(3, cplx_type(1., -1.));
        CplxSp mat = make_matrix(4, branches, y);
        BusGraph graph;
        graph.build(mat);
        std::vector<int> stranded;
        REQUIRE(tree_after(graph, mat, disconnect(branches, y, 1), stranded));
        CHECK(stranded == std::vector<int>{2, 3});
        CHECK(stranded == reference_after(mat, disconnect(branches, y, 1)));
    }
    SECTION("a long feeder off the root strands the root's small side"){
        // 0-7-8-0 is a small ring; 0-1-2-3-4-5-6 a longer feeder: opening (0, 1)
        // leaves {0, 7, 8} as the smaller piece, the side holding the DFS root
        const std::vector<Edge> branches = {{0, 7}, {7, 8}, {8, 0},
                                            {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}};
        const std::vector<cplx_type> y(branches.size(), cplx_type(0.5, -2.));
        CplxSp mat = make_matrix(9, branches, y);
        BusGraph graph;
        graph.build(mat);
        std::vector<int> stranded;
        REQUIRE(tree_after(graph, mat, disconnect(branches, y, 3), stranded));
        CHECK(stranded == std::vector<int>{0, 7, 8});
        CHECK(stranded == reference_after(mat, disconnect(branches, y, 3)));
    }
}

TEST_CASE("BusGraph: what the tree cannot settle is reported, not guessed", "[bus_graph]")
{
    SECTION("an edit that removes one direction only"){
        const std::vector<Edge> branches = {{0, 1}, {1, 2}};
        const std::vector<cplx_type> y(2, cplx_type(1., -1.));
        CplxSp mat = make_matrix(3, branches, y);
        std::vector<Edge> removed;
        const std::vector<Coeff> lopsided = {Coeff{0, 1, -y[0]}, Coeff{1, 0, cplx_type(0., 0.)}};
        CHECK_FALSE(BusGraph::removed_edges(mat, lopsided, BusGraph::default_threshold(), removed));
    }
    SECTION("a base graph that is not connected"){
        const std::vector<Edge> branches = {{0, 1}, {2, 3}};
        const std::vector<cplx_type> y(2, cplx_type(1., -1.));
        CplxSp mat = make_matrix(4, branches, y);
        BusGraph graph;
        graph.build(mat);
        CHECK_FALSE(graph.base_connected());
        CHECK(graph.cut(std::vector<Edge>{}).verdict == BusGraph::Verdict::Unknown);
    }
    SECTION("two tree edges of one cycle: together they split, so the tree says unknown"){
        const std::vector<Edge> branches = {{0, 1}, {1, 2}, {2, 0}};
        const std::vector<cplx_type> y(3, cplx_type(1., -1.));
        CplxSp mat = make_matrix(3, branches, y);
        BusGraph graph;
        graph.build(mat);
        int nb_unknown = 0;
        for(int a = 0; a < 3; ++a){
            for(int b = a + 1; b < 3; ++b){
                std::vector<Coeff> edits = disconnect(branches, y, a);
                const std::vector<Coeff> more = disconnect(branches, y, b);
                edits.insert(edits.end(), more.begin(), more.end());
                std::vector<int> stranded;
                if(!tree_after(graph, mat, edits, stranded)) ++nb_unknown;
                else CHECK(stranded == reference_after(mat, edits));
            }
        }
        // a triangle's DFS tree has two tree edges and one back edge: the pair made
        // of the two tree edges cannot be settled, the two pairs holding the back
        // edge and one tree edge cannot either
        CHECK(nb_unknown == 3);
    }
}

TEST_CASE("BusGraph: random graphs against the search", "[bus_graph]")
{
    for(unsigned seed = 1; seed <= 6; ++seed){
        std::mt19937 rng(seed);
        const int n = 40;
        std::vector<Edge> branches;
        // a spanning path so the base graph is connected, then a few chords, then a
        // few pendant feeders (bridges), then a couple of parallel branches
        for(int i = 1; i < 30; ++i) branches.emplace_back(i - 1, i);
        std::uniform_int_distribution<int> pick(0, 29);
        for(int c = 0; c < 12; ++c){
            const int a = pick(rng), b = pick(rng);
            if(a != b) branches.emplace_back(std::min(a, b), std::max(a, b));
        }
        for(int i = 30; i < n; ++i) branches.emplace_back(pick(rng), i);
        branches.push_back(branches[3]);
        branches.push_back(branches[7]);
        std::vector<cplx_type> y;
        std::uniform_real_distribution<real_type> mag(0.5, 5.);
        for(size_t k = 0; k < branches.size(); ++k) y.emplace_back(mag(rng), -mag(rng));
        INFO("seed " << seed);
        check_all(n, branches, y, 300, seed);
    }
}
