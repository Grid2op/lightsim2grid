// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LS2G_BUSGRAPH_H
#define LS2G_BUSGRAPH_H

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "Utils.hpp"

namespace ls2g {

namespace bus_graph_detail {
// a Coeff carries a complex value; applied to a real matrix (the DC Bbus) it is its
// real part that is subtracted -- the same rule as YbusPolicy::Contingency::readd_to_Ybus
template<typename T> struct CoeffAs { static T get(const cplx_type & v){ return v; } };
template<> struct CoeffAs<real_type> { static real_type get(const cplx_type & v){ return std::real(v); } };
}  // namespace bus_graph_detail

/**
 * The connectivity of a bus graph, prepared once so that "what does removing these
 * few edges do to it" is answered per query without walking the graph again.
 *
 * The graph is the sparsity pattern of an admittance matrix (Ybus or Bbus): an
 * edge between buses i and j wherever |M(i,j)| exceeds a threshold, the diagonal
 * ignored. That is exactly the graph the breadth-first searches of the batch
 * algorithms walk, so a verdict from here and a verdict from one of those agree
 * bus for bus -- and the fallback for the cases this cannot settle is one of them.
 *
 * One depth-first search from bus 0 gives everything a contingency needs:
 *   * the DFS tree, with each vertex's preorder rank and subtree size, so the
 *     subtree under a vertex is the contiguous rank range [pre, pre + size) and
 *     "is this bus below that vertex" is two integer compares;
 *   * the bridges (Tarjan's low-link), i.e. the tree edges whose removal splits the
 *     graph -- and then the two sides are the subtree of the child end and the rest.
 *
 * So for a contingency that removes ONE edge: a non-tree edge changes nothing, a
 * tree edge that is not a bridge changes nothing, a bridge strands the smaller
 * side, listed in time proportional to its size. A contingency that removes several
 * edges is settled only when none of them is a tree edge (the tree survives, so the
 * graph stays connected); otherwise it is reported as unknown and the caller falls
 * back to a search. That is the N-k case, and it is rare.
 *
 * Ties -- a bridge splitting the graph in two equal halves -- keep the side holding
 * bus 0, which is what a component labelling that starts from bus 0 and keeps the
 * first largest component does.
 *
 * Built from the solver-side matrix, so every bus id here is a SolverBusId's value.
 * The graph is read-only once built: several threads may query it at once.
 */
class BusGraph
{
    public:
        enum class Verdict {
            Connected,   // the graph stays in one piece
            Split,       // exactly one bridge went: `child` is the DFS child end of it
            Unknown      // cannot be settled from the tree alone: run a search
        };
        struct Cut {
            Verdict verdict;
            int child;   // meaningful for Split only
        };

        BusGraph() = default;

        /** the magnitude below which an admittance entry is no edge: the one the
            batch algorithms' searches have always used **/
        static real_type default_threshold() { return 1e-8; }

        /** Builds the graph of `mat`'s pattern (entries with |value| > threshold, off
            the diagonal) and runs the DFS. `mat` must be square and compressed, with
            a symmetric pattern (an admittance matrix is). **/
        template<typename T>
        void build(const Eigen::SparseMatrix<T> & mat, real_type threshold = 1e-8)
        {
            n_ = static_cast<int>(mat.cols());
            adj_start_.assign(static_cast<size_t>(n_) + 1, 0);
            adj_.clear();
            adj_.reserve(static_cast<size_t>(mat.nonZeros()));
            for(int col = 0; col < n_; ++col){
                for(typename Eigen::SparseMatrix<T>::InnerIterator it(mat, col); it; ++it){
                    const int row = static_cast<int>(it.row());
                    if(row == col) continue;
                    if(!(std::abs(it.value()) > threshold)) continue;
                    adj_.push_back(row);
                }
                adj_start_[static_cast<size_t>(col) + 1] = static_cast<int>(adj_.size());
            }
            _dfs();
        }

        bool built() const { return n_ >= 0; }
        int nb_bus() const { return n_; }
        /** whether the whole graph was reachable from bus 0 -- if not, every cut is
            Unknown: a labelling that starts disconnected is not a tree. **/
        bool base_connected() const { return connected_; }

        /**
         * The edges `edits` remove, for a set of coefficient edits applied to `mat`
         * (`value` is subtracted from entry (row_id, col_id), as
         * YbusPolicy::Contingency::remove_from_Ybus does). An edge is removed when
         * BOTH its entries fall to the threshold or below -- the pattern stays
         * symmetric. Returns false when an edit removes one direction and not the
         * other, which the tree cannot represent: report Unknown then.
         *
         * `out` receives each removed edge once, as (min bus, max bus).
         */
        template<typename T>
        static bool removed_edges(const Eigen::SparseMatrix<T> & mat,
                                  const std::vector<Coeff> & edits,
                                  real_type threshold,
                                  std::vector<std::pair<int, int> > & out)
        {
            out.clear();
            // the distinct off-diagonal entries edited (two parallel branches in one
            // contingency edit the same entry twice)
            std::vector<Coeff> total;
            for(const Coeff & c : edits){
                if(c.row_id == c.col_id) continue;
                bool found = false;
                for(const Coeff & t : total){
                    if(t.row_id == c.row_id && t.col_id == c.col_id){ found = true; break; }
                }
                if(!found) total.push_back(c);
            }
            // per directed entry, whether it goes below the threshold -- the edits
            // subtracted one by one, in their order, exactly as the edit itself does,
            // so the rounding is the same and the verdict agrees with a search of the
            // edited matrix
            std::vector<char> gone(total.size(), 0);
            for(size_t k = 0; k < total.size(); ++k){
                T after = mat.coeff(total[k].row_id, total[k].col_id);
                for(const Coeff & c : edits){
                    if(c.row_id == total[k].row_id && c.col_id == total[k].col_id){
                        after -= bus_graph_detail::CoeffAs<T>::get(c.value);
                    }
                }
                gone[k] = std::abs(after) > threshold ? 0 : 1;
            }
            // pair each direction with its mirror
            std::vector<char> paired(total.size(), 0);
            for(size_t k = 0; k < total.size(); ++k){
                if(paired[k]) continue;
                size_t mirror = total.size();
                for(size_t m = k + 1; m < total.size(); ++m){
                    if(total[m].row_id == total[k].col_id && total[m].col_id == total[k].row_id){ mirror = m; break; }
                }
                if(mirror == total.size()) return false;        // only one direction edited
                paired[k] = 1; paired[mirror] = 1;
                if(gone[k] != gone[mirror]) return false;        // asymmetric removal
                if(gone[k]){
                    const int a = static_cast<int>(total[k].row_id);
                    const int b = static_cast<int>(total[k].col_id);
                    out.emplace_back(std::min(a, b), std::max(a, b));
                }
            }
            return true;
        }

        /** The effect of removing `edges` (each as (i, j), i != j, present in the
            graph or not -- an edge that never existed changes nothing). **/
        Cut cut(const std::vector<std::pair<int, int> > & edges) const
        {
            Cut res;
            res.child = -1;
            if(!connected_){ res.verdict = Verdict::Unknown; return res; }
            // which of them are tree edges, and of those, bridges
            int nb_tree = 0;
            int bridge_child = -1;
            for(const auto & e : edges){
                const int child = _tree_child(e.first, e.second);
                if(child < 0) continue;              // not a tree edge: the tree survives it
                ++nb_tree;
                if(bridge_child_[static_cast<size_t>(child)]) bridge_child = child;
            }
            if(nb_tree == 0){ res.verdict = Verdict::Connected; return res; }
            if(edges.size() == 1){
                if(bridge_child < 0){ res.verdict = Verdict::Connected; return res; }
                res.verdict = Verdict::Split;
                res.child = bridge_child;
                return res;
            }
            // several edges, at least one of them a tree edge: whether the back edges
            // that covered it survive is not known here
            res.verdict = Verdict::Unknown;
            return res;
        }

        /** Split only: whether `bus` is on the side that is stranded (the smaller
            one; on a tie, the one without bus 0). **/
        bool is_stranded(const Cut & cut, int bus) const
        {
            return _in_subtree(cut.child, bus) == _subtree_is_stranded(cut.child);
        }

        /** Split only: appends the stranded side to `out`, in increasing bus order. **/
        void stranded_buses(const Cut & cut, std::vector<int> & out) const
        {
            const int c = cut.child;
            const int lo = pre_[static_cast<size_t>(c)];
            const int hi = lo + size_[static_cast<size_t>(c)];
            const size_t first = out.size();
            if(_subtree_is_stranded(c)){
                for(int k = lo; k < hi; ++k) out.push_back(order_[static_cast<size_t>(k)]);
            } else {
                for(int k = 0; k < lo; ++k) out.push_back(order_[static_cast<size_t>(k)]);
                for(int k = hi; k < n_; ++k) out.push_back(order_[static_cast<size_t>(k)]);
            }
            std::sort(out.begin() + static_cast<std::ptrdiff_t>(first), out.end());
        }

    private:
        // the DFS child end of tree edge (a, b), or -1 if (a, b) is not a tree edge
        int _tree_child(int a, int b) const
        {
            if(a < 0 || b < 0 || a >= n_ || b >= n_) return -1;
            if(parent_[static_cast<size_t>(b)] == a) return b;
            if(parent_[static_cast<size_t>(a)] == b) return a;
            return -1;
        }
        bool _in_subtree(int child, int bus) const
        {
            const int p = pre_[static_cast<size_t>(bus)];
            const int lo = pre_[static_cast<size_t>(child)];
            return p >= lo && p < lo + size_[static_cast<size_t>(child)];
        }
        // the subtree is the stranded side when it is the smaller one, or on a tie
        // (bus 0 is the root: it is never inside a proper subtree)
        bool _subtree_is_stranded(int child) const
        {
            const int sz = size_[static_cast<size_t>(child)];
            return sz <= n_ - sz;
        }

        // iterative DFS from bus 0: preorder ranks, subtree sizes, parents, and
        // Tarjan's low-link for the bridges. Iterative rather than recursive so a
        // long radial feeder cannot exhaust the stack.
        void _dfs()
        {
            const size_t n = static_cast<size_t>(n_);
            pre_.assign(n, -1);
            low_.assign(n, 0);
            parent_.assign(n, -1);
            size_.assign(n, 0);
            bridge_child_.assign(n, 0);
            order_.clear();
            order_.reserve(n);
            connected_ = true;
            if(n_ == 0) return;

            std::vector<std::pair<int, int> > stack;   // (vertex, next neighbour slot)
            stack.reserve(n);
            int counter = 0;
            pre_[0] = counter++;
            low_[0] = pre_[0];
            order_.push_back(0);
            stack.emplace_back(0, adj_start_[0]);
            while(!stack.empty()){
                const int v = stack.back().first;
                int & slot = stack.back().second;
                if(slot < adj_start_[static_cast<size_t>(v) + 1]){
                    const int w = adj_[static_cast<size_t>(slot)];
                    ++slot;
                    if(w == parent_[static_cast<size_t>(v)]) continue;   // no multi-edges in a matrix pattern
                    if(pre_[static_cast<size_t>(w)] < 0){
                        parent_[static_cast<size_t>(w)] = v;
                        pre_[static_cast<size_t>(w)] = counter++;
                        low_[static_cast<size_t>(w)] = pre_[static_cast<size_t>(w)];
                        order_.push_back(w);
                        stack.emplace_back(w, adj_start_[static_cast<size_t>(w)]);
                    } else {
                        low_[static_cast<size_t>(v)] = std::min(low_[static_cast<size_t>(v)], pre_[static_cast<size_t>(w)]);
                    }
                } else {
                    stack.pop_back();
                    size_[static_cast<size_t>(v)] = counter - pre_[static_cast<size_t>(v)];
                    const int u = parent_[static_cast<size_t>(v)];
                    if(u >= 0){
                        low_[static_cast<size_t>(u)] = std::min(low_[static_cast<size_t>(u)], low_[static_cast<size_t>(v)]);
                        if(low_[static_cast<size_t>(v)] > pre_[static_cast<size_t>(u)]) bridge_child_[static_cast<size_t>(v)] = 1;
                    }
                }
            }
            connected_ = (counter == n_);
        }

        int n_ = -1;
        bool connected_ = false;
        // CSR adjacency
        std::vector<int> adj_start_;
        std::vector<int> adj_;
        // the DFS tree
        std::vector<int> pre_;            // preorder rank of each bus
        std::vector<int> low_;            // Tarjan's low-link
        std::vector<int> parent_;         // -1 for the root
        std::vector<int> size_;           // subtree size
        std::vector<int> order_;          // bus at each preorder rank
        std::vector<char> bridge_child_;  // 1 iff the edge to the parent is a bridge
};

}  // namespace ls2g

#endif  // LS2G_BUSGRAPH_H
