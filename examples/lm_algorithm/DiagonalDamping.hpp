// Copyright (c) 2026
// Experimental plugin, not part of the lightsim2grid core distribution.
//
// DiagonalEntry: a minimal NRSystem extension that reserves one Jacobian
// feature slot per existing unknown column (i, i), so a diagonal load can
// later be added to J without changing its sparsity pattern -- the standard
// trick behind Levenberg-Marquardt / Tikhonov-regularized Newton: solving
// (J + lambda*D) dx = -F instead of J dx = -F to regularize an ill-conditioned
// or near-singular Newton direction (as opposed to Iwamoto/MaxVoltageChange,
// which only ever rescale the *same*, already-bad, direction).
//
// DiagonalEntry claims no rows/cols of its own (register_in is a no-op
// besides remembering the ledger) and its fill_feature_values is a
// deliberate no-op: the damping value is never written through the normal,
// once-per-outer-NR-iteration component protocol path. It is written
// directly, per LM trial, via NRSystem::update_trailing_feature_values (the
// small, generic core addition this example relies on) -- see LMNRAlgo.hpp
// for how. This keeps every LM-specific concept (lambda, accept/reject, the
// Levenberg-vs-Marquardt-scaled delta choice) out of core entirely; core
// only knows "some extension reserved n trailing feature slots that can be
// updated fast."
//
// Must be the LAST extension in the NRSystem<Base, Rest...> pack: its n
// reserved feature handles are exactly the trailing block of NRSystem's
// internal `sink_`/`feature_pos_`, which is what lets
// update_trailing_feature_values locate them with nothing more than a count
// (see with_diagonal_entry below, which enforces this by construction).

#ifndef LM_DIAGONAL_DAMPING_H
#define LM_DIAGONAL_DAMPING_H

#include <powerflow_algorithm/NRSystem.hpp>

namespace ls2g {

class DiagonalEntry
{
public:
    void update_state(
        const Base                       * /*nr_system_base_ptr*/,
        const LSGrid                     * /*lsgrid_ptr*/,
        const EigenRefConstCplxSpMat     & /*Ybus*/,
        const Eigen::Ref<const CplxVect> & /*Sbus*/,
        const Eigen::Ref<const RealVect> & /*slack_weights*/) {}

    void init_topology(
        const Eigen::Ref<const IntVect>  & /*slack_ids*/,
        const Eigen::Ref<const RealVect> & /*slack_weights*/,
        const Eigen::Ref<const IntVect>  & /*pv*/,
        const Eigen::Ref<const IntVect>  & /*pq*/) {}

    // Claims no rows/cols of its own -- just remembers the ledger so
    // declare_feature_entries (called strictly later, from build_J_sparsity,
    // after every component's register_in has run) can read its final size.
    void register_in(NRLedger& ledger) { ledger_ptr_ = &ledger; }

    // One feature slot per EXISTING unknown column (i, i). A pure numerical
    // regularization term, not tied to any bus/physical quantity -- it needs
    // nothing beyond the final unknown count.
    void declare_feature_entries(FeatureSink& sink) {
        const int n = static_cast<int>(ledger_ptr_->size());
        for (int i = 0; i < n; ++i) sink.add(i, i);
    }

    // No-op on purpose: see the file header comment. Damping values are
    // written directly by NRSystem::update_trailing_feature_values, not
    // through this once-per-outer-iteration path.
    void fill_feature_values(FeatureWriter& /*writer*/, const Eigen::Ref<const RealVect>& /*Va*/) const {}

    void adjust_mismatch(const Eigen::Ref<const CplxVect>& /*V_t*/,
                         const Eigen::Ref<const RealVect>& /*dx*/,
                         Eigen::Ref<CplxVect> /*mis*/) const {}
    void fill_custom_rows(Eigen::Ref<RealVect> /*res*/,
                          const Eigen::Ref<const RealVect>& /*Va*/,
                          const Eigen::Ref<const RealVect>& /*Vm*/,
                          const Eigen::Ref<const RealVect>& /*dx*/) const {}
    void apply_step(const Eigen::Ref<const RealVect>& /*dx*/) {}

    void clear() { ledger_ptr_ = nullptr; }

private:
    const NRLedger* ledger_ptr_ = nullptr;
};

// ---- Injects DiagonalEntry as the LAST extension of an existing NRSystem<Base, Rest...> ----
template <class T> struct with_diagonal_entry;
template <class... Rest>
struct with_diagonal_entry<NRSystem<Base, Rest...>> {
    using type = NRSystem<Base, Rest..., DiagonalEntry>;
};

} // namespace ls2g

#endif // LM_DIAGONAL_DAMPING_H
