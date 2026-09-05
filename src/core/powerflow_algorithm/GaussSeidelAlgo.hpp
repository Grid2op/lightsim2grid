// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GAUSSSEIDEL_ALGO_H
#define GAUSSSEIDEL_ALGO_H

#include "BaseAlgo.hpp"

namespace ls2g {

class LS2G_API GaussSeidelAlgo : public BaseAlgo
{
    public:
        GaussSeidelAlgo() noexcept :BaseAlgo(true) {};

        ~GaussSeidelAlgo() noexcept override = default;

        // compute_pf leaves BaseAlgo::mis_bus_ filled at the voltage it returns --
        // with the plain formula, this family has no extension state to fold in.
        // Inherited by GaussSeidelSynchAlgo, which shares compute_pf.
        static constexpr bool FILLS_BUS_MISMATCH = true;
        bool fills_bus_mismatch() const noexcept override { return FILLS_BUS_MISMATCH; }

        // todo  can be factorized
        Eigen::Ref<const Eigen::SparseMatrix<real_type> > get_J() const override {
            throw std::runtime_error("get_J: There is no jacobian in the Gauss Seidel method");
        }

        // todo change the name!
        bool compute_pf(
            const EigenRefConstCplxSpMat     & Ybus,
            const Eigen::Ref<const CplxVect> & V,
            const Eigen::Ref<const CplxVect> & Sbus,
            const Eigen::Ref<const IntVect>  & slack_ids,
            const Eigen::Ref<const RealVect> & slack_weights,  // currently unused
            const Eigen::Ref<const IntVect>  & pv,
            const Eigen::Ref<const IntVect>  & pq,
            int                              max_iter,
            real_type                        tol
        ) override;

        void reset() override {
            BaseAlgo::reset();
            // drop the Ybus caches below: they are rebuilt by the next compute_pf
            ybus_rowmajor_ = Eigen::SparseMatrix<cplx_type, Eigen::RowMajor>();
            ybus_diag_ = CplxVect();
        }

    protected:

        virtual
        void one_iter(
            Eigen::Ref<CplxVect>            tmp_Sbus,
            const EigenRefConstCplxSpMat    & Ybus,
            const Eigen::Ref<const IntVect> & pv,
            const Eigen::Ref<const IntVect> & pq
        );

        // Ybus arrives column-major (that is what EigenRefConstCplxSpMat is),
        // but a Gauss-Seidel sweep needs one ROW of it at a time -- and a row of
        // a column-major sparse matrix is not stored anywhere: reading it means
        // scanning every column. `Ybus.row(k) * V_`, called once per pq bus and
        // twice per pv bus, therefore made one sweep cost O(nb_bus * nnz)
        // instead of O(nnz). These two caches, rebuilt once per compute_pf
        // (O(nnz), against the hundreds or thousands of sweeps that follow),
        // give the sweep the row it needs directly:
        //   - ybus_rowmajor_ : the same coefficients, row-major, so row k is a
        //     contiguous InnerIterator walk;
        //   - ybus_diag_     : Ybus(k, k), which Ybus.coeff(k, k) was
        //     re-finding by binary search inside column k every time.
        // The row-major inner iterator visits a row's coefficients by
        // increasing column index, which is the order Eigen's sparse * dense
        // product summed them in, so the sums are unchanged to the last bit.
        void refresh_ybus_cache(const EigenRefConstCplxSpMat & Ybus);

        // row k of Ybus times the current voltages: sum_j Ybus(k, j) * V_(j)
        cplx_type ybus_row_times_V(int k) const {
            cplx_type res(0., 0.);
            for (Eigen::SparseMatrix<cplx_type, Eigen::RowMajor>::InnerIterator it(ybus_rowmajor_, k); it; ++it)
                res += it.value() * V_.coeff(it.col());
            return res;
        }

        Eigen::SparseMatrix<cplx_type, Eigen::RowMajor> ybus_rowmajor_;
        CplxVect                                        ybus_diag_;

    private:
        // no copy allowed
        GaussSeidelAlgo(const GaussSeidelAlgo&) = delete;
        GaussSeidelAlgo(GaussSeidelAlgo&&) = delete;
        GaussSeidelAlgo & operator=(GaussSeidelAlgo&&) = delete;
        GaussSeidelAlgo & operator=(const GaussSeidelAlgo&) = delete;

};


} // namespace ls2g

#endif // GAUSSSEIDEL_ALGO_H