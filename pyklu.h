#include <iostream>
#include <vector>
#include <stdio.h>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "cs.h"
#include "klu.h"

// eigen is necessary to easily pass data from numpy to c++ without any copy.

// class to handle the memory of the system
class KLUSolver
{
    public:
        KLUSolver():symbolic_(),numeric_(),common_(),n_(-1){
            klu_defaults(&common_);
        }

        ~KLUSolver()
         {
             klu_free_symbolic(&symbolic_, &common_);
             klu_free_numeric(&numeric_, &common_);
         }

         void analyze_old(int n,
                      Eigen::Ref<Eigen::VectorXi> Ap,
                      Eigen::Ref<Eigen::VectorXi> Ai){
            n_ = n;
            symbolic_ = klu_analyze(n, &Ap(0), &Ai(0), &common_);
         }

         void analyze(Eigen::SparseMatrix<double> & J){
            // default Eigen representation: column major, which is good for klu !
            // J is const here, even if it's not said in klu_analyze
            n_ = J.cols(); // should be equal to
            symbolic_ = klu_analyze(n_, J.outerIndexPtr(), J.innerIndexPtr(), &common_);
         }
         void solve_old(Eigen::Ref<Eigen::VectorXi> Ap,
                    Eigen::Ref<Eigen::VectorXi> Ai,
                    Eigen::Ref<Eigen::VectorXd> Ax,
                    Eigen::Ref<Eigen::VectorXd> b){
            numeric_ = klu_factor(&Ap(0), &Ai(0), &Ax(0), symbolic_, &common_);
            klu_solve(symbolic_, numeric_, n_, 1, &b(0), &common_);
          }

         void solve(Eigen::SparseMatrix<double> & J,
                    Eigen::Ref<Eigen::VectorXd> b
                    ){
            // solves (for x) the linear system J.x = b
            // supposes that the solver has been initialized (call klu_solver.analyze(M) before calling that)
            // J is const even if it does not compile if said const
            numeric_ = klu_factor(J.outerIndexPtr(), J.innerIndexPtr(), J.valuePtr(), symbolic_, &common_);
            klu_solve(symbolic_, numeric_, n_, 1, &b(0), &common_);
         }

        Eigen::SparseMatrix<std::complex< double > > _make_diagonal_matrix(Eigen::Ref<const Eigen::VectorXcd > diag_val){
            // TODO their might be a more efficient way to do that
            auto n = diag_val.size();
            Eigen::SparseMatrix<std::complex< double > > res(n,n);
            // first method, without a loop of mine
            res.setIdentity();  // segfault if attempt to use this function without this
            res.diagonal() = diag_val;

            // second method, with an "optimized" loop
            // res.reserve(Eigen::VectorXi::Constant(n,1)); // i reserve one number per columns (speed optim)
            // for(unsigned int i = 0; i < n; ++i){
            //    res.insert(i,i) = diag_val(i);
            //}
            return res;
        }

        void _dSbus_dV(Eigen::SparseMatrix<std::complex< double > > & dS_dVm,
                       Eigen::SparseMatrix<std::complex< double > > & dS_dVa,
                       const Eigen::Ref<const Eigen::SparseMatrix<std::complex< double > > > & Ybus,
                       const Eigen::Ref<const Eigen::VectorXcd > & V)
        {
            // "slow" implementation close to pypower, but with sparse matrix
            // TODO check i cannot optimize that with numba code in pandapower instead
            Eigen::VectorXcd Ibus = Ybus * V;
            Eigen::SparseMatrix<std::complex< double > > diagV = _make_diagonal_matrix(V);

            Eigen::VectorXcd Ibus_conj = Ibus.conjugate();
            Eigen::SparseMatrix<std::complex< double > > diagIbus_conj = _make_diagonal_matrix(Ibus_conj);

            Eigen::VectorXcd Vnorm = V.array() / V.array().abs();
            Eigen::SparseMatrix<std::complex< double > > diagVnorm = _make_diagonal_matrix(Vnorm);

            Eigen::SparseMatrix<std::complex< double > > tmp = Ybus * diagVnorm;
            tmp = tmp.conjugate();
            dS_dVm = diagV * tmp + diagIbus_conj * diagVnorm;

            // TODO this is the same code for this "tmp" and the previous one, except for "diagV" and "diagVnorm"
            tmp = Ybus * diagV;
            tmp = tmp.conjugate();
            auto tmp2 = diagIbus_conj - tmp;  // conj(diagIbus - Ybus * diagV)
            dS_dVa = 1.0i * diagV * tmp2;

            // python implementation
            // dS_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm
            // dS_dVa = 1j * diagV * conj(diagIbus - Ybus * diagV)
        }

        void _get_values_J(int & nb_obj_this_col,
                           std::vector<int> & inner_index,
                           std::vector<double> & values,
                           const Eigen::SparseMatrix<double> & mat,  // ex. dS_dVa_r
                           const std::vector<int> & index_row_inv, // ex. pvpq_inv
                           const Eigen::VectorXi & index_col, // ex. pvpq
                           int col_id,
                           int row_lag  // 0 for J11 for example, n_pvpq for J12
                           )
        {
            /**
            This function will fill the "inner_index" and "values" with the non zero values
            present in the matrix "mat" for the column of the J matrix with id "col_id"
            which corresponds to the column "index_col(col_id)" of the matrix mat.

            The rows need to be converted using another vector too. For example, row "j" of J
            need to be filled with element k of matrix "mat" with k such that "index_row[j] = k"
            Hence, we pass as the argument of this function the "inverse" of index_row, which is such
            that : "j = index_row_inv[k]" is easily computable given k.
            **/
            int col_id_mat = index_col(col_id);

            int start_id = mat.outerIndexPtr()[col_id_mat];
            int end_id = mat.outerIndexPtr()[col_id_mat+1];
            for(int obj_id = start_id; obj_id < end_id; ++obj_id)
            {
                int row_id_dS_dVa = mat.innerIndexPtr()[obj_id];
                // I add the value only if the rows was selected in the indexes
                int row_id = index_row_inv[row_id_dS_dVa];
                if(row_id >= 0)
                {
                    inner_index.push_back(row_id+row_lag);
                    values.push_back(mat.valuePtr()[obj_id]);
                    nb_obj_this_col++;
                }
            }
        }

        void do_newton(const Eigen::SparseMatrix<std::complex< double > > & Ybus,
                       const Eigen::VectorXcd & V,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       const Eigen::VectorXi & pvpq,
                       int max_iter){

            // initialize once and for all the "inverse" of these vector
            std::vector<int> pvpq_inv(dS_dVm.size(), -1);
            for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
            std::vector<int> pq_inv(dS_dVm.size(), -1);
            for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;

            nr_iter = 0; //current step
            

        }
        Eigen::SparseMatrix<double>
             create_jacobian_matrix_test(const Eigen::SparseMatrix<std::complex< double > > & Ybus,
                                         const Eigen::VectorXcd & V,
                                         const Eigen::VectorXi & pq,
                                         const Eigen::VectorXi & pvpq
                                         ){
            /**
            J has the shape
            | J11 | J12 |               | (pvpq, pvpq) | (pvpq, pq) |
            | --------- | = dimensions: | ------------------------- |
            | J21 | J22 |               |  (pq, pvpq)  | (pvpq, pq) |
            python implementation:
            J11 = dS_dVa[array([pvpq]).T, pvpq].real
            J12 = dS_dVm[array([pvpq]).T, pq].real
            J21 = dS_dVa[array([pq]).T, pvpq].imag
            J22 = dS_dVm[array([pq]).T, pq].imag
            **/

            // DO NOT USE, FOR DEBUG ONLY!

            // might also be accelerated in some cases, it's one of the "slow" implementation of pandapower
            // TODO apparently, following matrix have also a fixed shape, so i can be faster !
            Eigen::SparseMatrix<std::complex< double > > dS_dVm;
            Eigen::SparseMatrix<std::complex< double > > dS_dVa;
            _dSbus_dV(dS_dVm, dS_dVa, Ybus, V);
            const int n_pvpq = pvpq.size();
            const int n_pq = pq.size();

            const int size_j = n_pvpq + n_pq;
            // TODO here: don't start from scratch each time, reuse the former J to replace it's coefficient
            // TODO this is rather slow, see if i cannot be smarter than this! (look at pandapower code)
            Eigen::SparseMatrix<double> J_(size_j,size_j);
            // pre allocate a large enough matrix
            J_.reserve(2*(dS_dVa.nonZeros()+dS_dVm.nonZeros()));
            // from an experiment, outerIndexPtr is inialized, with the number of columns
            // innerIndexPtr and valuePtr are not.

            dS_dVa.makeCompressed();
            Eigen::SparseMatrix<double> dS_dVa_r = dS_dVa.real();
            Eigen::SparseMatrix<double> dS_dVa_i = dS_dVa.imag();
            Eigen::SparseMatrix<double> dS_dVm_r = dS_dVm.real();
            Eigen::SparseMatrix<double> dS_dVm_i = dS_dVm.imag();
            // TODO not sure it's mandatory
            dS_dVa_r.makeCompressed();
            dS_dVa_i.makeCompressed();
            dS_dVm_r.makeCompressed();
            dS_dVm_i.makeCompressed();

            //TODO Check the stuff bellow, and execute it once and for all only for all iterations!
            std::vector<int> pvpq_inv(dS_dVm.size(), -1);
            for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;

            std::vector<int> pq_inv(dS_dVm.size(), -1);
            for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;

            // i fill the buffer columns per columns
            int nb_obj_this_col = 0;
            std::vector<int> inner_index;
            std::vector<double> values;

            for(int col_id=0; col_id < n_pvpq; ++col_id){
                // reset from the previous column
                nb_obj_this_col = 0;
                inner_index.clear();
                values.clear();

                // fill with the first column with the column of dS_dVa[:,pvpq[col_id]]
                // and check the row order !
                _get_values_J(nb_obj_this_col, inner_index, values,
                              dS_dVa_r,
                              pvpq_inv, pvpq,
                              col_id, 0);
                // fill the rest of the rows with the first column of dS_dVa_imag[:,pq[col_id]]
                _get_values_J(nb_obj_this_col, inner_index, values,
                              dS_dVa_i,
                              pq_inv, pvpq,
                              col_id, n_pvpq
                              );

                // "efficient" insert of the element in the matrix
                J_.reserve(Eigen::VectorXi::Constant(size_j,nb_obj_this_col));
                for(int in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
                    int row_id = inner_index[in_ind];
                    J_.insert(row_id, col_id) = values[in_ind];
                }
            }

            //TODO make same for the second part (have a funciton for previous loop)
            for(int col_id=0; col_id < n_pq; ++col_id){
                // reset from the previous column
                nb_obj_this_col = 0;
                inner_index.clear();
                values.clear();

                // fill with the first column with the column of dS_dVa[:,pvpq[col_id]]
                // and check the row order !
                _get_values_J(nb_obj_this_col, inner_index, values,
                              dS_dVm_r,
                              pvpq_inv, pq,
                              col_id, 0);

                // fill the rest of the rows with the first column of dS_dVa_imag[:,pq[col_id]]
                _get_values_J(nb_obj_this_col, inner_index, values,
                              dS_dVm_i,
                              pq_inv, pq,
                              col_id, n_pvpq
                              );

                // "efficient" insert of the element in the matrix
                J_.reserve(Eigen::VectorXi::Constant(size_j,nb_obj_this_col));
                for(int in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
                    int row_id = inner_index[in_ind];
                    J_.insert(row_id, col_id + n_pvpq) = values[in_ind];
                }
            }
            return J_;
        }

         // TODO re add the references here for the last stuff
         std::tuple<Eigen::VectorXd, Eigen::VectorXcd> one_iter_test(Eigen::SparseMatrix<double> J,
                                  Eigen::Ref<Eigen::VectorXd> F,
                                  Eigen::VectorXi pv,
                                  Eigen::VectorXi pq,
                                  Eigen::VectorXcd V,
                                  Eigen::SparseMatrix<std::complex< double > >  Ybus,
                                  Eigen::VectorXcd Sbus
                                  ){
            //TODO do not use, for DEBUG only!!!
            // get the sizes for convenience
            auto npv = pv.size();
            auto npq = pq.size();

            // supposes that "klu_analyze" has been already called ! klu_solver.analyze(J) should have been called
            numeric_ = klu_factor(J.outerIndexPtr(), J.innerIndexPtr(), J.valuePtr(), symbolic_, &common_);
            klu_solve(symbolic_, numeric_, n_, 1, &F(0), &common_);
            auto dx = -1.0*F;

            // update voltage (this should be done consistently with "klu_solver._evaluate_Fx")
            Vm_ = V.array().abs();  // update Vm and Va again in case
            Va_ = V.array().arg();  // we wrapped around with a negative Vm

            if (npv > 0) Va_(pv) += dx.segment(0,npv);
            if (npq > 0){
                Va_(pq) += dx.segment(npv,npq);
                Vm_(pq) += dx.segment(npv+npq, npq);
            }

            // TODO change here for not having to cast all the time ...
            V = Vm_.array() * (Va_.array().cos().cast<std::complex< double > >() + 1.0i * Va_.array().sin().cast<std::complex< double > >() );

            F = _evaluate_Fx(Ybus, V, Sbus, pv, pq);
            return std::tuple<Eigen::VectorXd, Eigen::VectorXcd>(F, V);
          }

        Eigen::VectorXd _evaluate_Fx(Eigen::SparseMatrix<std::complex< double > > &  Ybus,
                                     Eigen::VectorXcd & V,
                                     Eigen::VectorXcd & Sbus,
                                     Eigen::VectorXi & pv,
                                     Eigen::VectorXi & pq){

            auto npv = pv.size();
            auto npq = pq.size();

            // compute the mismatch
            Eigen::VectorXcd tmp = Ybus * V;  // this is a vector
            tmp = tmp.array().conjugate();  // i take the conjugate
            auto mis = V.array() * tmp.array() - Sbus.array();
            auto real_ = mis.real();
            auto imag_ = mis.imag();

            // build and fill the result
            Eigen::VectorXd res(npv + 2*npq);
            res.segment(0,npv) = real_(pv);
            res.segment(npv,npq) = real_(pq);
            res.segment(npv+npq, npq) = imag_(pq);
            return res;
          }

         bool _check_for_convergence( Eigen::Ref<Eigen::VectorXd> F, double tol){
            return F.lpNorm<Eigen::Infinity>()  < tol;
        }

    private:
        // solver initialization
        klu_symbolic* symbolic_;
        klu_numeric* numeric_;
        klu_common common_;
        int n_;

        // solution of the problem
        Eigen::VectorXd Vm_;
        Eigen::VectorXd Va_;
        Eigen::SparseMatrix<double> J_;
        int nr_iter;

        // no copy allowed
        KLUSolver( const KLUSolver & ) ;
        KLUSolver & operator=( const KLUSolver & ) ;
};

class Dog{
    public:
        Dog(){};
        void bark(){};
};

