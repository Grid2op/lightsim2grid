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
            res.setIdentity();  // segfault if attempt to use it
            res.diagonal() = diag_val;

            // second method, with an "optimized" loop
            // res.reserve(Eigen::VectorXi::Constant(n,1)); // i reserve one number per columns (speed optim)
            // for(unsigned int i = 0; i < n; ++i){
            //    res.insert(i,i) = diag_val(i);
            //}
            return res;
        }

        Eigen::SparseMatrix<double>
             create_jacobian_matrix(const Eigen::SparseMatrix<std::complex< double > > Ybus,
                                    const Eigen::VectorXcd V,
                                    const Eigen::VectorXi pq,
                                    const Eigen::VectorXi pvpq
//                                    const Eigen::Ref<const Eigen::SparseMatrix<std::complex< double > > > & Ybus,
//                                    const Eigen::Ref<const Eigen::VectorXcd > & V,
//                                    const Eigen::Ref<const Eigen::VectorXi > & pq,
//                                    const Eigen::Ref<const Eigen::VectorXi > & pvpq
                                    ){
            // might also be accelerated in some cases, it's one of the "slow" implementation of pandapower
            // TODO apparently, following matrix have also a fixed shape, so i can be faster !
            Eigen::SparseMatrix<std::complex< double > > dS_dVm;
            Eigen::SparseMatrix<std::complex< double > > dS_dVa;
            _dSbus_dV(dS_dVm, dS_dVa, Ybus, V);
            auto n_pvpq = pvpq.size();
            auto n_pq = pq.size();
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

            int size_j = n_pvpq + n_pq;
            // TODO here: don't start from scratch each time, reuse the former J to replace it's coefficient
            // TODO this is rather slow, see if i cannot be smarter than this! (look at pandapower code)
            Eigen::SparseMatrix<double> J(size_j,size_j);
            // from an experiment, outerIndexPtr is inialized, with the number of columns
            // innerIndexPtr and valuePtr are not.

            // i fill the buffer columns per columns
            int nb_obj_this_col = 0;
            int start_id, end_id;
            std::vector<int> inner_index;
            std::vector<double> values;

            Eigen::SparseMatrix<double> dS_dVa_r = dS_dVa.real();
            Eigen::SparseMatrix<double> dS_dVa_i = dS_dVa.imag();
            Eigen::SparseMatrix<double> dS_dVm_r = dS_dVm.real();
            Eigen::SparseMatrix<double> dS_dVm_i = dS_dVm.imag();
            dS_dVa_r.makeCompressed();
            dS_dVa_i.makeCompressed();
            dS_dVa.makeCompressed();
            dS_dVm_r.makeCompressed();
            dS_dVm_i.makeCompressed();
            //TODO Check the stuff bellow, and execute it once and for all only for all iterations!
            std::vector<int> pvpq_inv(n_pvpq, -1);
            for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
            std::vector<int> pq_inv(n_pq, -1);
            for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;

//            double test = J.coeffRef(size_j -1 , size_j -1);
//            test = test + 1;
//            int toto = J.outerIndexPtr()[size_j-1];
//            std::cout << "test " << toto << std::endl;
            for(int col_id=0; col_id < n_pvpq; ++col_id){
                nb_obj_this_col = 0;

                // fill with the first column with the column of dS_dVa[:,pvpq[col_id]]
                // and check the row order !
                int col_id_dS_dVa = pvpq[col_id];
                start_id = dS_dVa_r.outerIndexPtr()[col_id_dS_dVa];
                end_id = dS_dVa_r.outerIndexPtr()[col_id_dS_dVa+1];
//                std::cout << "\tstart_id " << start_id << " end_id " << end_id <<std::endl;
                // TODO try to optimize that here push back is meeeeh
                for(int obj_id = start_id; obj_id < end_id; ++obj_id)
                {
//                    std::cout << "check " << obj_id << std::endl;
                    int row_id_dS_dVa = dS_dVa_r.valuePtr()[obj_id];
                    // I add the value only if the rows was selected in the indexes
                    int row_id = pvpq_inv[row_id_dS_dVa];
                    if(row_id > 0)
                    {
                        inner_index.push_back(pvpq_inv[row_id_dS_dVa]);
                        values.push_back(dS_dVa_r.valuePtr()[obj_id]);
                        nb_obj_this_col++;
                    }
                }

                // fill the rest of the rows with the first column of dS_dVa_imag[:,pq[col_id]]
                // TODO refactor with previous loop

                // indicate number of elements put in the matrix
//                std::cout << "col_id " << col_id << " size_j " << size_j <<std::endl;
                J.outerIndexPtr()[col_id] = nb_obj_this_col;
            }

            //TODO make same for the second part (have a funciton for previous loop)

            // and now set the value computed in J
            std::cout << "before "<<std::endl;
            double * valuePtr = new double[values.size()];  // deletion is ensured by the J sparse matrix
            memcpy(valuePtr, &values[0], values.size());
            std::cout << "valuePtr " << valuePtr <<std::endl;
            J.valuePtr() = valuePtr;

//            int * innerIndexPtr = new int[values.size()]; // deletion is ensured by the J sparse matrix
//            memcpy(innerIndexPtr, &inner_index[0], inner_index.size());
//            J.innerIndexPtr() = innerIndexPtr;
            std::cout << "after "<<std::endl;


//            J.setIdentity();
//            std::cout << "J.valuePtr() " << J.valuePtr() << std::endl;
//            std::cout << "J.innerIndexPtr() " << J.innerIndexPtr() << std::endl;
//            std::cout << "J.outerIndexPtr() " << J.outerIndexPtr() << std::endl;
            // this is returned to check that i got the same results as pandapower
            return J;
            // Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permforJ11(pvpq);
            // Eigen::SparseMatrix<double> J11 = dS_dVm.real() * permforJ11;
            //reorder the row / coumn of the matrix

//            Eigen::Block<Eigen::SparseMatrix<double>, Eigen::Dynamic, Eigen::Dynamic > J11 = J.block(0,0, n_pvpq, n_pvpq);

//            auto & J12 = J.block(0,n_pvpq, n_pvpq, n_pq);
//            auto & J21 = J.block(n_pvpq, 0, n_pq, n_pvpq);
//            auto & J22 = J.block(n_pvpq, n_pvpq, n_pq, n_pq);

//            for(unsigned int i_ = 0; i_ < n_pvpq; ++i_){
//                int i = pvpq(i_);
//                J.block(0,0, n_pvpq, n_pvpq).col(i) = dS_dVa.col(i).real();
//                J.coeffRef(i,i) = dS_dVa(i,i).real();
//                J.coeffRef(i,i + n_pvpq) = dS_dVm(i,i).imag(i);
//            }
        }
        void _dSbus_dV(Eigen::SparseMatrix<std::complex< double > > & dS_dVm,
                       Eigen::SparseMatrix<std::complex< double > > & dS_dVa,
                       const Eigen::Ref<const Eigen::SparseMatrix<std::complex< double > > > & Ybus,
                       const Eigen::Ref<const Eigen::VectorXcd > & V)
        {
            // "slow" implementation close to pypower, but with sparse matrix
            // TODO check i cannot optimize that with numba code in pandapower
            Eigen::VectorXcd Ibus = Ybus * V;  // vector
            Eigen::VectorXcd _1_vabs = 1. / V.array().abs();

            Eigen::SparseMatrix<std::complex< double > > diagV = _make_diagonal_matrix(V);
            auto Ibus_conj = Ibus.conjugate();
            Eigen::SparseMatrix<std::complex< double > > diagIbus_conj = _make_diagonal_matrix(Ibus_conj);
            auto Vnorm = V * _1_vabs;
            Eigen::SparseMatrix<std::complex< double > > diagVnorm = _make_diagonal_matrix(Vnorm);

            Eigen::SparseMatrix<std::complex< double > > tmp = Ybus * diagVnorm;
            tmp = tmp.conjugate();
            dS_dVm = diagV * tmp + Ibus_conj * diagVnorm;  //diagIbus_conj * diagVnorm;

            Eigen::SparseMatrix<std::complex< double > > tmp2 = diagIbus_conj - tmp;  // conj(diagIbus - Ybus * diagV)
            dS_dVa = 1.0i * diagV * tmp2;

            // python implementation
            // dS_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm
            // dS_dVa = 1j * diagV * conj(diagIbus - Ybus * diagV)
//            return std::tuple<Eigen::SparseMatrix<std::complex< double > >,
//                              Eigen::SparseMatrix<std::complex< double > > >(dS_dVm, dS_dVa);
        }


         // TODO re add the references here for the last stuff
         std::tuple<Eigen::VectorXd, Eigen::VectorXcd> one_iter(Eigen::SparseMatrix<double> J,
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
            // TODO init Va and Vm only once at the end
            Eigen::VectorXd Vm = V.array().abs();  // update Vm and Va again in case
            Eigen::VectorXd Va = V.array().arg();  // we wrapped around with a negative Vm
            if (npv > 0) Va(pv) += dx.segment(0,npv);
            if (npq > 0){
                Va(pq) += dx.segment(npv,npq);
                Vm(pq) += dx.segment(npv+npq, npq);
            }

            // TODO change here for not having to cast all the time ...
            V = Vm.array() * (Va.array().cos().cast<std::complex< double > >() + 1.0i * Va.array().sin().cast<std::complex< double > >() );

            F = _evaluate_Fx(Ybus, V, Sbus, pv, pq);
            return std::tuple<Eigen::VectorXd, Eigen::VectorXcd>(F, V);
          }

        // TODO re add the references here for the last stuff
        Eigen::VectorXd _evaluate_Fx(Eigen::SparseMatrix<std::complex< double > > &  Ybus,
                                        Eigen::VectorXcd V,
                                        Eigen::VectorXcd Sbus,
                                        Eigen::VectorXi pv,
                                        Eigen::VectorXi pq){

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
        klu_symbolic* symbolic_;
        klu_numeric* numeric_;
        klu_common common_;
        int n_;

        // no copy allowed
        KLUSolver( const KLUSolver & ) ;
        KLUSolver & operator=( const KLUSolver & ) ;
};

class Dog{
    public:
        Dog(){};
        void bark(){};
};

