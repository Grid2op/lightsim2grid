#include <vector>
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

         void analyze(Eigen::SparseMatrix<double> J){
         // default Eigen representation: column major, which is good for klu !
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
         void solve(Eigen::SparseMatrix<double> J,
                    Eigen::Ref<Eigen::VectorXd> b
                    ){
            // solves the linear system M.x = b
            // supposes that the solver has been initialized (call klu_solver.analyze(M) before calling that)
            numeric_ = klu_factor(J.outerIndexPtr(), J.innerIndexPtr(), J.valuePtr(), symbolic_, &common_);
            klu_solve(symbolic_, numeric_, n_, 1, &b(0), &common_);
         }

        std::tuple<Eigen::SparseMatrix<std::complex< double > >,
                   Eigen::SparseMatrix<std::complex< double > > >
                   _dSbus_dV(Eigen::SparseMatrix<std::complex< double > > Ybus,
                            Eigen::VectorXcd V)
        {
            // "slow" implementation close to pypower, but with sparse matrix
            Eigen::VectorXcd Ibus = Ybus * V;  // vector
            Eigen::VectorXcd _1_vabs = 1. / V.array().abs();

            Eigen::DiagonalMatrix<std::complex< double >, Eigen::Dynamic > diagV(V);
            auto Ibus_conj = Ibus.conjugate();
            Eigen::DiagonalMatrix<std::complex< double >, Eigen::Dynamic > diagIbus_conj(Ibus_conj);
            Eigen::DiagonalMatrix<std::complex< double >, Eigen::Dynamic > diagVnorm(V * _1_vabs);

            Eigen::SparseMatrix<std::complex< double > > tmp = Ybus * diagV;
            tmp = tmp.conjugate();
            Eigen::SparseMatrix<std::complex< double > > dS_dVm = diagV * tmp *_1_vabs + Ibus_conj * V * _1_vabs;  //diagIbus_conj * diagVnorm;

            Eigen::SparseMatrix<std::complex< double > > tmp2 = diagIbus_conj - tmp;  // conj(diagIbus - Ybus * diagV)
            Eigen::SparseMatrix<std::complex< double > > dS_dVa = 1.0i * diagV * tmp2;

            // python implementation
            // dS_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm
            // dS_dVa = 1j * diagV * conj(diagIbus - Ybus * diagV)
            return std::tuple<Eigen::SparseMatrix<std::complex< double > >,
                              Eigen::SparseMatrix<std::complex< double > > >(dS_dVm, dS_dVa);
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
            res.segment(npv,npq)  = real_(pq);
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

