#include <iostream>
#include <vector>
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <complex>      // std::complex, std::conj

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

// import klu package
extern "C" {
    #include "cs.h"
    #include "klu.h"
}

typedef std::complex<double> cdouble;




class CustTimer{
    public:
        CustTimer():start_(std::chrono::system_clock::now()){
            end_ = start_;
        };

        double duration(){
            end_ = std::chrono::system_clock::now();
            std::chrono::duration<double> res = end_ - start_;
            return res.count();
        }
    private:
        std::chrono::time_point<std::chrono::system_clock> start_;
        std::chrono::time_point<std::chrono::system_clock> end_;
//        void start(){
//            start_ = std::chrono::system_clock::now();
//        }
//        void end(){
//            end_ = std::chrono::system_clock::now();
//        }

};

/**
class to handle the solver using newton-raphson method, using KLU algorithm and sparse matrices.

As long as the admittance matrix of the sytem does not change, you can reuse the same solver.
Reusing the same solver is possible, but "reset" method must be called.

Otherwise, unexpected behaviour might follow, including "segfault".

**/
class KLUSolver
{
    public:
        KLUSolver():symbolic_(),numeric_(),common_(),n_(-1),need_factorize_(true),err_(-1),
                    timer_Fx_(0.){
            klu_defaults(&common_);
            timer_Fx_ = 0.;
            timer_solve_ = 0.;
            timer_initialize_ = 0.;
            timer_check_ = 0.;
            timer_dSbus_ = 0.;
            timer_fillJ_ = 0.;
            timer_total_nr_ = 0.;
        }

        ~KLUSolver()
         {
             klu_free_symbolic(&symbolic_, &common_);
             klu_free_numeric(&numeric_, &common_);
         }

        bool do_newton(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       ){
            // TODO check what can be checked: no voltage at 0, Ybus is square, Sbus same size than V and
            // TODO Ybus (nrow or ncol), pv and pq have value that are between 0 and nrow etc.
            if(err_ > 0) return false; // i don't do anything if there were a problem at the initialization
            auto timer = CustTimer();
            // initialize once and for all the "inverse" of these vectors
            int n_pv = pv.size();
            int n_pq = pq.size();
            Eigen::VectorXi pvpq(n_pv + n_pq);
            pvpq << pv, pq;
            int n_pvpq = pvpq.size();
            std::vector<int> pvpq_inv(V.size(), -1);
            for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
            std::vector<int> pq_inv(V.size(), -1);
            for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;

            // first check, if the problem is already solved, i stop there
            Eigen::VectorXd F = _evaluate_Fx(Ybus, V, Sbus, pv, pq);
            bool converged = _check_for_convergence(F, tol);
            nr_iter_ = 0; //current step
            bool res = true;  // have i converged or not
            bool has_just_been_inialized = false;  // to avoid a call to klu_refactor follow a call to klu_factor in the same loop
            while ((!converged) & (nr_iter_ < max_iter)){
                nr_iter_++;
                fill_jacobian_matrix(Ybus, V, pq, pvpq, pq_inv, pvpq_inv);
                if(need_factorize_){
                    initialize();
                    if(err_ != 0){
                        // I got an error during the initialization of the linear system, i need to stop here
                        res = false;
                        break;
                    }
                    has_just_been_inialized = true;
                }
                //TODO refactorize is called uselessly at the first iteration
                solve(F, has_just_been_inialized);
                has_just_been_inialized = false;
                if(err_ != 0){
                    // I got an error during the solving of the linear system, i need to stop here
                    res = false;
                    break;
                }
                auto dx = -1.0*F;

                // update voltage (this should be done consistently with "klu_solver._evaluate_Fx")
                Vm_ = V.array().abs();  // update Vm and Va again in case
                Va_ = V.array().arg();  // we wrapped around with a negative Vm

                if (n_pv > 0) Va_(pv) += dx.segment(0,n_pv);
                if (n_pq > 0){
                    Va_(pq) += dx.segment(n_pv,n_pq);
                    Vm_(pq) += dx.segment(n_pv+n_pq, n_pq);
                }

                // TODO change here for not having to cast all the time ... maybe
                V = Vm_.array() * (Va_.array().cos().cast<cdouble>() + 1.0i * Va_.array().sin().cast<cdouble>() );

                F = _evaluate_Fx(Ybus, V, Sbus, pv, pq);
                converged = _check_for_convergence(F, tol);
            }
            if(!converged){
                err_ = 4;
                res = false;
            }
            timer_total_nr_ += timer.duration();
            return res;
        }

        Eigen::SparseMatrix<double> get_J(){
            return J_;
        }
        Eigen::Ref<Eigen::VectorXd> get_Va(){
            return Va_;
        }
        Eigen::Ref<Eigen::VectorXd> get_Vm(){
            return Vm_;
        }
        int get_error(){
            return err_;
        }
        int get_nb_iter(){
            return nr_iter_;
        }
        void reset(){
            klu_free_symbolic(&symbolic_, &common_);
            klu_free_numeric(&numeric_, &common_);
            n_ = -1;
            common_ = klu_common();

            Vm_ = Eigen::VectorXd();  // voltage magnitude
            Va_= Eigen::VectorXd();  // voltage angle
            J_ = Eigen::SparseMatrix<double>();  // the jacobian matrix
            dS_dVm_ = Eigen::SparseMatrix<cdouble>();
            dS_dVa_ = Eigen::SparseMatrix<cdouble>();
            need_factorize_ = true;
            nr_iter_ = 0;  // number of iteration performs by the Newton Raphson algorithm
            err_ = -1; //error message:

            // reset timers
            timer_Fx_ = 0.;
            timer_solve_ = 0.;
            timer_initialize_ = 0.;
            timer_check_ = 0.;
            timer_dSbus_ = 0.;
            timer_fillJ_ = 0.;
            timer_total_nr_ = 0.;

        }
        bool converged(){
            return err_ == 0;
        }
         std::tuple<double, double, double, double, double, double, double> get_timers()
         {
            auto res = std::tuple<double, double, double, double, double, double, double>(
              timer_Fx_, timer_solve_, timer_initialize_, timer_check_, timer_dSbus_, timer_fillJ_, timer_total_nr_);
            return res;
         }

        std::tuple<Eigen::SparseMatrix<cdouble> , Eigen::SparseMatrix<cdouble> >
                    _get_ds_test(Eigen::SparseMatrix<cdouble> & Ybus,
                                Eigen::VectorXcd & V){
            _dSbus_dV(Ybus, V);
            auto res = std::tuple<Eigen::SparseMatrix<cdouble> , Eigen::SparseMatrix<cdouble> >(dS_dVm_, dS_dVa_);
            return res;
        }

    protected:
        void initialize(){
            // default Eigen representation: column major, which is good for klu !
            // J is const here, even if it's not said in klu_analyze
            auto timer = CustTimer();
            n_ = J_.cols(); // should be equal to J_.nrows()
            err_ = 0; // reset error message
            common_ = klu_common();
            symbolic_ = klu_analyze(n_, J_.outerIndexPtr(), J_.innerIndexPtr(), &common_);
            numeric_ = klu_factor(J_.outerIndexPtr(), J_.innerIndexPtr(), J_.valuePtr(), symbolic_, &common_);
            if (common_.status != KLU_OK) {
                err_ = 1;
            }
            need_factorize_ = false;
            timer_solve_ += timer.duration();
        }

        void solve(Eigen::VectorXd & b, bool has_just_been_inialized){
            // solves (for x) the linear system J.x = b
            // supposes that the solver has been initialized (call klu_solver.analyze() before calling that)
            // J is const even if it does not compile if said const
            auto timer = CustTimer();
            int ok;
            bool stop = false;
            if(!has_just_been_inialized){
                // if the call to "klu_factor" has been made this iteration, there is no need
                // to re factor again the matrix
                // i'm in the case where it has not
                ok = klu_refactor(J_.outerIndexPtr(), J_.innerIndexPtr(), J_.valuePtr(), symbolic_, numeric_, &common_);
                if (ok != 1) {
                    err_ = 2;
                    stop = true;
                }
            }
            if(!stop){
                ok = klu_solve(symbolic_, numeric_, n_, 1, &b(0), &common_);
                if (ok != 1) {
                    err_ = 3;
                }
            }
            timer_solve_ += timer.duration();
        }

        Eigen::SparseMatrix<cdouble>
            _make_diagonal_matrix(const Eigen::Ref<const Eigen::VectorXcd > & diag_val){
            // TODO their might be a more efficient way to do that
            auto n = diag_val.size();
            Eigen::SparseMatrix<cdouble> res(n,n);
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

        void _dSbus_dV(Eigen::SparseMatrix<cdouble> & dS_dVm,
                           Eigen::SparseMatrix<cdouble> & dS_dVa,
                           const Eigen::Ref<const Eigen::SparseMatrix<cdouble> > & Ybus,
                           const Eigen::Ref<const Eigen::VectorXcd > & V)
        {
            // "slow" implementation close to pypower, but with sparse matrix
            // TODO check i cannot optimize that with numba code in pandapower instead
            auto timer = CustTimer();
            Eigen::VectorXcd Ibus = Ybus * V;
            Eigen::SparseMatrix<cdouble> diagV = _make_diagonal_matrix(V);

            Eigen::VectorXcd Ibus_conj = Ibus.conjugate();
            Eigen::SparseMatrix<cdouble> diagIbus_conj = _make_diagonal_matrix(Ibus_conj);

            Eigen::VectorXcd Vnorm = V.array() / V.array().abs();
            Eigen::SparseMatrix<cdouble> diagVnorm = _make_diagonal_matrix(Vnorm);

            Eigen::SparseMatrix<cdouble> tmp = Ybus * diagVnorm;
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
            timer_dSbus_ += timer.duration();
        }

        void _dSbus_dV(const Eigen::Ref<const Eigen::SparseMatrix<cdouble> > & Ybus,
                       const Eigen::Ref<const Eigen::VectorXcd > & V){
            auto timer = CustTimer();
            auto size_dS = V.size();
            Eigen::VectorXcd Vnorm = V.array() / V.array().abs();
            Eigen::VectorXcd Ibus = Ybus * V;

            // TODO see if i can reuse previous values, i am not sure
            dS_dVm_ = 1.0 * Ybus;
            dS_dVa_ = 1.0 * Ybus;

            // i fill the buffer columns per columns
            for (int k=0; k < size_dS; ++k){
                for (Eigen::SparseMatrix<cdouble>::InnerIterator it(dS_dVm_,k); it; ++it)
                {
                    it.valueRef() *= Vnorm(it.col());  // dS_dVm[k] *= Vnorm[Yj[k]]
                    it.valueRef() = std::conj(it.valueRef()) * V(it.row());  // dS_dVm[k] = conj(dS_dVm[k]) * V[r]
                    if(it.col() == it.row()){
                        // diagonal element
                        it.valueRef() += std::conj(Ibus(it.row())) * Vnorm(it.row()); // dS_dVm[k] += buffer[r] # buffer being conj(Ibus * Vnorm)
                    }
                }
            }

            for (int k=0; k < size_dS; ++k){
                for (Eigen::SparseMatrix<cdouble>::InnerIterator it(dS_dVa_,k); it; ++it)
                {
                    it.valueRef() *= V(it.col());  // dS_dVa[k] *= V[Yj[k]]
                    if(it.col() == it.row()){
                        // diagonal element
                        it.valueRef() -= Ibus(it.row());  // dS_dVa[k] = -Ibus[r] + dS_dVa[k]
                    }
                    cdouble tmp = static_cast<cdouble>(1.0i) * V(it.row());
                    it.valueRef() = std::conj(-it.valueRef()) * tmp;  // dS_dVa[k] = conj(-dS_dVa[k]) * (1j * V[r])
                }
            }
            dS_dVa_.makeCompressed();
            dS_dVm_.makeCompressed();
            timer_dSbus_ += timer.duration();
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

        void fill_jacobian_matrix(const Eigen::SparseMatrix<cdouble> & Ybus,
                                  const Eigen::VectorXcd & V,
                                  const Eigen::VectorXi & pq,
                                  const Eigen::VectorXi & pvpq,
                                  const std::vector<int> & pq_inv,
                                  const std::vector<int> & pvpq_inv
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

            auto timer = CustTimer();
            // might also be accelerated in some cases, it's one of the "slow" implementation of pandapower
            // TODO apparently, following matrix have also a fixed shape, so i can be faster !
//            Eigen::SparseMatrix<std::complex< double > > dS_dVm;
//            Eigen::SparseMatrix<std::complex< double > > dS_dVa;
//            _dSbus_dV(dS_dVm_, dS_dVa_, Ybus, V);
            _dSbus_dV(Ybus, V);
            Eigen::SparseMatrix<double> dS_dVa_r = dS_dVa_.real();
            Eigen::SparseMatrix<double> dS_dVa_i = dS_dVa_.imag();
            Eigen::SparseMatrix<double> dS_dVm_r = dS_dVm_.real();
            Eigen::SparseMatrix<double> dS_dVm_i = dS_dVm_.imag();
            // TODO not sure it's mandatory
//            dS_dVa_r.makeCompressed();
//            dS_dVa_i.makeCompressed();
//            dS_dVm_r.makeCompressed();
//            dS_dVm_i.makeCompressed();

            const int n_pvpq = pvpq.size();
            const int n_pq = pq.size();

            const int size_j = n_pvpq + n_pq;

            bool need_insert = false;  // i optimization: i don't need to insert the coefficient in the matrix
            if(J_.cols() != size_j)
            {
                // optim : if the matrix was already computed, i don't initalize it, i instead reuse as much as i can
                // i can do that because the matrix will ALWAYS have the same non zero coefficients.
                // in this if, i allocate it in a "large enough" place to avoid copy when first filling it
                need_insert = true;
                J_ = Eigen::SparseMatrix<double>(size_j,size_j);
                // pre allocate a large enough matrix
                J_.reserve(2*(dS_dVa_.nonZeros()+dS_dVm_.nonZeros()));
                // from an experiment, outerIndexPtr is inialized, with the number of columns
                // innerIndexPtr and valuePtr are not.
            }

            // i fill the buffer columns per columns
            int nb_obj_this_col = 0;
            std::vector<int> inner_index;
            std::vector<double> values;

            // TODO use the loop provided above (in dS) if J is already initialized
            // fill n_pvpq leftmost columns
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
                for(int in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
                    int row_id = inner_index[in_ind];
                    if(need_insert) J_.insert(row_id, col_id) = values[in_ind];
                    else J_.coeffRef(row_id, col_id) = values[in_ind];
                }
            }

            //TODO make same for the second part (have a funciton for previous loop)
            // fill the remaining n_pq columns
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
                for(int in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
                    int row_id = inner_index[in_ind];
                    if(need_insert) J_.insert(row_id, col_id + n_pvpq) = values[in_ind];
                    else J_.coeffRef(row_id, col_id + n_pvpq) = values[in_ind];
                }
            }
            J_.makeCompressed();
            timer_fillJ_ += timer.duration();
        }

        Eigen::VectorXd _evaluate_Fx(const Eigen::SparseMatrix<cdouble> &  Ybus,
                                     const Eigen::VectorXcd & V,
                                     const Eigen::VectorXcd & Sbus,
                                     const Eigen::VectorXi & pv,
                                     const Eigen::VectorXi & pq){
            auto timer = CustTimer();
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
            timer_Fx_ += timer.duration();
            return res;
        }

         bool _check_for_convergence(const Eigen::VectorXd & F,
                                     double tol){
            auto timer = CustTimer();
            return F.lpNorm<Eigen::Infinity>()  < tol;
            timer_check_ += timer.duration();
         }

    private:
        // solver initialization
        klu_symbolic* symbolic_;
        klu_numeric* numeric_;
        klu_common common_;
        int n_;

        // solution of the problem
        Eigen::VectorXd Vm_;  // voltage magnitude
        Eigen::VectorXd Va_;  // voltage angle
        Eigen::SparseMatrix<double> J_;  // the jacobian matrix
        Eigen::SparseMatrix<cdouble> dS_dVm_;
        Eigen::SparseMatrix<cdouble> dS_dVa_;
        bool need_factorize_;
        int nr_iter_;  // number of iteration performs by the Newton Raphson algorithm
        int err_; //error message:
        // -1 : the solver has not been initialized (call initialize in this case)
        // 0 everything ok
        // 1: i can't factorize the matrix (klu_factor)
        // 2: i can't refactorize the matrix (klu_refactor)
        // 3: i can't solve the system (klu_solve)
        // 4: end of possible iterations (divergence because nr_iter_ >= max_iter

        // timers
         double timer_Fx_;
         double timer_solve_;
         double timer_initialize_;
         double timer_check_;
         double timer_dSbus_;
         double timer_fillJ_;
         double timer_total_nr_;

        // no copy allowed
        KLUSolver( const KLUSolver & ) ;
        KLUSolver & operator=( const KLUSolver & ) ;

        // debug func i don't want to remove yet
        void analyze_old(int n,
                      Eigen::Ref<Eigen::VectorXi> Ap,
                      Eigen::Ref<Eigen::VectorXi> Ai){
            n_ = n;
            symbolic_ = klu_analyze(n, &Ap(0), &Ai(0), &common_);
        }
        void solve_old(Eigen::Ref<Eigen::VectorXi> Ap,
                    Eigen::Ref<Eigen::VectorXi> Ai,
                    Eigen::Ref<Eigen::VectorXd> Ax,
                    Eigen::Ref<Eigen::VectorXd> b){
            numeric_ = klu_factor(&Ap(0), &Ai(0), &Ax(0), symbolic_, &common_);
            klu_solve(symbolic_, numeric_, n_, 1, &b(0), &common_);
        }

                 // TODO re add the references here for the last stuff
        std::tuple<Eigen::VectorXd, Eigen::VectorXcd> one_iter_test(Eigen::SparseMatrix<double> J,
                              Eigen::Ref<Eigen::VectorXd> F,
                              Eigen::VectorXi pv,
                              Eigen::VectorXi pq,
                              Eigen::VectorXcd V,
                              Eigen::SparseMatrix<cdouble>  Ybus,
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
            V = Vm_.array() * (Va_.array().cos().cast<cdouble>() + 1.0i * Va_.array().sin().cast<cdouble>() );

            F = _evaluate_Fx(Ybus, V, Sbus, pv, pq);
            return std::tuple<Eigen::VectorXd, Eigen::VectorXcd>(F, V);
        }

        Eigen::SparseMatrix<double>
             create_jacobian_matrix_test(const Eigen::SparseMatrix<cdouble> & Ybus,
                                         const Eigen::VectorXcd & V,
                                         const Eigen::VectorXi & pq,
                                         const Eigen::VectorXi & pvpq
                                         ){

            // DO NOT USE, FOR DEBUG ONLY!
            int n_pvpq = pvpq.size();
            int n_pq = pvpq.size();
            std::vector<int> pvpq_inv(V.size(), -1);
            for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
            std::vector<int> pq_inv(V.size(), -1);
            for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;
            fill_jacobian_matrix(Ybus, V, pq, pvpq, pq_inv, pvpq_inv);
            return J_;
        }

        bool initialize_test(Eigen::SparseMatrix<double > & J){
            // default Eigen representation: column major, which is good for klu !
            // J is const here, even if it's not said in klu_analyze
            int n = J.cols(); // should be equal to J_.nrows()
            err_ = 0; // reset error message
            klu_common common = klu_common();
            bool res = true;
            klu_symbolic* symbolic = klu_analyze(n, J.outerIndexPtr(), J.innerIndexPtr(), &common);
//            auto numeric = klu_factor(J.outerIndexPtr(), J.innerIndexPtr(), J.valuePtr(), symbolic, &common);
            klu_factor(J.outerIndexPtr(), J.innerIndexPtr(), J.valuePtr(), symbolic, &common);
            if (common_.status != KLU_OK) {
                err_ = 1;
                res = false;
            }
            return res;
        }

};

