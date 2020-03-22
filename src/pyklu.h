#include <iostream>
#include <vector>
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <complex>      // std::complex, std::conj
#include <cmath>  // for PI

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
            /**
            This method uses the newton raphson algorithm to compute voltage angles and magnitudes at each bus
            of the system.
            If the Ybus matrix changed, please uses the appropriate method to recomptue it!
            **/
            // TODO check what can be checked: no voltage at 0, Ybus is square, Sbus same size than V and
            // TODO Ybus (nrow or ncol), pv and pq have value that are between 0 and nrow etc.
            reset_timer();
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
            reset_timer();
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
        void reset_timer(){
            timer_Fx_ = 0.;
            timer_solve_ = 0.;
            timer_initialize_ = 0.;
            timer_check_ = 0.;
            timer_dSbus_ = 0.;
            timer_fillJ_ = 0.;
            timer_total_nr_ = 0.;
        }
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
            dS_dVm_ = Ybus;
            dS_dVa_ = Ybus;

            // i fill the buffer columns per columns
            for (int k=0; k < size_dS; ++k){
                for (Eigen::SparseMatrix<cdouble>::InnerIterator it(dS_dVm_,k); it; ++it)
                {
                    it.valueRef() *= Vnorm(it.col());  // dS_dVm[k] *= Vnorm[Yj[k]]
                    it.valueRef() = std::conj(it.valueRef()) * V(it.row());  // dS_dVm[k] = conj(dS_dVm[k]) * V[r]
                    if(it.col() == it.row()){
                        // diagonal element
                        it.valueRef() += std::conj(Ibus(it.row())) * Vnorm(it.row()); // dS_dVm[k] += buffer[r] # buffer being conj(Ibus) * Vnorm
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
//                if(need_insert){
//                    for(int in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
//                        int row_id = inner_index[in_ind];
//                        J_.insert(row_id, col_id) = values[in_ind];
//                        // else J_.coeffRef(row_id, col_id) = values[in_ind];
//                    }
//                }else{
//                    int in_ind=0;
//                    for (Eigen::SparseMatrix<double>::InnerIterator it(J_,col_id); it; ++it, ++in_ind)
//                    {
////                        int row_id = inner_index[it.row()];
//                        it.valueRef() = values[it.row()];
//                    }
//                }
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
//                if(need_insert){
//                    for(int in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
//                        int row_id = inner_index[in_ind];
//                        J_.insert(row_id, col_id + n_pvpq) = values[in_ind];
//                        // else J_.coeffRef(row_id, col_id) = values[in_ind];
//                    }
//                }else{
//                    int in_ind=0;
//                    for (Eigen::SparseMatrix<double>::InnerIterator it(J_, col_id + n_pvpq); it; ++it, ++in_ind)
//                    {
//                        int row_id = inner_index[it.row()];
//                        it.valueRef() = values[row_id];
//                    }
//                }
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

class DataModel{
    public:
        DataModel():sn_mva_(0.),f_hz_(0.){};
        void set_f_hz(double f_hz) { f_hz_ = f_hz;}
        void set_sn_mva(double sn_mva) { sn_mva_ = sn_mva;}

        void init_bus(Eigen::VectorXd bus_vn_kv, int nb_line, int nb_trafo){
            /**
            initialize the bus_vn_kv_ member
            and
            initialize the Ybus_ matrix at the proper shape
            **/
            int nb_bus = bus_vn_kv.size();
            bus_vn_kv_ = bus_vn_kv;
            Ybus_ = Eigen::SparseMatrix<cdouble>(nb_bus, nb_bus);
            Ybus_.reserve(nb_bus + 2*nb_line + 2*nb_trafo);

            // per unit conversion
            bus_pu_ = Eigen::VectorXd::Constant(nb_bus, 1.0 / sn_mva_);
            bus_pu_.array() *= bus_vn_kv_.array() * bus_vn_kv_.array(); // np.square(base_kv) / net.sn_mva
        }

        Eigen::SparseMatrix<cdouble> get_Ybus(){
            return Ybus_;
        }

        void init_powerlines(const Eigen::VectorXd & branch_r,
                             const Eigen::VectorXd & branch_x,
                             const Eigen::VectorXd & branch_c,
                             const Eigen::VectorXd & branch_g,
                             const Eigen::VectorXi & branch_from_id,
                             const Eigen::VectorXi & branch_to_id
                             )
        {
            /**
            This method initialize the Ybus matrix from the branch matrix.
            It has to be called once when the solver is initialized. Afterwards, a call to
            "updateYbus" should be made for performance optimiaztion instead. //TODO
            **/

            // TODO check what can be checked: branch_* have same size, no id in branch_to_id that are
            // TODO not in [0, .., buv_vn_kv.size()] etc.

            int nb_bus = bus_vn_kv_.size();
            int nb_line = branch_r.size();

            //diagonal coefficients
            Eigen::VectorXcd diag_vect = Eigen::VectorXcd::Constant(nb_bus, 0.);

            // stuff with half susceptance h, or b or whatever
            // b = line["c_nf_per_km"].values * length_km * parallel *  (2 * net.f_hz * math.pi * 1e-9 * baseR)
            // g = line["g_us_per_km"].values * length_km * parallel * 1e-6 * baseR
            // branch[f:t, BR_B] = b - g * 1j
            Eigen::VectorXcd branch_bb = Eigen::VectorXcd::Constant(nb_line, 2.0 * f_hz_ * M_PI * 1e-9);
            branch_bb.array() *= branch_c.array();
            Eigen::VectorXcd branch_gg = Eigen::VectorXcd::Constant(nb_line, 1e-6);
            branch_gg.array() *= branch_g.array();
//            Eigen::VectorXcd branch_b = branch_cc.array() - 1.0i * branch_gg.array();

            // TODO save g and b here (z = g + j.b) and second "b" too :-) [I hope you realize i'm doing magic here]
            // fill the matrix
            for(int line_id =0; line_id < nb_line; ++line_id){
                // compute from / to
                int from_id = branch_from_id(line_id);
                int to_id = branch_to_id(line_id);

                cdouble base_for_pu_from = bus_pu_(from_id);
                cdouble base_for_pu_to = bus_pu_(to_id);
                cdouble z = branch_r(line_id) + 1.0i * branch_x(line_id);
                cdouble y = 1.0 / z;  // y = 1/((r + 1j * x) / baseR)

                // add the b, why not :-/
                cdouble h = branch_bb(line_id); // yes it's the correct one
                h = static_cast<cdouble>(1.0i) * 0.5 * h;

                //TODO parrallel powerlines!
                // fill non diagonal coefficient
                Ybus_.insert(from_id, to_id) = -y * base_for_pu_from;
                Ybus_.insert(to_id, from_id) = -y * base_for_pu_to;

                // fill diagonal coefficient
                cdouble tmp = y + h;
                diag_vect(from_id) += tmp * base_for_pu_from;
                diag_vect(to_id) += tmp * base_for_pu_to;
            }

            for(int bus_id=0; bus_id < nb_bus; ++bus_id){
                // assign diagonal coefficient
                Ybus_.insert(bus_id, bus_id) = diag_vect(bus_id);
            }
        }

        void init_shunt(const Eigen::VectorXd & shunt_p_mw,
                        const Eigen::VectorXd & shunt_q_mvar,
                        const Eigen::VectorXi & shunt_bus_id)
        {
            /**
            supposes the matrix Ybus_ has already been initialized
            //TODO move that in the Sbus vector instead, more flexible imho
            **/
            int nb_shunt = shunt_q_mvar.size();
            cdouble tmp;
            int bus_id;
            for(int shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
                // assign diagonal coefficient
                tmp = shunt_p_mw(shunt_id) + 1.0i * shunt_q_mvar(shunt_id);
                bus_id = shunt_bus_id(shunt_id);
                Ybus_.coeffRef(bus_id, bus_id) -= tmp;
            }
        }
        void init_trafo(const Eigen::VectorXd & trafo_tap_step_pct,
                        const Eigen::VectorXd & trafo_tap_step_degree,
                        const Eigen::VectorXd & trafo_tap_step_percent,
                        const Eigen::VectorXd & trafo_tap_pos,
                        const Eigen::Vector<bool, Eigen::Dynamic> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                        const Eigen::VectorXd & trafo_vn_hv,
                        const Eigen::VectorXd & trafo_vn_lv,
                        const Eigen::VectorXd & trafo_vk_percent,
                        const Eigen::VectorXd & trafo_vkr_percent,
                        const Eigen::VectorXd & trafo_sn_trafo_mva,
                        const Eigen::VectorXd & trafo_pfe_kw,
                        const Eigen::VectorXd & trafo_i0_pct,
                        const Eigen::VectorXi & trafo_hv_id,
                        const Eigen::VectorXi & trafo_lv_id
                        )
        {
            /**
            This part is purely magical as of now
            **/
            //TODO "parrallel" in the pandapower dataframe, like for lines, are not handled. Handle it python side!

            //TODO only for "trafo model = t"
            //TODO supposes that the step start at 0 for "no ratio"
            int nb_trafo = trafo_tap_step_pct.size();
            Eigen::VectorXd ratio = 1.0 + 0.01 * trafo_tap_step_percent.array() * trafo_tap_pos.array();
            Eigen::VectorXd vn_trafo_lv = trafo_vn_lv;
            Eigen::VectorXd vn_lv = Eigen::VectorXd(nb_trafo);  //bus_vn_kv_[trafo_lv_id]; // TODO check if it compiles
            for(int i = 0; i<nb_trafo; ++i) {vn_lv(i) = bus_vn_kv_(trafo_lv_id(i));}  //TODO optimize that

            // compute r and x
            Eigen::VectorXd tmp = vn_trafo_lv.array() / vn_lv.array();
            tmp = tmp.array() * tmp.array();
            Eigen::VectorXd tap_lv = tmp * sn_mva_;
            Eigen::VectorXd _1_sn_trafo_mva = 1.0 / trafo_sn_trafo_mva.array();
            Eigen::VectorXd z_sc = 0.01 * trafo_vk_percent.array() * _1_sn_trafo_mva.array() * tap_lv.array();
            Eigen::VectorXd r_sc = 0.01 * trafo_vkr_percent.array() * _1_sn_trafo_mva.array() * tap_lv.array();
            Eigen::VectorXd tmp2 = z_sc.array()*z_sc.array() - r_sc.array() * r_sc.array();
            Eigen::VectorXd x_sc = z_sc.cwiseSign().array() * tmp2.cwiseSqrt().array();

            // compute h, the subsceptance
            Eigen::VectorXd baseR = Eigen::VectorXd(nb_trafo);
            for(int i = 0; i<nb_trafo; ++i) {baseR(i) = bus_pu_(trafo_lv_id(i));}  //TODO optimize that
            Eigen::VectorXd pfe =  trafo_pfe_kw.array() * 1e-3;

            // Calculate subsceptance ###
            Eigen::VectorXd vnl_squared = trafo_vn_lv.array() * trafo_vn_lv.array();
            Eigen::VectorXd b_real = pfe.array() / vnl_squared.array() * baseR.array();
            tmp2 = (trafo_i0_pct.array() * 0.01 * trafo_sn_trafo_mva.array());
            Eigen::VectorXd b_img =  tmp2.array() * tmp2.array() - pfe.array() * pfe.array();

            for(int i = 0; i<nb_trafo; ++i) {if (b_img(i) < 0.)  b_img(i) = 0.;}
            b_img = b_img.cwiseSqrt();
            b_img.array() *= baseR.array() / vnl_squared.array();
            Eigen::VectorXcd y = - 1.0i * b_real.array().cast<cdouble>() - b_img.array().cast<cdouble>() * trafo_i0_pct.cwiseSign().array();
            Eigen::VectorXcd h_sc = y.array() / tmp.array();

            //transform trafo from t model to pi model, of course...
            cdouble my_i = 1.0i;
            for(int i = 0; i<nb_trafo; ++i){
                if(h_sc(i) == 0.) continue;
                cdouble za_star = 0.5 * (r_sc(i) +my_i * x_sc(i));
                cdouble zc_star = - my_i / h_sc(i);
                cdouble zSum_triangle = za_star * za_star + 2.0 * za_star * zc_star;
                cdouble zab_triangle = zSum_triangle / zc_star;
                cdouble zbc_triangle = zSum_triangle / za_star;

                r_sc(i) = zab_triangle.real();
                x_sc(i) = zab_triangle.imag();
                h_sc(i) = -2.0 * my_i / zbc_triangle;
            }


            // add it to the admittance matrix
            //TODO

        }

    protected:
        // member of the grid
        double sn_mva_;  // access with net.sn_mva
        double f_hz_;

        // powersystem representation
        // 1. bus
        Eigen::VectorXd bus_vn_kv_;
        Eigen::VectorXd bus_pu_;
        // 2. powerline

        // as matrix, for the solver
        Eigen::SparseMatrix<cdouble> Ybus_;

        // to solve the newton raphson
        KLUSolver _solver;

};

