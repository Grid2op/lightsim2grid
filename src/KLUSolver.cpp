#include "KLUSolver.h"

bool KLUSolver::do_newton(const Eigen::SparseMatrix<cdouble> & Ybus,
                          Eigen::VectorXcd & V,
                          const Eigen::VectorXcd & Sbus,
                          const Eigen::VectorXi & pv,
                          const Eigen::VectorXi & pq,
                          int max_iter,
                          double tol
                          )
{
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
    V_ = V;

    // first check, if the problem is already solved, i stop there
    Eigen::VectorXd F = _evaluate_Fx(Ybus, V, Sbus, pv, pq);
    bool converged = _check_for_convergence(F, tol);
    nr_iter_ = 0; //current step
    bool res = true;  // have i converged or not
    bool has_just_been_inialized = false;  // to avoid a call to klu_refactor follow a call to klu_factor in the same loop
    while ((!converged) & (nr_iter_ < max_iter)){
        nr_iter_++;
        fill_jacobian_matrix(Ybus, V_, pq, pvpq, pq_inv, pvpq_inv);
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
        Vm_ = V_.array().abs();  // update Vm and Va again in case
        Va_ = V_.array().arg();  // we wrapped around with a negative Vm

        if (n_pv > 0) Va_(pv) += dx.segment(0,n_pv);
        if (n_pq > 0){
            Va_(pq) += dx.segment(n_pv,n_pq);
            Vm_(pq) += dx.segment(n_pv+n_pq, n_pq);
        }

        // TODO change here for not having to cast all the time ... maybe
        V_ = Vm_.array() * (Va_.array().cos().cast<cdouble>() + 1.0i * Va_.array().sin().cast<cdouble>() );

        F = _evaluate_Fx(Ybus, V_, Sbus, pv, pq);
        converged = _check_for_convergence(F, tol);
    }
    if(!converged){
        err_ = 4;
        res = false;
    }
    timer_total_nr_ += timer.duration();
    return res;
}

void KLUSolver::reset(){
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

void KLUSolver::initialize(){
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

void KLUSolver::solve(Eigen::VectorXd & b, bool has_just_been_inialized){
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

void KLUSolver::_dSbus_dV(const Eigen::Ref<const Eigen::SparseMatrix<cdouble> > & Ybus,
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

void KLUSolver::_get_values_J(int & nb_obj_this_col,
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

void KLUSolver::fill_jacobian_matrix(const Eigen::SparseMatrix<cdouble> & Ybus,
                                     const Eigen::VectorXcd & V,
                                     const Eigen::VectorXi & pq,
                                     const Eigen::VectorXi & pvpq,
                                     const std::vector<int> & pq_inv,
                                     const std::vector<int> & pvpq_inv
                                     )
{
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
    _dSbus_dV(Ybus, V);
    Eigen::SparseMatrix<double> dS_dVa_r = dS_dVa_.real();
    Eigen::SparseMatrix<double> dS_dVa_i = dS_dVa_.imag();
    Eigen::SparseMatrix<double> dS_dVm_r = dS_dVm_.real();
    Eigen::SparseMatrix<double> dS_dVm_i = dS_dVm_.imag();

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

Eigen::VectorXd KLUSolver::_evaluate_Fx(const Eigen::SparseMatrix<cdouble> &  Ybus,
                                        const Eigen::VectorXcd & V,
                                        const Eigen::VectorXcd & Sbus,
                                        const Eigen::VectorXi & pv,
                                        const Eigen::VectorXi & pq)
{
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