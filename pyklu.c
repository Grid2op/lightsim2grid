/* a simple wrapper for easy use in PyKLU's python source code */

#include "cs.h"
#include "klu.h"


int solve_linear_system(int n, int* Ap, int* Ai, double* Ax, double* b) {
    //this only at the beginning
    klu_symbolic* Symbolic;
    klu_numeric* Numeric;
    klu_common Common;
    klu_defaults(&Common);
    Symbolic = klu_analyze(n, Ap, Ai, &Common);

    // I need to do that at each iteration
    Numeric = klu_factor(Ap, Ai, Ax, Symbolic, &Common);
    klu_solve(Symbolic, Numeric, n, 1, b, &Common);

    // this only at the end
    klu_free_symbolic(&Symbolic, &Common);
    klu_free_numeric(&Numeric, &Common);
    return 0;
}
//void do_newton(double *Ybus, double *V,
//               int pvpq, //TODO
//               int pq, //TODO
//               int createJ, //TODO
//               int pvpq_lookup, //TODO
//               int npv,  // number of pv nodes
//               int npq //number of pq nodes
//               ){
//    int j1 = 0;
//    int j2 = npv ; // j1:j2 - V angle of pv buses
//    int j3 = j2;
//    int j4 = j2 + npq;  // j3:j4 - V angle of pq buses
//    int j5 = j4;
//    int j6 = j4 + npq;  // j5:j6 - V mag of pq buses
//
//}
//
//void _evaluate_Fx(cs* F, cs* Ybus, cs *V, double *Sbus, int pv, int pq){
//    // evalute F(x)
//    F = cs_multiply(Ybus, V);
//    // mis = V * conj(Ybus * V) - Sbus
//    // F = r_[mis[pv].real,
//    //       mis[pq].real,
//    //       mis[pq].imag]
//}