#include <vector>
#include "Eigen/Dense"
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

         void analyze(int n,
                      Eigen::Ref<Eigen::VectorXi> Ap,
                      Eigen::Ref<Eigen::VectorXi> Ai){
            n_ = n;
            symbolic_ = klu_analyze(n, &Ap(0), &Ai(0), &common_);
         }

         void solve(Eigen::Ref<Eigen::VectorXi> Ap,
                    Eigen::Ref<Eigen::VectorXi> Ai,
                    Eigen::Ref<Eigen::VectorXd> Ax,
                    Eigen::Ref<Eigen::VectorXd> b){
//            std::cout << b[0] << std::endl;
            numeric_ = klu_factor(&Ap(0), &Ai(0), &Ax(0), symbolic_, &common_);
            klu_solve(symbolic_, numeric_, n_, 1, &b(0), &common_);
//            std::cout << b[0] << std::endl;
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

