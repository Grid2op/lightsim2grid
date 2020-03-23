#ifndef DATAMODEL_H
#define DATAMODEL_H

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
#include "Eigen/SparseLU"
// import klu package
#include "KLUSolver.h"


class DataModel{
    public:
        DataModel():sn_mva_(0.),f_hz_(0.){};
        void set_f_hz(double f_hz) { f_hz_ = f_hz;}
        void set_sn_mva(double sn_mva) { sn_mva_ = sn_mva;}

        void init_bus(const Eigen::VectorXd & bus_vn_kv, int nb_line, int nb_trafo);

        Eigen::SparseMatrix<cdouble> get_Ybus(){
            return Ybus_;
        }
        Eigen::VectorXcd get_Sbus(){
            return Sbus_;
        }

        void init_powerlines(const Eigen::VectorXd & branch_r,
                             const Eigen::VectorXd & branch_x,
                             const Eigen::VectorXd & branch_c,
//                             const Eigen::VectorXd & branch_g,
                             const Eigen::VectorXi & branch_from_id,
                             const Eigen::VectorXi & branch_to_id
                             );
        void init_shunt(const Eigen::VectorXd & shunt_p_mw,
                        const Eigen::VectorXd & shunt_q_mvar,
                        const Eigen::VectorXi & shunt_bus_id);
        void init_trafo(const Eigen::VectorXd & trafo_r,
                        const Eigen::VectorXd & trafo_x,
                        const Eigen::VectorXcd & trafo_b,
                        const Eigen::VectorXd & trafo_tap_step_pct,
//                        const Eigen::VectorXd & trafo_tap_step_degree,
                        const Eigen::VectorXd & trafo_tap_pos,
                        const Eigen::Vector<bool, Eigen::Dynamic> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                        const Eigen::VectorXi & trafo_hv_id,
                        const Eigen::VectorXi & trafo_lv_id
                        );
        void init_generators(const Eigen::VectorXd & generators_p,
                             const Eigen::VectorXd & generators_v,
                             const Eigen::VectorXi & generators_bus_id);
        void init_loads(const Eigen::VectorXd & loads_p,
                        const Eigen::VectorXd & loads_q,
                        const Eigen::VectorXi & loads_bus_id);

        void add_slackbus(int slack_bus_id){
            slack_bus_id_ = slack_bus_id;
        }

        // compute admittance matrix
        void init_Ybus();

        // dc powerflow
        void init_dcY(Eigen::SparseMatrix<double> & dcYbus);
        Eigen::SparseMatrix<double> dc_pf(const Eigen::VectorXd & p);

        // ac powerflows
        void compute_newton();

        // data converters
        std::tuple<Eigen::VectorXd,
           Eigen::VectorXd,
           Eigen::VectorXcd>
           get_trafo_param(const Eigen::VectorXd & trafo_vn_hv,
                           const Eigen::VectorXd & trafo_vn_lv,
                           const Eigen::VectorXd & trafo_vk_percent,
                           const Eigen::VectorXd & trafo_vkr_percent,
                           const Eigen::VectorXd & trafo_sn_trafo_mva,
                           const Eigen::VectorXd & trafo_pfe_kw,
                           const Eigen::VectorXd & trafo_i0_pct,
                           const Eigen::VectorXi & trafo_lv_id);
//    protected:
    // add method to change topology, change ratio of transformers, change

    protected:
        // member of the grid
        double sn_mva_;  // access with net.sn_mva
        double f_hz_;

        // powersystem representation
        // 1. bus
        Eigen::VectorXd bus_vn_kv_;
        Eigen::VectorXd bus_pu_;

        // 2. powerline
        // have the r, x, and h
        Eigen::VectorXd powerlines_r_;
        Eigen::VectorXd powerlines_x_;
        Eigen::VectorXcd powerlines_h_;
        Eigen::VectorXi powerlines_bus_or_id_;
        Eigen::VectorXi powerlines_bus_ex_id_;

        // 3. shunt
        // have the p_mw and q_mvar
        Eigen::VectorXd shunts_p_mw_;
        Eigen::VectorXd shunts_q_mvar_;
        Eigen::VectorXi shunts_bus_id_;

        // 4. transformers
        // have the r, x, h and ratio
        // ratio is computed from the tap, so maybe store tap num and tap_step_pct
        Eigen::VectorXd transformers_r_;
        Eigen::VectorXd transformers_x_;
        Eigen::VectorXcd transformers_h_;
        Eigen::VectorXd transformers_ratio_;
        Eigen::VectorXi transformers_bus_hv_id_;
        Eigen::VectorXi transformers_bus_lv_id_;

        // 5. generators
        Eigen::VectorXd generators_p_;
        Eigen::VectorXd generators_v_;
        Eigen::VectorXi generators_bus_id_;

        // 6. loads
        Eigen::VectorXd loads_p_;
        Eigen::VectorXd loads_q_;
        Eigen::VectorXi loads_bus_id_;

        // 7. slack bus
        int slack_bus_id_;


        // as matrix, for the solver
        Eigen::SparseMatrix<cdouble> Ybus_;

        // to solve the newton raphson
        Eigen::VectorXcd Sbus_;
        KLUSolver _solver;

    protected:

        void fillYbusBranch(Eigen::SparseMatrix<cdouble> & res, bool ac);
        void fillYbusShunt(Eigen::SparseMatrix<cdouble> & res, bool ac);
        void fillYbusTrafo(Eigen::SparseMatrix<cdouble> & res, bool ac);

};

#endif  //DATAMODEL_H
