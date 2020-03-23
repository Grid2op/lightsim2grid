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

// import klu package
#include "KLUSolver.h"


class DataModel{
    public:
        DataModel():sn_mva_(0.),f_hz_(0.){};
        void set_f_hz(double f_hz) { f_hz_ = f_hz;}
        void set_sn_mva(double sn_mva) { sn_mva_ = sn_mva;}

        void init_bus(const Eigen::VectorXd & bus_vn_kv, int nb_line, int nb_trafo){
            /**
            initialize the bus_vn_kv_ member
            and
            initialize the Ybus_ matrix at the proper shape
            **/
            int nb_bus = bus_vn_kv.size();
            bus_vn_kv_ = bus_vn_kv;  // base_kv
            Ybus_ = Eigen::SparseMatrix<cdouble>(nb_bus, nb_bus);
            Ybus_.reserve(nb_bus + 2*nb_line + 2*nb_trafo);

            // per unit conversion
            bus_pu_ = Eigen::VectorXd::Constant(nb_bus, 1.0 / sn_mva_);
            bus_pu_.array() *= bus_vn_kv_.array() * bus_vn_kv_.array(); // np.square(base_kv) / net.sn_mva

            // init diagonal coefficients
            for(int bus_id=0; bus_id < nb_bus; ++bus_id){
                // assign diagonal coefficient
                Ybus_.insert(bus_id, bus_id) = 0.;
            }
        }

        Eigen::SparseMatrix<cdouble> get_Ybus(){
            return Ybus_;
        }

        void init_powerlines(const Eigen::VectorXd & branch_r,
                             const Eigen::VectorXd & branch_x,
                             const Eigen::VectorXd & branch_c,
//                             const Eigen::VectorXd & branch_g,
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

            int nb_line = branch_r.size();

            powerlines_bus_or_id_ = branch_from_id;
            powerlines_bus_ex_id_ = branch_to_id;

            Eigen::VectorXd bus_pu_from = Eigen::VectorXd(nb_line); // (from_id);
            for(int i = 0; i < nb_line; ++i) {bus_pu_from(i) = bus_pu_(branch_from_id(i));}

            powerlines_r_ = branch_r.array() / bus_pu_from.array();
            powerlines_x_ = branch_x.array() / bus_pu_from.array();

            powerlines_h_ = Eigen::VectorXcd::Constant(nb_line, 2.0 * f_hz_ * M_PI * 1e-9);
            powerlines_h_.array() *= branch_c.array();
            powerlines_h_.array() *=  bus_pu_from.cast<cdouble>().array();

//            fillYbusBranch(powerlines_r_, powerlines_x_, powerlines_h_, powerlines_bus_or_id_, powerlines_bus_ex_id_);

        }

        void init_shunt(const Eigen::VectorXd & shunt_p_mw,
                        const Eigen::VectorXd & shunt_q_mvar,
                        const Eigen::VectorXi & shunt_bus_id)
        {
            /**
            supposes the matrix Ybus_ has already been initialized
            //TODO move that in the Sbus vector instead, more flexible imho
            **/

            shunts_p_mw_ = shunt_p_mw;
            shunts_q_mvar_ = shunt_q_mvar;
            shunts_bus_id_ = shunt_bus_id;
        }

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
                                   const Eigen::VectorXi & trafo_lv_id)
        {
            //TODO only for "trafo model = t"
            //TODO supposes that the step start at 0 for "no ratio"
            int nb_trafo = trafo_vn_lv.size();

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
            Eigen::VectorXcd b_sc = y.array() / tmp.array();

            //transform trafo from t model to pi model, of course...
            // (remove that if trafo model is not t, but directly pi)
            cdouble my_i = 1.0i;
            for(int i = 0; i<nb_trafo; ++i){
                if(b_sc(i) == 0.) continue;
                cdouble za_star = 0.5 * (r_sc(i) + my_i * x_sc(i));
                cdouble zc_star = - my_i / b_sc(i);
                cdouble zSum_triangle = za_star * za_star + 2.0 * za_star * zc_star;
                cdouble zab_triangle = zSum_triangle / zc_star;
                cdouble zbc_triangle = zSum_triangle / za_star;

                r_sc(i) = zab_triangle.real();
                x_sc(i) = zab_triangle.imag();
                b_sc(i) = -2.0 * my_i / zbc_triangle;
            }

            std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXcd> res =
                std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXcd>(std::move(r_sc), std::move(x_sc), std::move(b_sc));
            return res;
        }

        void init_trafo(const Eigen::VectorXd & trafo_r,
                        const Eigen::VectorXd & trafo_x,
                        const Eigen::VectorXcd & trafo_b,
                        const Eigen::VectorXd & trafo_tap_step_pct,
//                        const Eigen::VectorXd & trafo_tap_step_degree,
                        const Eigen::VectorXd & trafo_tap_pos,
                        const Eigen::Vector<bool, Eigen::Dynamic> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                        const Eigen::VectorXi & trafo_hv_id,
                        const Eigen::VectorXi & trafo_lv_id
                        )
        {
            /**
            INPUT DATA ARE ALREADY PAIR UNIT !!
            DOES NOT WORK WITH POWERLINES
            **/
            //TODO "parrallel" in the pandapower dataframe, like for lines, are not handled. Handle it python side!
            Eigen::VectorXd ratio = 1.0 + 0.01 * trafo_tap_step_pct.array() * trafo_tap_pos.array() * (2*trafo_tap_hv.array().cast<double>() - 1.0);

            transformers_r_ = trafo_r;
            transformers_x_ = trafo_x;
            transformers_h_ = trafo_b;
            transformers_ratio_ = ratio;
            transformers_bus_hv_id_ = trafo_hv_id;
            transformers_bus_lv_id_ = trafo_lv_id;
        }

        void compute_newton(){
            init_Ybus();
        };

        void init_Ybus(){
            /**
            Supposes that the powerlines, shunt and transformers are initialized.
            And it fills the Ybus matrix.
            **/
            fillYbusBranch();
            fillYbusShunt();
            fillYbusTrafo();
        }
    protected:
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


        // as matrix, for the solver
        Eigen::SparseMatrix<cdouble> Ybus_;

        // to solve the newton raphson
        KLUSolver _solver;

    protected:

        void fillYbusBranch()
        {
            // fill the matrix
            int nb_line = powerlines_r_.size();
            int nb_bus = bus_vn_kv_.size();
            cdouble my_i = 1.0i;

            //diagonal coefficients
            Eigen::VectorXcd diag_vect = Eigen::VectorXcd::Constant(nb_bus, 0.);

            for(int line_id =0; line_id < nb_line; ++line_id){
                // get the from / to bus id
                int from_id = powerlines_bus_or_id_(line_id);
                int to_id = powerlines_bus_ex_id_(line_id);

                // compute the y
                cdouble z = powerlines_r_(line_id) + my_i * powerlines_x_(line_id);
                cdouble y = 1.0 / z;

                // convert subsceptance to half subsceptance, applied on each ends
                cdouble h = powerlines_h_(line_id); // yes it's the correct one
                h = my_i * 0.5 * h;

                // fill non diagonal coefficient
                Ybus_.coeffRef(from_id, to_id) -= y; // * base_for_pu_from;
                Ybus_.coeffRef(to_id, from_id) -= y; // * base_for_pu_to;

                // fill diagonal coefficient
                cdouble tmp = y + h;
                diag_vect(from_id) += tmp; // * base_for_pu_from;
                diag_vect(to_id) += tmp; //  * base_for_pu_to;
            }

            for(int bus_id=0; bus_id < nb_bus; ++bus_id){
                // assign diagonal coefficient
                Ybus_.coeffRef(bus_id, bus_id) += diag_vect(bus_id);
            }

        }

        void fillYbusShunt(){
            int nb_shunt = shunts_q_mvar_.size();
            cdouble tmp;
            int bus_id;
            for(int shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
                // assign diagonal coefficient
                tmp = shunts_p_mw_(shunt_id) + 1.0i * shunts_q_mvar_(shunt_id);
                bus_id = shunts_bus_id_(shunt_id);
                Ybus_.coeffRef(bus_id, bus_id) -= tmp;
            }
        }

        void fillYbusTrafo(){
            //TODO merge that with fillYbusBranch!
            int nb_trafo = transformers_bus_hv_id_.size();
            cdouble my_i = 1.0i;
            for(int trafo_id =0; trafo_id < nb_trafo; ++trafo_id){
                // compute from / to
                int hv_id = transformers_bus_hv_id_(trafo_id);
                int lv_id = transformers_bus_lv_id_(trafo_id);

                // get the transformers ratio
                double r = transformers_ratio_(trafo_id);

                cdouble z = transformers_r_(trafo_id) + my_i * transformers_x_(trafo_id);
                cdouble y = 1.0 / z;

                // subsecptance
                cdouble h = transformers_h_(trafo_id);
                h = my_i * 0.5 * h;

                // fill non diagonal coefficient
                cdouble tmp = y / r;
                Ybus_.coeffRef(hv_id, lv_id) -= tmp ;
                Ybus_.coeffRef(lv_id, hv_id) -= tmp;

                // fill diagonal coefficient
                Ybus_.coeffRef(hv_id, hv_id) += (tmp + h)  / r ; //TODO is it (tmp / r + h) or (tmp + h) /r ???
                Ybus_.coeffRef(lv_id, lv_id) += (tmp + h) * r;
            }
        }

};

#endif  //DATAMODEL_H
