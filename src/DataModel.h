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

typedef std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> tuple3d;
typedef std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> tuple4d;

class DataModel{
    public:
        DataModel():sn_mva_(0.),f_hz_(0.){};
        void set_f_hz(double f_hz) { f_hz_ = f_hz;}
        void set_sn_mva(double sn_mva) { sn_mva_ = sn_mva;}

        void init_bus(const Eigen::VectorXd & bus_vn_kv, int nb_line, int nb_trafo);

        // get some internal information
        Eigen::SparseMatrix<cdouble> get_Ybus(){
            return Ybus_;
        }
        Eigen::VectorXcd get_Sbus(){
            return Sbus_;
        }
        Eigen::VectorXi get_pv(){
            return bus_pv_;
        }
        Eigen::VectorXi get_pq(){
            return bus_pq_;
        }
        Eigen::Ref<Eigen::VectorXd> get_Va(){
            return _solver.get_Va();
        }
        Eigen::Ref<Eigen::VectorXd> get_Vm(){
            return _solver.get_Vm();
        }
        // All methods to init this data model
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

        // All results access
        tuple3d get_loads_res() const {return tuple3d(res_load_p_, res_load_q_, res_load_v_);}
        tuple3d get_shunts_res() const {return tuple3d(res_shunt_p_, res_shunt_q_, res_shunt_v_);}
        tuple3d get_gen_res() const {return tuple3d(res_gen_p_, res_gen_q_, res_gen_v_);}
        tuple4d get_lineor_res() const {return tuple4d(res_powerline_por_, res_powerline_qor_, res_powerline_vor_, res_powerline_aor_);}
        tuple4d get_lineex_res() const {return tuple4d(res_powerline_pex_, res_powerline_qex_, res_powerline_vex_, res_powerline_aex_);}
        tuple4d get_trafoor_res() const {return tuple4d(res_trafo_por_, res_trafo_qor_, res_trafo_vor_, res_trafo_aor_);}
        tuple4d get_trafoex_res() const {return tuple4d(res_trafo_pex_, res_trafo_qex_, res_trafo_vex_, res_trafo_aex_);}

        // compute admittance matrix
        void init_Ybus();

        // dc powerflow
        void init_dcY(Eigen::SparseMatrix<double> & dcYbus);
        Eigen::VectorXcd dc_pf(const Eigen::VectorXd & p, const Eigen::VectorXcd Va0);

        // ac powerflows
        bool compute_newton(Eigen::VectorXcd & V,
                            int max_iter,
                            double tol);

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

        // results
        /**
        Compute the results vector from the Va, Vm post powerflow
        **/
        void compute_results();
        /**
        reset the results in case of divergence of the powerflow.
        **/
        void reset_results();

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
        Eigen::VectorXi bus_pv_;
        Eigen::VectorXi bus_pq_;

        // to solve the newton raphson
        Eigen::VectorXcd Sbus_;
        KLUSolver _solver;

        // results of the powerflow
        Eigen::VectorXd res_load_p_; // in MW
        Eigen::VectorXd res_load_q_; // in MVar
        Eigen::VectorXd res_load_v_; // in kV

        Eigen::VectorXd res_gen_p_;  // in MW
        Eigen::VectorXd res_gen_q_;  // in MVar
        Eigen::VectorXd res_gen_v_;  // in kV

        Eigen::VectorXd res_powerline_por_;  // in MW
        Eigen::VectorXd res_powerline_qor_;  // in MVar
        Eigen::VectorXd res_powerline_vor_;  // in kV
        Eigen::VectorXd res_powerline_aor_;  // in kA
        Eigen::VectorXd res_powerline_pex_;  // in MW
        Eigen::VectorXd res_powerline_qex_;  // in MVar
        Eigen::VectorXd res_powerline_vex_;  // in kV
        Eigen::VectorXd res_powerline_aex_;  // in kA

        Eigen::VectorXd res_trafo_por_;  // in MW
        Eigen::VectorXd res_trafo_qor_;  // in MVar
        Eigen::VectorXd res_trafo_vor_;  // in kV
        Eigen::VectorXd res_trafo_aor_;  // in kA
        Eigen::VectorXd res_trafo_pex_;  // in MW
        Eigen::VectorXd res_trafo_qex_;  // in MVar
        Eigen::VectorXd res_trafo_vex_;  // in kV
        Eigen::VectorXd res_trafo_aex_;  // in kA

        Eigen::VectorXd res_shunt_p_;  // in MW
        Eigen::VectorXd res_shunt_q_;  // in MVar
        Eigen::VectorXd res_shunt_v_;  // in kV

    protected:

        void fillYbusBranch(Eigen::SparseMatrix<cdouble> & res, bool ac);
        void fillYbusShunt(Eigen::SparseMatrix<cdouble> & res, bool ac);
        void fillYbusTrafo(Eigen::SparseMatrix<cdouble> & res, bool ac);
        /**
        This method will compute the results for both the powerlines and the trafos
        **/
        void res_powerlines(const Eigen::Ref<Eigen::VectorXd> & Va,
                            const Eigen::Ref<Eigen::VectorXd> & Vm,
                            int nb_element,
                            const Eigen::VectorXd & el_r,
                            const Eigen::VectorXd & el_x,
                            const Eigen::VectorXi & bus_or_id_,
                            const Eigen::VectorXi & bus_ex_id_,
                            Eigen::VectorXd & por,  // in MW
                            Eigen::VectorXd & qor,  // in MVar
                            Eigen::VectorXd & vor,  // in kV
                            Eigen::VectorXd & aor,  // in kA
                            Eigen::VectorXd & pex,  // in MW
                            Eigen::VectorXd & qex,  // in MVar
                            Eigen::VectorXd & vex,  // in kV
                            Eigen::VectorXd & aex  // in kA
                            );

        /**
        compute the amps from the p, the q and the v (v should NOT be pair unit)
        **/
        void _get_amps(Eigen::VectorXd & a, const Eigen::VectorXd & p, const Eigen::VectorXd & q, const Eigen::VectorXd & v);

        /**
        This method will compute the results for the shunt and the loads FOR THE VOLTAGE ONLY
        **/
        void res_loads(const Eigen::Ref<Eigen::VectorXd> & Va,
                       const Eigen::Ref<Eigen::VectorXd> & Vm,
                       int nb_element,
                       const Eigen::VectorXi & bus_id,
                       Eigen::VectorXd & v);

};

#endif  //DATAMODEL_H
