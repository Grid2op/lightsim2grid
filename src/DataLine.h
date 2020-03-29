#ifndef DATALINE_H
#define DATALINE_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "DataGeneric.h"

class DataLine : public DataGeneric
{
    public:
    DataLine() {};

    void init(const Eigen::VectorXd & branch_r,
              const Eigen::VectorXd & branch_x,
              const Eigen::VectorXcd & branch_h,
              const Eigen::VectorXi & branch_from_id,
              const Eigen::VectorXi & branch_to_id
              );

    int nb() { return powerlines_r_.size(); }

    void deactivate(int powerline_id, bool & need_reset) {_deactivate(powerline_id, status_, need_reset);}
    void reactivate(int powerline_id, bool & need_reset) {_reactivate(powerline_id, status_, need_reset);}
    void change_bus_or(int powerline_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(powerline_id, new_bus_id, bus_or_id_, need_reset, nb_bus);}
    void change_bus_ex(int powerline_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(powerline_id, new_bus_id, bus_ex_id_, need_reset, nb_bus);}

    void fillYbus(Eigen::SparseMatrix<cdouble> & res, bool ac, const std::vector<int> & id_grid_to_solver);

    void compute_results(const Eigen::Ref<Eigen::VectorXd> & Va,
                         const Eigen::Ref<Eigen::VectorXd> & Vm,
                         const Eigen::Ref<Eigen::VectorXcd> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const Eigen::VectorXd & bus_vn_kv);
    void reset_results();

    tuple4d get_lineor_res() const {return tuple4d(res_powerline_por_, res_powerline_qor_, res_powerline_vor_, res_powerline_aor_);}
    tuple4d get_lineex_res() const {return tuple4d(res_powerline_pex_, res_powerline_qex_, res_powerline_vex_, res_powerline_aex_);}

    protected:
        // physical properties
        Eigen::VectorXd powerlines_r_;
        Eigen::VectorXd powerlines_x_;
        Eigen::VectorXcd powerlines_h_;

        // input data
        Eigen::VectorXi bus_or_id_;
        Eigen::VectorXi bus_ex_id_;
        std::vector<bool> status_;

        //output data
        Eigen::VectorXd res_powerline_por_;  // in MW
        Eigen::VectorXd res_powerline_qor_;  // in MVar
        Eigen::VectorXd res_powerline_vor_;  // in kV
        Eigen::VectorXd res_powerline_aor_;  // in kA
        Eigen::VectorXd res_powerline_pex_;  // in MW
        Eigen::VectorXd res_powerline_qex_;  // in MVar
        Eigen::VectorXd res_powerline_vex_;  // in kV
        Eigen::VectorXd res_powerline_aex_;  // in kA
};

#endif  //DATALINE_H
