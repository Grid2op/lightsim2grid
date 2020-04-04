#ifndef DATASHUNT_H
#define DATASHUNT_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "DataGeneric.h"

class DataShunt : public DataGeneric
{
    public:
    DataShunt() {};

    void init(const Eigen::VectorXd & shunt_p_mw,
                     const Eigen::VectorXd & shunt_q_mvar,
                     const Eigen::VectorXi & shunt_bus_id
              );

    int nb() { return p_mw_.size(); }

    void deactivate(int shunt_id, bool & need_reset) {_deactivate(shunt_id, status_, need_reset);}
    void reactivate(int shunt_id, bool & need_reset) {_reactivate(shunt_id, status_, need_reset);}
    void change_bus(int shunt_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(shunt_id, new_bus_id, bus_id_, need_reset, nb_bus);}
    void change_p(int shunt_id, double new_p, bool & need_reset);
    void change_q(int shunt_id, double new_q, bool & need_reset);
    int get_bus(int shunt_id) {return _get_bus(shunt_id, status_, bus_id_);}

    void fillYbus(std::vector<Eigen::Triplet<cdouble> > & res, bool ac, const std::vector<int> & id_grid_to_solver);
    void fillYbus(Eigen::SparseMatrix<cdouble> & res, bool ac, const std::vector<int> & id_grid_to_solver);

    void compute_results(const Eigen::Ref<Eigen::VectorXd> & Va,
                         const Eigen::Ref<Eigen::VectorXd> & Vm,
                         const Eigen::Ref<Eigen::VectorXcd> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const Eigen::VectorXd & bus_vn_kv);
    void reset_results();
    virtual double get_p_slack(int slack_bus_id);

    tuple3d get_res() const {return tuple3d(res_p_, res_q_, res_v_);}
    const std::vector<bool>& get_status() const {return status_;}

    protected:
        // physical properties

        // input data
        Eigen::VectorXd p_mw_;
        Eigen::VectorXd q_mvar_;
        Eigen::VectorXi bus_id_;
        std::vector<bool> status_;

        //output data
        Eigen::VectorXd res_p_;  // in MW
        Eigen::VectorXd res_q_;  // in MVar
        Eigen::VectorXd res_v_;  // in kV
};

#endif  //DATASHUNT_H
