#ifndef DATAGEN_H
#define DATAGEN_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "DataGeneric.h"

class DataGen: public DataGeneric
{
    public:
    DataGen() {};

    void init(const Eigen::VectorXd & generators_p,
                     const Eigen::VectorXd & generators_v,
                     const Eigen::VectorXi & generators_bus_id
              );

    int nb() { return p_mw_.size(); }

    void deactivate(int gen_id, bool & need_reset) {_deactivate(gen_id, status_, need_reset);}
    void reactivate(int gen_id, bool & need_reset) {_reactivate(gen_id, status_, need_reset);}
    void change_bus(int gen_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(gen_id, new_bus_id, bus_id_, need_reset, nb_bus);}
    int get_bus(int gen_id) {return _get_bus(gen_id, status_, bus_id_);}
    void change_p(int gen_id, double new_p, bool & need_reset);
    void change_v(int gen_id, double new_v_pu, bool & need_reset);

    void fillSbus(Eigen::VectorXcd & Sbus, bool ac, const std::vector<int> & id_grid_to_solver);
    void fillpv(std::vector<int>& bus_pv,
                std::vector<bool> & has_bus_been_added,
                int slack_bus_id_solver,
                const std::vector<int> & id_grid_to_solver);

    void compute_results(const Eigen::Ref<Eigen::VectorXd> & Va,
                         const Eigen::Ref<Eigen::VectorXd> & Vm,
                         const Eigen::Ref<Eigen::VectorXcd> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const Eigen::VectorXd & bus_vn_kv);
    void reset_results();
    int get_slack_bus_id(int gen_id);
    virtual void set_p_slack(int slack_bus_id, double p_slack);

    void get_vm_for_dc(Eigen::VectorXd & Vm);
    /**
    this functions makes sure that the voltage magnitude of every connected bus is properly used to initialize
    the ac powerflow
    **/
    void set_vm(Eigen::VectorXcd & V, const std::vector<int> & id_grid_to_solver);

    tuple3d get_res() const {return tuple3d(res_p_, res_q_, res_v_);}
    const std::vector<bool>& get_status() const {return status_;}

    protected:
        // physical properties

        // input data
        Eigen::VectorXd p_mw_;
        Eigen::VectorXd vm_pu_;
        Eigen::VectorXi bus_id_;
        std::vector<bool> status_;

        //output data
        Eigen::VectorXd res_p_;  // in MW
        Eigen::VectorXd res_q_;  // in MVar
        Eigen::VectorXd res_v_;  // in kV
};

#endif  //DATAGEN_H
