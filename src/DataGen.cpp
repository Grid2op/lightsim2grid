#include "DataGen.h"

void DataGen::init(const Eigen::VectorXd & generators_p,
                     const Eigen::VectorXd & generators_v,
                     const Eigen::VectorXi & generators_bus_id)
{
    p_mw_ = generators_p;
    vm_pu_ = generators_v;
    bus_id_ = generators_bus_id;
    status_ = std::vector<bool>(generators_p.size(), true);
}

void DataGen::fillSbus(Eigen::VectorXcd & Sbus, bool ac, const std::vector<int> & id_grid_to_solver){
    int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    double tmp;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the load is disconnected
        if(!status_[gen_id]) continue;

        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            //TODO improve error message with the gen_id
            throw std::runtime_error("One generator is connected to a disconnected bus.");
        }
        tmp = p_mw_(gen_id);
        Sbus.coeffRef(bus_id_solver) += tmp;
    }
}

void DataGen::fillpv(std::vector<int> & bus_pv,
                     std::vector<bool> & has_bus_been_added,
                     int slack_bus_id_solver,
                     const std::vector<int> & id_grid_to_solver)
{
    int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;

        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            //TODO improve error message with the gen_id
            throw std::runtime_error("One generator is connected to a disconnected bus.");
        }
        if(bus_id_solver == slack_bus_id_solver) continue;  // slack bus is not PV
        if(has_bus_been_added[bus_id_solver]) continue; // i already added this bus
        bus_pv.push_back(bus_id_solver);
        has_bus_been_added[bus_id_solver] = true;  // don't add it a second time
    }
}

void DataGen::compute_results(const Eigen::Ref<Eigen::VectorXd> & Va,
                               const Eigen::Ref<Eigen::VectorXd> & Vm,
                               const Eigen::Ref<Eigen::VectorXcd> & V,
                               const std::vector<int> & id_grid_to_solver,
                               const Eigen::VectorXd & bus_vn_kv)
{
    int nb_gen = nb();
    v_kv_from_vpu(Va, Vm, status_, nb_gen, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
    res_p_ = p_mw_;
    // res_q_ = q_mvar_;
}

void DataGen::reset_results(){
    res_p_ = Eigen::VectorXd();  // in MW
    res_q_ = Eigen::VectorXd();  // in MVar
    res_v_ = Eigen::VectorXd();  // in kV
}

void DataGen::get_vm_for_dc(Eigen::VectorXd & Vm){
    int nb_gen = nb();
    int bus_id_me;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;
        bus_id_me = bus_id_(gen_id);
        double tmp = vm_pu_(bus_id_me);
        if(tmp != 0.) Vm(bus_id_me) = tmp;
    }
}

void DataGen::change_p(int gen_id, double new_p, bool & need_reset)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if(!my_status) throw std::runtime_error("Impossible to change the active value of a disconnected generator");
    p_mw_(gen_id) = new_p;
}

void DataGen::change_v(int gen_id, double new_v_pu, bool & need_reset)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if(!my_status) throw std::runtime_error("Impossible to change the voltage setpoint of a disconnected generator");
    vm_pu_(gen_id) = new_v_pu;
}

void DataGen::set_vm(Eigen::VectorXcd & V, const std::vector<int> & id_grid_to_solver)
{
    int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;

        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            //TODO improve error message with the gen_id
            throw std::runtime_error("One generator is connected to a disconnected bus.");
        }
        double tmp = std::abs(V(bus_id_solver));
        if(tmp == 0.) tmp = 1.0;
        tmp = 1.0 / tmp;
        tmp *= vm_pu_(gen_id);
        V(bus_id_solver) *= tmp;
    }
}
