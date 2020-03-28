#include "DataGeneric.h"

const int DataGeneric::_deactivated_bus_id = -1;

// TODO all functions bellow are generic ! Make a base class for that
void DataGeneric::_get_amps(Eigen::VectorXd & a, const Eigen::VectorXd & p, const Eigen::VectorXd & q, const Eigen::VectorXd & v){
    const double _1_sqrt_3 = 1.0 / std::sqrt(3.);
    Eigen::VectorXd p2q2 = p.array() * p.array() + q.array() * q.array();
    p2q2 = p2q2.array().cwiseSqrt();

    // modification in case of disconnected powerlines
    // because i don't want to divide by 0. below
    Eigen::VectorXd v_tmp = v;
    for(auto & el: v_tmp){
        if(el == 0.) el = 1.0;
    }
    a = p2q2.array() * _1_sqrt_3 / v_tmp.array();
}
void DataGeneric::_reactivate(int el_id, std::vector<bool> & status, bool & need_reset){
    bool val = status.at(el_id);
    if(!val) need_reset = true;  // I need to recompute the grid, if a status has changed
    status.at(el_id) = true;  //TODO why it's needed to do that again
}
void DataGeneric::_deactivate(int el_id, std::vector<bool> & status, bool & need_reset){
    bool val = status.at(el_id);
    if(val) need_reset = true;  // I need to recompute the grid, if a status has changed
    status.at(el_id) = false;  //TODO why it's needed to do that again
}
void DataGeneric::_change_bus(int el_id, int new_bus_me_id, Eigen::VectorXi & el_bus_ids, bool & need_reset, int nb_bus){
    // bus id here "me_id" and NOT "solver_id"
    // throw error: object id does not exist
    if(el_id >= el_bus_ids.size()) throw std::out_of_range("change_bus: id too high");
    if(el_id < 0) throw std::out_of_range("change_bus: id too low");
    // throw error: bus id does not exist
    if(new_bus_me_id >= nb_bus) throw std::out_of_range("change_bus: not enough bus to connect object to bus id");
    if(new_bus_me_id < 0) throw std::out_of_range("change_bus: negative bus id");
    int & bus_me_id = el_bus_ids(el_id);
    if(bus_me_id != new_bus_me_id) need_reset = true;  // in this case i changed the bus, i need to recompute the jacobian and reset the solver
    bus_me_id = new_bus_me_id;
}