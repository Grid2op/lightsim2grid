#ifndef DATAGENERIC_H
#define DATAGENERIC_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"

/**
Base class for every object that can be manipulated
**/
class DataGeneric
{
    protected:
        static const int _deactivated_bus_id;
        static const cdouble my_i;

        /**
        activation / deactivation of elements
        **/
        void _reactivate(int el_id, std::vector<bool> & status, bool & need_reset);
        void _deactivate(int el_id, std::vector<bool> & status, bool & need_reset);
        void _change_bus(int el_id, int new_bus_me_id, Eigen::VectorXi & el_bus_ids, bool & need_reset, int nb_bus);

        /**
        compute the amps from the p, the q and the v (v should NOT be pair unit)
        **/
        void _get_amps(Eigen::VectorXd & a, const Eigen::VectorXd & p, const Eigen::VectorXd & q, const Eigen::VectorXd & v);

        /**
        **/
        void v_kv_from_vpu(const Eigen::Ref<Eigen::VectorXd> & Va,
                           const Eigen::Ref<Eigen::VectorXd> & Vm,
                           const std::vector<bool> & status,
                           int nb_element,
                           const Eigen::VectorXi & bus_me_id,
                           const std::vector<int> & id_grid_to_solver,
                           const Eigen::VectorXd & bus_vn_kv,
                           Eigen::VectorXd & v);

};

#endif // DATAGENERIC_H

