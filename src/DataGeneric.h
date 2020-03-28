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

        /**
        activation / deactivation of elements
        **/
        void _reactivate(int el_id, std::vector<bool> & status, bool & need_reset);
        void _deactivate(int el_id, std::vector<bool> & status, bool & need_reset);
        void _change_bus(int el_id, int new_bus_me_id, Eigen::VectorXi & el_bus_ids, bool & need_reset);

        /**
        compute the amps from the p, the q and the v (v should NOT be pair unit)
        **/
        void _get_amps(Eigen::VectorXd & a, const Eigen::VectorXd & p, const Eigen::VectorXd & q, const Eigen::VectorXd & v);

};

#endif // DATAGENERIC_H

