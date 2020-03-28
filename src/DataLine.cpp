#include "DataLine.h"

void DataLine::init(const Eigen::VectorXd & branch_r,
                    const Eigen::VectorXd & branch_x,
                    const Eigen::VectorXcd & branch_h,
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

    //TODO consistency with trafo: have a converter methods to convert this value into pu, and store the pu
    // in this method

    bus_or_id_ = branch_from_id;
    bus_ex_id_ = branch_to_id;
    powerlines_h_ = branch_h;
    powerlines_r_ = branch_r;
    powerlines_x_ = branch_x;
    status_ = std::vector<bool>(branch_r.size(), true); // by default everything is connected
}

void DataLine::fillYbus(Eigen::SparseMatrix<cdouble> & res, bool ac, const std::vector<int> & id_grid_to_solver)
{
    // fill the matrix
    //TODO template here instead of "if" for ac / dc
    int nb_line = powerlines_r_.size();
    cdouble my_i = 1.0i;

    //diagonal coefficients
    for(int line_id =0; line_id < nb_line; ++line_id){
        // i only add this if the powerline is connected
        if(!status_[line_id]) continue;

        // get the from / to bus id
        // compute from / to
        int bus_or_id_me = bus_or_id_(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataLine::fillYbusBranch: A line is connected (or) to a disconnected bus.");
        }
        int bus_ex_id_me = bus_ex_id_(line_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataLine::fillYbusBranch: A line is connected (or) to a disconnected bus.");
        }

        // convert subsceptance to half subsceptance, applied on each ends
        cdouble h = 0.;
        if(ac){
            h = powerlines_h_(line_id); // yes it's the correct one
            h = my_i * 0.5 * h;
        }

        // compute the admittance y
        cdouble y = 0.;
        cdouble z = powerlines_x_(line_id);
        if(ac){
            z *= my_i;
            z += powerlines_r_(line_id);
        }
        if (z !=0. ) y = 1.0 / z;

        // fill non diagonal coefficient
        res.coeffRef(bus_or_solver_id, bus_ex_solver_id) -= y; // * base_for_pu_from;
        res.coeffRef(bus_ex_solver_id, bus_or_solver_id) -= y; // * base_for_pu_to;

        // fill diagonal coefficient
        cdouble tmp = y + h;
        res.coeffRef(bus_or_solver_id, bus_or_solver_id) += tmp;
        res.coeffRef(bus_ex_solver_id, bus_ex_solver_id) += tmp;
    }
}

void DataLine::reset_results()
{
    res_powerline_por_ = Eigen::VectorXd();  // in MW
    res_powerline_qor_ = Eigen::VectorXd();  // in MVar
    res_powerline_vor_ = Eigen::VectorXd();  // in kV
    res_powerline_aor_ = Eigen::VectorXd();  // in kA
    res_powerline_pex_ = Eigen::VectorXd();  // in MW
    res_powerline_qex_ = Eigen::VectorXd();  // in MVar
    res_powerline_vex_ = Eigen::VectorXd();  // in kV
    res_powerline_aex_ = Eigen::VectorXd();  // in kA
}


void DataLine::compute_results(const Eigen::Ref<Eigen::VectorXd> & Va,
                               const Eigen::Ref<Eigen::VectorXd> & Vm,
                               const Eigen::Ref<Eigen::VectorXcd> & V,
                               const std::vector<int> & id_grid_to_solver,
                               const Eigen::VectorXd & bus_vn_kv)
{
    // it needs to be initialized at 0.
    int nb_element = nb();
    res_powerline_por_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MW
    res_powerline_qor_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MVar
    res_powerline_vor_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kV
    res_powerline_aor_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kA
    res_powerline_pex_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MW
    res_powerline_qex_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MVar
    res_powerline_vex_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kV
    res_powerline_aex_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kA
    cdouble my_i = 1.0i;
    for(int line_id = 0; line_id < nb_element; ++line_id){
        // don't do anything if the element is disconnected
        if(!status_[line_id]) continue;

        //physical properties
        double r = powerlines_r_(line_id);
        double x = powerlines_x_(line_id);
        cdouble h = my_i * 0.5 * powerlines_h_(line_id);
        cdouble y = 1.0 / (r + my_i * x);

        // connectivity
        int bus_or_id_me = bus_or_id_(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataModel::res_powerlines: A powerline or a trafo is connected (or) to a disconnected bus.");
        }
        int bus_ex_id_me = bus_ex_id_(line_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataModel::res_powerlines: A powerline or a trafo is connected (ex) to a disconnected bus.");
        }

        // results of the powerflow
        cdouble Eor = V(bus_or_solver_id);
        cdouble Eex = V(bus_ex_solver_id);

        // powerline equations
        cdouble I_orex = (y + h) * Eor - y * Eex;
        cdouble I_exor = (y + h) * Eex - y * Eor;

        I_orex = std::conj(I_orex);
        I_exor = std::conj(I_exor);
        cdouble s_orex = Eor * I_orex;
        cdouble s_exor = Eex * I_exor;

        res_powerline_por_(line_id) = std::real(s_orex);
        res_powerline_qor_(line_id) = std::imag(s_orex);
        res_powerline_pex_(line_id) = std::real(s_exor);
        res_powerline_qex_(line_id) = std::imag(s_exor);

        // retrieve voltages magnitude in kv instead of pu
        double v_or = Vm(bus_or_solver_id);
        double v_ex = Vm(bus_ex_solver_id);
        double bus_vn_kv_or = bus_vn_kv(bus_or_id_me);
        double bus_vn_kv_ex = bus_vn_kv(bus_ex_id_me);
        res_powerline_vor_(line_id) = v_or * bus_vn_kv_or;
        res_powerline_vex_(line_id) = v_ex * bus_vn_kv_ex;
    }
    _get_amps(res_powerline_aor_, res_powerline_por_, res_powerline_qor_, res_powerline_vor_);
    _get_amps(res_powerline_aex_, res_powerline_pex_, res_powerline_qex_, res_powerline_vex_);
}

// for powerline
/**
void DataLine::deactivate(int powerline_id, bool & need_reset)
{
    _deactivate(powerline_id, status_, need_reset);
}
void DataLine::reactivate(int powerline_id, bool & need_reset)
{
    _reactivate(powerline_id, status_, need_reset);
}
void DataLine::change_bus_or(int powerline_id, int new_bus_id, bool & need_reset)
{
    _change_bus(powerline_id, new_bus_id, bus_or_id_, need_reset);
}
void DataLine::change_bus_ex(int powerline_id, int new_bus_id, bool & need_reset)
{
    _change_bus(powerline_id, new_bus_id, bus_ex_id_, need_reset);
}
**/
