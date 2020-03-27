#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
//#include <pybind11/stl_bind.h>
//PYBIND11_MAKE_OPAQUE(std::vector<int>);
//PYBIND11_MAKE_OPAQUE(std::vector<double>);
#include <pybind11/eigen.h>

#include "KLUSolver.h"
#include "DataConverter.h"
#include "DataModel.h"

namespace py = pybind11;

PYBIND11_MODULE(pyklu_cpp, m) {
//    py::bind_vector<std::vector<int> >(m, "VectorInt");
//    py::bind_vector<std::vector<double> >(m, "VectorDouble");
    py::class_<KLUSolver>(m, "KLUSolver")
        .def(py::init<>())
//        .def("initialize_test", &KLUSolver::initialize_test)  // initialize the solver (DO NOT USE)
//        .def("get_ds_test", &KLUSolver::_get_ds_test)  //test function to test if the partial derivatives are properly computed
//        .def("get_Ybus", &KLUSolver::get_Ybus)  //test function to test if the partial derivatives are properly computed
        .def("get_J", &KLUSolver::get_J)  // (get the jacobian matrix, sparse csc matrix)
        .def("get_Va", &KLUSolver::get_Va)  // get the voltage angle vector (vector of double)
        .def("get_Vm", &KLUSolver::get_Vm)  // get the voltage magnitude vector (vector of double)
        .def("get_error", &KLUSolver::get_error)  // get the error message, see the definition of "err_" for more information
        .def("get_nb_iter", &KLUSolver::get_nb_iter)  // return the number of iteration performed at the last optimization
        .def("reset", &KLUSolver::reset)  // reset the solver to its original state
        .def("converged", &KLUSolver::converged)  // whether the solver has converged
        .def("do_newton", &KLUSolver::do_newton, py::call_guard<py::gil_scoped_release>())  // perform the newton raphson optimization
        .def("get_timers", &KLUSolver::get_timers)  // returns the timers corresponding to times the solver spent in different part
        .def("solve", &KLUSolver::do_newton, py::call_guard<py::gil_scoped_release>() );  // perform the newton raphson optimization


    // converters
    py::class_<PandaPowerConverter>(m, "PandaPowerConverter")
        .def(py::init<>())
        .def("set_f_hz", &PandaPowerConverter::set_f_hz)
        .def("set_sn_mva", &PandaPowerConverter::set_sn_mva)
        .def("get_line_param", &PandaPowerConverter::get_line_param)
        .def("get_trafo_param", &PandaPowerConverter::get_trafo_param);

    py::class_<DataModel>(m, "DataModel")
        .def(py::init<>())
        // general parameters

        // init the grid
        .def("init_bus", &DataModel::init_bus)
        .def("init_powerlines", &DataModel::init_powerlines)
        .def("init_shunt", &DataModel::init_shunt)
        .def("init_trafo", &DataModel::init_trafo)
        .def("init_generators", &DataModel::init_generators)
        .def("init_loads", &DataModel::init_loads)
        .def("add_slackbus", &DataModel::add_slackbus)

        // modify the grid
        .def("deactivate_bus", &DataModel::deactivate_bus)
        .def("reactivate_bus", &DataModel::reactivate_bus)

        .def("deactivate_powerline", &DataModel::deactivate_powerline)
        .def("reactivate_powerline", &DataModel::reactivate_powerline)
        .def("change_bus_powerline_or", &DataModel::change_bus_powerline_or)
        .def("change_bus_powerline_ex", &DataModel::change_bus_powerline_ex)

        .def("deactivate_trafo", &DataModel::deactivate_trafo)
        .def("reactivate_trafo", &DataModel::reactivate_trafo)
        .def("change_bus_trafo_hv", &DataModel::change_bus_trafo_hv)
        .def("change_bus_trafo_lv", &DataModel::change_bus_trafo_lv)

        .def("deactivate_load", &DataModel::deactivate_load)
        .def("reactivate_load", &DataModel::reactivate_load)
        .def("change_bus_load", &DataModel::change_bus_load)

        .def("deactivate_gen", &DataModel::deactivate_gen)
        .def("reactivate_gen", &DataModel::reactivate_gen)
        .def("change_bus_gen", &DataModel::change_bus_gen)

        .def("deactivate_shunt", &DataModel::deactivate_shunt)
        .def("reactivate_shunt", &DataModel::reactivate_shunt)
        .def("change_bus_shunt", &DataModel::change_bus_shunt)

        // get back the results
        .def("get_Va", &DataModel::get_Va)
        .def("get_Vm", &DataModel::get_Vm)
        .def("get_loads_res", &DataModel::get_loads_res)
        .def("get_shunts_res", &DataModel::get_shunts_res)
        .def("get_gen_res", &DataModel::get_gen_res)
        .def("get_lineor_res", &DataModel::get_lineor_res)
        .def("get_lineex_res", &DataModel::get_lineex_res)
        .def("get_trafohv_res", &DataModel::get_trafohv_res)
        .def("get_trafolv_res", &DataModel::get_trafolv_res)

        // do something with the grid
        // .def("init_Ybus", &DataModel::init_Ybus) // temporary
        .def("get_Ybus", &DataModel::get_Ybus)
        .def("get_Sbus", &DataModel::get_Sbus)
        .def("get_pv", &DataModel::get_pv)
        .def("get_pq", &DataModel::get_pq)
        .def("dc_pf", &DataModel::dc_pf)
        .def("compute_newton", &DataModel::compute_newton)
        ;

}