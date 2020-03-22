#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
//#include <pybind11/stl_bind.h>
//PYBIND11_MAKE_OPAQUE(std::vector<int>);
//PYBIND11_MAKE_OPAQUE(std::vector<double>);
#include <pybind11/eigen.h>

#include "pyklu.h"

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

    py::class_<DataModel>(m, "DataModel")
        .def(py::init<>())
        .def("set_f_hz", &DataModel::set_f_hz)
        .def("set_sn_mva", &DataModel::set_sn_mva)
        .def("init_bus", &DataModel::init_bus)
        .def("get_Ybus", &DataModel::get_Ybus)
        .def("init_powerlines", &DataModel::init_powerlines)
        .def("init_shunt", &DataModel::init_shunt)
        .def("init_trafo", &DataModel::init_trafo);

}