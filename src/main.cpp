#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
//#include <pybind11/stl_bind.h>
//PYBIND11_MAKE_OPAQUE(std::vector<int>);
//PYBIND11_MAKE_OPAQUE(std::vector<double>);
#include <pybind11/eigen.h>

#include "pyklu.h"

namespace py = pybind11;

PYBIND11_MODULE(pyklu_package, m) {
//    py::bind_vector<std::vector<int> >(m, "VectorInt");
//    py::bind_vector<std::vector<double> >(m, "VectorDouble");
    py::class_<KLUSolver>(m, "KLUSolver")
        .def(py::init<>())
        .def("analyze", &KLUSolver::analyze)
        .def("solve", &KLUSolver::solve)
        .def("_evaluate_Fx", &KLUSolver::_evaluate_Fx)
        .def("one_iter", &KLUSolver::one_iter)
        .def("_check_for_convergence", &KLUSolver::_check_for_convergence);
//    py::class_<Dog>(m, "Dog")
//        .def(py::init<>())
//        .def("bark", &Dog::bark);
}