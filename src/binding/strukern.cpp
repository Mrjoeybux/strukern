#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "../kernels/basekernel.h"

using namespace Eigen;
namespace py = pybind11;
using namespace std;
using namespace pybind11::literals;

PYBIND11_MODULE(strukern, m) {
    py::class_<Kernel<string>>(m, "AbstractStringKernel")
        .def_static("normaliseHilbert", &Kernel<string>::normaliseHilbert)
        .def_static("normaliseKrein", &Kernel<string>::normaliseKrein);

    py::class_<DiracKernel<string>, Kernel<string>>(m, "DiracStringKernel")
        .def(py::init<>())
        .def("dot", &DiracKernel<string>::dot)
        .def("computeKernelMatrix", &DiracKernel<string>::computeKernelMatrix);
    
    /*py::class_<Kernel<double>>(m, "AbstractDoubleKernel");
    py::class_<DiracKernel<double>, Kernel<double>>(m, "DiracDoubleKernel")
        .def(py::init<>())
        .def("dot", &DiracKernel<double>::dot);*/
};