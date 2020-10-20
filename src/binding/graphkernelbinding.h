#ifndef GRAPHKERNELBINDING
#define GRAPHKERNELBINDING
#include "../kernels/basekernel.h"
#include "../kernels/gedkernel.h"
#include "../kernels/graphkernels.cpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <string>
namespace py = pybind11;

void declare_graphkernels(py::module &m) {
  py::class_<GEDKernel, Kernel<string>>(m, "GEDKernel").def(py::init<const GEDEditCosts &, const GEDMethods &, const string &, const int &>())
  .def("setID", &GEDKernel::set_ID);
  py::class_<OptimalMolecularAssignmentKernel<vector<double>, vector<double>>, Kernel<UndirectedGraph<vector<double>, vector<double>>>>(m, "OptimalMolecularAssignmentKernel").def(py::init<>());
}
#endif /* GRAPHKERNELBINDING */
