#ifndef GRAPHKERNELBINDING
#define GRAPHKERNELBINDING
#include "../kernels/basekernel.h"
#include "../kernels/gedkernel.h"
#include "../kernels/graphkernels.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <string>
namespace py = pybind11;

void declare_graphkernels(py::module &m) {
  py::class_<GEDKernel, Kernel<string>>(m, "GEDKernel").def(py::init<const GEDEditCosts &, const GEDMethods &, const string &, const int &>())
  .def("setID", &GEDKernel::set_ID);
}
#endif /* GRAPHKERNELBINDING */
