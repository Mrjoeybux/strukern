#ifndef STRINGKERNELBINDING
#define STRINGKERNELBINDING
#include "../kernels/basekernel.h"
#include "../kernels/stringkernels.cpp"
#include "../kernels/stringkernels.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <string>
namespace py = pybind11;

void declare_stringkernels(py::module &m) {
  py::class_<Kernel<string, int>>(m, "AbstractStringKernel")
      .def("computeKernelMatrix", &Kernel<string, int>::computeKernelMatrix)
      .def_static("normaliseHilbert", &Kernel<string, int>::normaliseHilbert)
      .def_static("normaliseKrein", &Kernel<string, int>::normaliseKrein)
      .def("dot", &Kernel<string, int>::dot);

  py::class_<AbstractCompressionKernel<string, int>, Kernel<string, int>>(m, "AbstractCompressionKernel")
      .def("compress", &AbstractCompressionKernel<string, int>::compress);

  py::class_<StringCompressionKernel, AbstractCompressionKernel<string, int>>(m, "StringCompressionKernel");

  py::class_<ZlibCompressionKernel, StringCompressionKernel>(m, "ZlibCompressionKernel").def(py::init<>());
}

#endif /* STRINGKERNELBINDING */
