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
  py::class_<StringCompressionKernel, AbstractCompressionKernel<string>>(m, "StringCompressionKernel");

  py::class_<ZlibCompressionKernel, StringCompressionKernel>(m, "ZlibCompressionKernel").def(py::init<>());
}

#endif /* STRINGKERNELBINDING */
