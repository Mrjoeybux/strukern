#ifndef IMAGEKERNELBINDING
#define IMAGEKERNELBINDING
#include "../kernels/basekernel.h"
#include "../kernels/imagekernels.cpp"
#include "../kernels/imagekernels.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <string>
namespace py = pybind11;

void declare_imagekernels(py::module &m) {
  py::class_<Kernel<JPEGImageMat, int>>(m, "AbstractImageKernel")
      .def("computeKernelMatrix", &Kernel<JPEGImageMat, int>::computeKernelMatrix)
      .def_static("normaliseHilbert", &Kernel<JPEGImageMat, int>::normaliseHilbert)
      .def_static("normaliseKrein", &Kernel<JPEGImageMat, int>::normaliseKrein)
      .def("dot", &Kernel<JPEGImageMat, int>::dot);

  py::class_<AbstractCompressionKernel<JPEGImageMat, int>, Kernel<JPEGImageMat, int>>(m, "AbstractImageCompressionKernel")
      .def("compress", &AbstractCompressionKernel<JPEGImageMat, int>::compress);

  py::class_<vImageCompressionKernel, AbstractCompressionKernel<JPEGImageMat, int>>(m, "ImageCompressionKernel");

  py::class_<vJPEGCompressionKernel, vImageCompressionKernel>(m, "vJPEGCompressionKernel").def(py::init<>());
}

#endif /* IMAGEKERNELBINDING */
