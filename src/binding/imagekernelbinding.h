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
  py::class_<ImageCompressionKernel, AbstractCompressionKernel<JPEGImageMat>>(m, "ImageCompressionKernel");

  py::class_<JPEGCompressionKernel, ImageCompressionKernel>(m, "JPEGCompressionKernel").def(py::init<const CompressionMethod &>());
}

#endif /* IMAGEKERNELBINDING */