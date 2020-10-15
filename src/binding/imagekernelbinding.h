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
  /*m.def("compressColour", &compressColour);
  m.def("concatVertical", &concatVertical);
    m.def("concatHorizontal", &concatHorizontal);*/

  py::class_<ImageCompressionKernel, AbstractCompressionKernel<ImageMat>>(m, "ImageCompressionKernel");

  py::class_<JPEGCompressionKernel, ImageCompressionKernel>(m, "JPEGCompressionKernel").def(py::init<const CompressionMethod &, const ImageType &>());
}

#endif /* IMAGEKERNELBINDING */
