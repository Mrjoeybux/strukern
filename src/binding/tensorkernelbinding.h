//
// Created by Mrjoeybux on 23/09/2020.
//
#include "../kernels/tensorkernel.h"
#include "../kernels/tensorkernel.cpp"
#ifndef STRUKERN_TENSORKERNELBINDING_H
#define STRUKERN_TENSORKERNELBINDING_H

void declare_tensorkernels(py::module &m) {
  /*m.def("compressColour", &compressColour);
  m.def("concatVertical", &concatVertical);
    m.def("concatHorizontal", &concatHorizontal);*/

  py::class_<CubicTensorCompressionKernel, AbstractCompressionKernel<Tensor>>(m, "CubicTensorCompressionKernel");

  py::class_<TThreshCompressionKernel, CubicTensorCompressionKernel>(m, "TThreshTensorCompressionKernel")
      .def(py::init<const CompressionDistanceMeasure &>());
}

#endif //STRUKERN_TENSORKERNELBINDING_H
