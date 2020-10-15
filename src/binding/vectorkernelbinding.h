//
// Created by Mrjoeybux on 13/05/2020.
//

#ifndef STRUKERN_VECTORKERNELBINDING_H
#define STRUKERN_VECTORKERNELBINDING_H
#include "../kernels/basekernel.h"
#include "../kernels/vectorkernels.h"
#include "../kernels/vectorkernels.cpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <string>
namespace py = pybind11;

void declare_vectorkernels(py::module &m) {

  py::class_<RBFKernel, Kernel<VectorXd>>(m, "RBFKernel").def(py::init<>());

  py::class_<Polynomial, Kernel<VectorXd>>(m, "PolynomialKernel").def(py::init<>());
}
#endif //STRUKERN_VECTORKERNELBINDING_H
