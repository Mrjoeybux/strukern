#ifndef STRUKERN
#define STRUKERN
#define GXL_GEDLIB_SHARED

#include "../kernels/basekernel.h"
#include "./graphbinding.h"
#include "./graphkernelbinding.h"
//#include "./imagekernelbinding.h"
#include "./stringkernelbinding.h"
#include "./vectorkernelbinding.h"
#include "./tensorkernelbinding.h"
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Eigen;
namespace py = pybind11;
using namespace std;
using namespace pybind11::literals;

template<typename T>
void declare_abstractbasekernel(py::module &m, string &name) {
  string classname = "Abstract" + name + "Kernel";
  py::class_<Kernel<T>>(m, classname.c_str())
      .def("computeKernelMatrix", &Kernel<T>::computeKernelMatrix)
      .def("computeRectangularKernelMatrix", &Kernel<T>::computeRectangularKernelMatrix)
      .def_static("normaliseHilbert", &Kernel<T>::normaliseHilbert)
      .def_static("normaliseKrein", &Kernel<T>::normaliseKrein)
      .def("dot", &Kernel<T>::dot);
}

template<typename T>
void declare_abstractcompressionkernel(py::module &m, string &name) {
  string classname = "Abstract" + name + "CompressionKernel";
  py::class_<AbstractCompressionKernel<T>, Kernel<T>>(m, classname.c_str())
      .def("compress", &AbstractCompressionKernel<T>::compress)
      .def("dot", &AbstractCompressionKernel<T>::dot)
      .def("concat", &AbstractCompressionKernel<T>::concat)
      .def("computeKernelMatrix", &AbstractCompressionKernel<T>::computeKernelMatrix)
      .def("computeRectangularKernelMatrix", &AbstractCompressionKernel<T>::computeRectangularKernelMatrix);
      //.def("dot", &Kernel<T>::dot);
}

void declare_structs(py::module &m) {
  py::class_<KernelParams>(m, "KernelParams")
      .def(py::init<>())
      .def_readwrite("StringCompressionLevel", &KernelParams::StringCompressionLevel)
      .def_readwrite("MABTS", &KernelParams::MABTS)
      .def_readwrite("JPEGCompressionQuality", &KernelParams::JPEGCompressionQuality)
      .def_readwrite("LocalityImproved", &KernelParams::LocalityImproved)
      .def_readwrite("MolecularEditCosts", &KernelParams::MolecularEditCosts)
      .def_readwrite("RBFSigma", &KernelParams::RBFSigma)
      .def_readwrite("TensorCompressionLevel", &KernelParams::TensorCompressionLevel)
      .def_readwrite("TensorConcatDim", &KernelParams::TensorConcatDim);

  py::enum_<CompressionDistanceMeasure>(m, "CompressionDistanceMeasure", py::arithmetic())
      .value("NCD", CompressionDistanceMeasure::NCD, "Normalised Compression Distance.")
      .value("DIFF", CompressionDistanceMeasure::DIFF, "C(x) + C(y) - C(xy) - C(yx)");

  py::enum_<MultiInstanceMethod>(m, "MultiInstanceMethod", py::arithmetic(),
                                 "Method to compute a multi-instance kernel.")
      .value("Sum", MultiInstanceMethod::SUM, "Sum over all instances in both bags.")
      .value("OA", MultiInstanceMethod::OA, "Optimal assignment of instances in both bags.")
      .value("Max", MultiInstanceMethod::MAX, "Maximal attained value over all instances in both bags.")
      .value("Min", MultiInstanceMethod::MIN, "Minimal attained value over all instances in both bags.");

  /*py::enum_<CompressionMethod>(m, "ImageCompressionMethod", py::arithmetic(),
                               "Method to compute image compression kernel.")
      .value("Vertical", CompressionMethod::Vertical, "Vertical compression kernel.")
      .value("Horizontal", CompressionMethod::Horizontal, "Horizontal compression kernel.")
      .value("Both", CompressionMethod::Both, "Vertical and Horizontal compression kernel.");

  py::enum_<ImageType>(m, "ImageType", py::arithmetic(), "Black and white or RGB images.")
      .value("BW", ImageType::BW, "Black and white images.")
      .value("RGB", ImageType::RGB, "RGB images.");*/

  py::enum_<GEDEditCosts>(m, "GEDEditCosts", py::arithmetic(), "Edit costs used by GED.")
      .value("CatalystEuclid", GEDEditCosts::CatalystEUCLID, "Catalyst features with Euclidean edit costs.")
      .value("CatalystGauss", GEDEditCosts::CatalystGAUSS, "Catalyst features with Gaussian edit costs.")
      .value("Constant", GEDEditCosts::Constant, "Constant edit costs.")
      .value("ChemicallyInformed", GEDEditCosts::Constant, "Chemically informed edit costs.")
      .value("Dirac", GEDEditCosts::Constant, "Dirac edit costs.")
      .value("Mutagenicity", GEDEditCosts::Mutagenicity, "Mutagenicity edit costs.");

  py::enum_<GEDMethods>(m, "GEDMethods", py::arithmetic(), "Method used to compute GED.")
      .value("Bipartite", GEDMethods::BIPARTITE, "BIPARTITE.")
      .value("IPFP", GEDMethods::IPFP, "IPFP");
}

template<typename T>
void declare_multiinstance(py::module &m, string &name) {
  string classname = name + "MultiInstance", abstractclassname = "Abstract" + classname;
  py::class_<Kernel<vector<T>>>(m, abstractclassname.c_str())
      .def("computeKernelMatrix", &Kernel<vector<T>>::computeKernelMatrix)
      .def("computeRectangularKernelMatrix", &Kernel<vector<T>>::computeRectangularKernelMatrix)
      .def_static("normaliseHilbert", &Kernel<vector<T>>::normaliseHilbert)
      .def_static("normaliseKrein", &Kernel<vector<T>>::normaliseKrein)
      .def("dot", &Kernel<vector<T>>::dot);
  py::class_<MultiInstanceKernel<T>, Kernel<vector<T>>>(m, classname.c_str())
      .def(py::init<const Kernel<T> *, const MultiInstanceMethod &>());
}


PYBIND11_MODULE(strukern, m) {
  auto dt = m.def_submodule("datastructures");
  string i = "i", d = "d", s = "s";
  declare_graph<int, int>(dt, i, i);
  declare_graph<double, double>(dt, d, d);
  declare_graph<string, string>(dt, s, s);
  declare_structs(dt);

  auto bk = m.def_submodule("abstractbasekernels");
  string str = "String", img = "Image", vec = "Vector", tens = "Tensor";
  declare_abstractbasekernel<string>(bk, str);
  //declare_abstractbasekernel<ImageMat>(bk, img);
  declare_abstractbasekernel<VectorXd>(bk, vec);
  declare_abstractbasekernel<Tensor>(bk, tens);
  declare_abstractcompressionkernel<string>(bk, str);
  //declare_abstractcompressionkernel<ImageMat>(bk, img);
  declare_abstractcompressionkernel<Tensor>(bk, tens);

  auto vk = m.def_submodule("vectorkernels");
  declare_vectorkernels(vk);
  declare_multiinstance<VectorXd>(vk, vec);

  auto sk = m.def_submodule("stringkernels");
  declare_stringkernels(sk);

  //auto ik = m.def_submodule("imagekernels");
  //declare_imagekernels(ik);

  auto tk = m.def_submodule("tensorkernels");
  declare_tensorkernels(tk);

  auto gk = m.def_submodule("graphkernels");
  string gedname = "GED";
  declare_graphkernels(gk);
  declare_multiinstance<string>(gk, gedname);
};


#endif /* STRUKERN */
