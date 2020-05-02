#ifndef STRUKERN
#define STRUKERN
#define GXL_GEDLIB_SHARED

#include "../kernels/basekernel.h"
#include "./graphbinding.h"
#include "./graphkernelbinding.h"
#include "./imagekernelbinding.h"
#include "./stringkernelbinding.h"
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
            .def_static("normaliseHilbert", &Kernel<T>::normaliseHilbert)
            .def_static("normaliseKrein", &Kernel<T>::normaliseKrein)
            .def("dot", &Kernel<T>::dot);
}

template<typename T>
void declare_abstractcompressionkernel(py::module &m, string &name) {
    string classname = "Abstract" + name + "CompressionKernel";
    py::class_<AbstractCompressionKernel<T>, Kernel<T>>(m, classname.c_str())
            .def("compress", &AbstractCompressionKernel<T>::compress);
}

void declare_structs(py::module &m) {
    py::class_<KernelParams>(m, "KernelParams")
            .def(py::init<>())
            .def_readwrite("ZlibCompressionLevel", &KernelParams::ZlibCompressionLevel)
            .def_readwrite("MABTS", &KernelParams::MABTS)
            .def_readwrite("JPEGCompressionQuality", &KernelParams::JPEGCompressionQuality);

    py::enum_<MultiInstanceMethod>(m, "MultiInstanceMethod", py::arithmetic(),
                                   "Method to compute a multi-instance kernel.")
            .value("Sum", MultiInstanceMethod::SUM, "Sum over all instances.")
            .value("Max", MultiInstanceMethod::MAX, "Calculate optimal assignment over all instances.");

    py::enum_<CompressionMethod>(m, "ImageCompressionMethod", py::arithmetic(),
                                 "Method to compute image compression kernel.")
            .value("Vertical", CompressionMethod::Vertical, "Vertical compression kernel.")
            .value("Horizontal", CompressionMethod::Horizontal, "Horizontal compression kernel.")
            .value("Both", CompressionMethod::Both, "Vertical and Horizontal compression kernel.");

    py::enum_<GEDEditCosts>(m, "GEDEditCosts", py::arithmetic(), "Edit costs used by GED.").value("MABTS",
                                                                                                  GEDEditCosts::MABTS,
                                                                                                  "MABTS.");
    py::enum_<GEDMethods>(m, "GEDMethods", py::arithmetic(), "Method used to compute GED.")
            .value("BIPARTITE", GEDMethods::BIPARTITE, "BIPARTITE.");
}

template<typename T>
void declare_multiinstance(py::module &m, string &name) {
    string classname = name + "MultiInstance", abstractclassname = "Abstract" + classname;
    py::class_<Kernel<vector<T>>>(m, abstractclassname.c_str());
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
    string str = "String", img = "Image";
    declare_abstractbasekernel<string>(bk, str);
    declare_abstractbasekernel<JPEGImageMat>(bk, img);
    declare_abstractcompressionkernel<string>(bk, str);
    declare_abstractcompressionkernel<JPEGImageMat>(bk, img);

    auto sk = m.def_submodule("stringkernels");
    declare_stringkernels(sk);

    auto ik = m.def_submodule("imagekernels");
    declare_imagekernels(ik);

    auto gk = m.def_submodule("graphkernels");
    string gedname = "GED";
    declare_graphkernels(gk);
    declare_multiinstance<string>(gk, gedname);
};


#endif /* STRUKERN */
