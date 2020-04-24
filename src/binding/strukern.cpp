#include "./graphbinding.h"
#include "./imagekernelbinding.h"
#include "./stringkernelbinding.h"
//#include <iostream>
//
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Eigen;
namespace py = pybind11;
using namespace std;
using namespace pybind11::literals;

PYBIND11_MODULE(strukern, m) {

  string i = "i", d = "d", s = "s";
  declare_graph<int, int>(m, i, i);
  declare_graph<double, double>(m, d, d);
  declare_graph<string, string>(m, s, s);

  // auto gk = m.def_submodule("graphkernels");
  auto sk = m.def_submodule("stringkernels");
  auto ik = m.def_submodule("imagekernels");
  // auto ik = m.def_submodule("imagekernels");
  declare_stringkernels(sk);
  declare_imagekernels(ik);
};