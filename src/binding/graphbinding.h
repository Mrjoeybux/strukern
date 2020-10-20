#ifndef GRAPHBINDING
#define GRAPHBINDING
#include "../domain/graph.h"
#include <pybind11/pybind11.h>
#include <string>
namespace py = pybind11;

template <typename node_type, typename edge_type> void declare_graph(py::module &m, std::string &typestr1, std::string &typestr2) {
  std::string pyclass_name = std::string("Graph") + typestr1 + typestr2;
  py::class_<UndirectedGraph<node_type, edge_type>>(m, pyclass_name.c_str())
      .def(py::init<>())
      .def("add_node", &UndirectedGraph<node_type, edge_type>::add_node)
      .def("add_edge", &UndirectedGraph<node_type, edge_type>::add_edge)
      .def("get_neighbours", &UndirectedGraph<node_type, edge_type>::get_neighbours);
  //.def("__repr__", &Graph<N, E>::AsString);
}
#endif /* GRAPHBINDING */
