#ifndef GRAPHBINDING
#define GRAPHBINDING
#include "../domain/graph.h"
#include <pybind11/pybind11.h>
#include <string>
namespace py = pybind11;

template <typename N, typename E> void declare_graph(py::module &m, std::string &typestr1, std::string &typestr2) {
  std::string pyclass_name = std::string("Graph") + typestr1 + typestr2;
  py::class_<Graph<N, E>>(m, pyclass_name.c_str())
      .def(py::init<const unordered_map<uint, vector<uint>> &, const unordered_map<uint, vector<N>> &,
                    const unordered_map<pair<uint, uint>, vector<E>, pair_hash> &, const bool &>(),
           py::arg("AdjacencyList"), py::arg("NodeAttributes") = unordered_map<uint, vector<E>>(),
           py::arg("EdgeAttributes") = unordered_map<pair<uint, uint>, vector<E>, pair_hash>(), py::arg("Directed") = false)

      .def("getID", &Graph<N, E>::getID)
      .def("NumNodes", &Graph<N, E>::NumNodes, "A function which returns the total number of nodes in the graph.")
      .def("NumEdges", &Graph<N, E>::NumEdges, "A function which returns the total number of edges in the graph.")
      .def("getUniqueEdges", &Graph<N, E>::getUniqueEdges)
      .def("getNodeAttribute", &Graph<N, E>::getNodeAttribute, "A function which returns the attribute of a given node.")
      .def("getEdgeAttribute", &Graph<N, E>::getEdgeAttribute, "A function which returns the attribute of a given edge.")
      .def("getNumNeighbours", &Graph<N, E>::getNumNeighbours,
           "A function which returns the total number neighbours of a given "
           "node.")
      .def("getNeighbours", &Graph<N, E>::getNeighbours,
           "A function which returns the indices of the neighbours of a given "
           "node.");
  //.def("__repr__", &Graph<N, E>::AsString);
}
#endif /* GRAPHBINDING */
