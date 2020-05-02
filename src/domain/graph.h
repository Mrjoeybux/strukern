#ifndef GRAPH
#define GRAPH
#include "common.h"
#include "src/env/ged_env.hpp"
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

template <typename N = int, typename E = int> class Graph {
private:
  unordered_map<uint, vector<uint>> adjlist;
  unordered_map<uint, vector<N>> n_attr;
  unordered_map<pair<uint, uint>, vector<E>, pair_hash> e_attr;
  bool directed;
  string ID;

public:
  Graph(
      const unordered_map<uint, vector<uint>> &AdjacencyList,
      const unordered_map<uint, vector<N>> &NodeAttributes = unordered_map<uint, vector<N>>(),
      const unordered_map<pair<uint, uint>, vector<E>, pair_hash> &EdgeAttributes = unordered_map<pair<uint, uint>, vector<E>, pair_hash>(),
      const bool &Directed = false, const string &ID = "");

  Graph();

  ~Graph(){};

  string getID() const;

  uint NumNodes() const;

  uint NumEdges() const;

  vector<pair<uint, uint>> getUniqueEdges() const;

  vector<uint> getNeighbours(const uint &index) const;

  uint getNumNeighbours(const uint &index) const;

  optional<vector<N>> getNodeAttribute(const uint &index) const;

  optional<vector<E>> getEdgeAttribute(const uint &index1, const uint &index2) const;

  // string AsString() const;

  ged::GEDGraph::GraphID AddToGEDLIB(ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                     const vector<string> &NodeAttributeNames, const vector<string> &EdgeAttributeNames) const;

  ged::GXLLabel createNodeLabel(const uint &index, const vector<string> &NodeAttributeNames) const;

  ged::GXLLabel createEdgeLabel(const uint &index1, const uint &index2, const vector<string> &EdgeAttributeNames) const;
};

template <typename N, typename E> Graph<N, E>::Graph(){};

template <typename N, typename E>
Graph<N, E>::Graph(const unordered_map<uint, vector<uint>> &AdjacencyList, const unordered_map<uint, vector<N>> &NodeAttributes,
                   const unordered_map<pair<uint, uint>, vector<E>, pair_hash> &EdgeAttributes, const bool &Directed, const string &ID) {

  this->directed = Directed;
  this->adjlist = AdjacencyList;
  this->n_attr = NodeAttributes;
  this->e_attr = EdgeAttributes;
  this->ID = ID;
};

template <typename N, typename E> string Graph<N, E>::getID() const { return this->ID; }

template <typename N, typename E> uint Graph<N, E>::NumNodes() const { return this->adjlist.size(); }

template <typename N, typename E> uint Graph<N, E>::NumEdges() const {
  uint num_edges = 0;
  for (int i = 0; i < this->NumNodes(); i++) {
    num_edges += this->adjlist.at(i).size();
  }
  if (this->directed) {
    return num_edges;
  }
  return num_edges / 2;
}

template <typename N, typename E> vector<pair<uint, uint>> Graph<N, E>::getUniqueEdges() const {
  vector<pair<uint, uint>> edges;
  for (const auto &node_id : this->adjlist) {
    for (const auto &neigh_id : node_id.second) {
      pair<uint, uint> forward_pair = make_pair(node_id.first, neigh_id), backward_pair = make_pair(neigh_id, node_id.first);
      if ((find(edges.begin(), edges.end(), forward_pair) == edges.end()) &&
          (find(edges.begin(), edges.end(), backward_pair) == edges.end())) {
        edges.emplace_back(forward_pair);
      } else {
        continue;
      }
    }
  }
  return edges;
}

template <typename N, typename E> vector<uint> Graph<N, E>::getNeighbours(const uint &index) const { return this->adjlist.at(index); }

template <typename N, typename E> uint Graph<N, E>::getNumNeighbours(const uint &index) const { return this->adjlist.at(index).size(); }

template <typename N, typename E> optional<vector<N>> Graph<N, E>::getNodeAttribute(const uint &index) const {
  if (this->n_attr.empty()) {
    return {};
  }
  return this->n_attr.at(index);
}

template <typename N, typename E> optional<vector<E>> Graph<N, E>::getEdgeAttribute(const uint &index1, const uint &index2) const {
  if (this->e_attr.empty()) {
    return {};
  }
  if (this->e_attr.find(pair<uint, uint>(index1, index2)) != this->e_attr.end()) {
    return this->e_attr.at(pair<uint, uint>(index1, index2));
  }
  if (this->e_attr.find(pair<uint, uint>(index2, index1)) != this->e_attr.end()) {
    return this->e_attr.at(pair<uint, uint>(index2, index1));
  }
  return {};
}

template <typename N, typename E>
ged::GXLLabel Graph<N, E>::createNodeLabel(const uint &index, const vector<string> &NodeAttributeNames) const {
  ged::GXLLabel label;
  optional<vector<E>> attr = this->getNodeAttribute(index);
  if (attr.has_value()) {
    for (uint i = 0; i < NodeAttributeNames.size(); i++) {
      label[NodeAttributeNames[i]] = attr.value()[i];
    }
  }
  return label;
}

template <typename N, typename E>
ged::GXLLabel Graph<N, E>::createEdgeLabel(const uint &index1, const uint &index2, const vector<string> &EdgeAttributeNames) const {
  ged::GXLLabel label;
  optional<vector<E>> attr = this->getEdgeAttribute(index1, index2);
  if (attr.has_value()) {
    for (uint i = 0; i < EdgeAttributeNames.size(); i++) {
      label[EdgeAttributeNames[i]] = attr.value()[i];
    }
  }
  return label;
}

template <typename N, typename E>
ged::GEDGraph::GraphID Graph<N, E>::AddToGEDLIB(ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                                const vector<string> &NodeAttributeNames, const vector<string> &EdgeAttributeNames) const {
  ged::GEDGraph::GraphID id = env.add_graph(this->getID(), "");
  for (uint i = 0; i < this->NumNodes(); i++) {
    ged::GXLLabel label = this->createNodeLabel(i, NodeAttributeNames);
    env.add_node(id, to_string(i), label);
  }
  vector<pair<uint, uint>> unique_edges = this->getUniqueEdges();
  for (uint i = 0; i < unique_edges.size(); i++) {
    ged::GXLLabel label = this->createEdgeLabel(unique_edges[i].first, unique_edges[i].second, EdgeAttributeNames);
    env.add_edge(id, to_string(unique_edges[i].first), to_string(unique_edges[i].second), label);
  }
  return id;
}

template <typename N, typename E>
void save_graphs_as_GXL(const vector<Graph<N, E>> &graphs, const vector<string> &NodeNames, const vector<string> &EdgeNames,
                        const vector<string> &FileNames = {}, const string &save_dir = "") {
  uint n = graphs.size();
  string name, filename = save_dir;
  ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
  for (uint i = 0; i < n; i++) {
    if (FileNames.empty()) {
      filename += "graph_" + to_string(i) + ".gxl";
    } else {
      filename += FileNames[i];
    }
    ged::GEDGraph::GraphID id = graphs[i].AddToGEDLIB(env, NodeNames, EdgeNames);
    env.save_as_gxl_graph(id, filename);
  }
}

/*template <typename N, typename E> string Graph<N, E>::AsString() const {
  stringstream ss;
  ss << "Graph object with adjacency list as\n";
  for (uint i = 0; i < this->NumNodes(); i++) {
    ss << "Node " << i << ": [";
    uint n = this->getNumNeighbours(i), count = 1;
    for (auto neigh : this->getNeighbours(i)) {
      if (count == n) {
        ss << neigh << "]\n";
      } else {
        ss << neigh << ", ";
      }
      count++;
    }
  }
  return ss.str();
}*/
#endif /* GRAPH */
