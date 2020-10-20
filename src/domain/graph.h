#ifndef GRAPH
#define GRAPH
#include <dlib/graph.h>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "common.h"
using namespace dlib;
using namespace std;

template <typename node_type>
struct Node{
    Node(const uint &index, const node_type &attribute){
      this->index = index;
      this->attribute = attribute;
    };

    void add_neighbour(const uint &index){
      this->neighbours.insert(index);
    }
    node_type attribute;
    uint index;
    unordered_set<uint> neighbours;
};

template <typename edge_type>
struct Edge{
public:
    Edge(const uint &start_index, const uint &end_index, const edge_type &attribute){
      this->start_index = start_index;
      this->end_index = end_index;
      this->attribute = attribute;
    };
    edge_type attribute;
    uint start_index;
    uint end_index;
};

template <typename node_type, typename edge_type>
class UndirectedGraph{
public:
    unordered_map<uint, Node<node_type>> nodes;
    unordered_map<pair<uint, uint>, Edge<edge_type>, pair_hash> edges;
    void add_node(const uint &node_index, const node_type &attribute);
    void add_edge(const uint &start_index, const uint &end_index, const edge_type &attribute);
    uint number_of_nodes() const;
    uint number_of_edges() const;
    void add_neighbour(const uint node_index, const uint neighbour_index);
    Node<node_type> get_node(const uint &node_index) const;
    Edge<edge_type> get_edge(const uint &start_index, const uint &end_index) const;
    unordered_set<uint> get_neighbours(const uint &node_index) const;
    uint num_neighbours(const uint &index) const;


};

template <typename node_type, typename edge_type>
uint UndirectedGraph<node_type, edge_type>::number_of_nodes() const{
  return this->nodes.size();
}

template <typename node_type, typename edge_type>
uint UndirectedGraph<node_type, edge_type>::number_of_edges() const{
  return this->edges.size();
}


template <typename node_type, typename edge_type>
void UndirectedGraph<node_type, edge_type>::add_node(const uint &node_index, const node_type &attribute){
  Node<node_type> node_to_add(node_index, attribute);
  this->nodes.insert(make_pair(node_index, node_to_add));
}

template <typename node_type, typename edge_type>
void UndirectedGraph<node_type, edge_type>::add_neighbour(const uint node_index, const uint neighbour_index) {
  this->nodes.at(node_index).add_neighbour(neighbour_index);
}

template <typename node_type, typename edge_type>
void UndirectedGraph<node_type, edge_type>::add_edge(const uint &start_index, const uint &end_index, const edge_type &attribute) {
  Edge<edge_type> edge_to_add(start_index, end_index, attribute);
  this->add_neighbour(start_index, end_index);
  this->add_neighbour(end_index, start_index);
  this->edges.insert(make_pair(make_pair(start_index, end_index), edge_to_add));
}

template <typename node_type, typename edge_type>
Node<node_type> UndirectedGraph<node_type, edge_type>::get_node(const uint &node_index) const{
  return this->nodes.at(node_index);
}

template <typename node_type, typename edge_type>
Edge<edge_type> UndirectedGraph<node_type, edge_type>::get_edge(const uint &start_index, const uint &end_index) const{
  if(this->edges.count(make_pair(start_index, end_index)) > 0){
    return this->edges.at(make_pair(start_index, end_index));
  };
  return this->edges.at(make_pair(end_index, start_index));
}

template <typename node_type, typename edge_type>
unordered_set<uint> UndirectedGraph<node_type, edge_type>::get_neighbours(const uint &node_index) const{
  return this->nodes.at(node_index).neighbours;
}

template <typename node_type, typename edge_type>
uint UndirectedGraph<node_type, edge_type>::num_neighbours(const uint &index) const{
  return this->nodes.at(index).neighbours.size();
}

#endif /* GRAPH */
