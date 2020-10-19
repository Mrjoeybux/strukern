#ifndef GRAPH
#define GRAPH
#include <dlib/graph.h>
#include <vector>
using namespace dlib;
using namespace std;

template <typename node_type, typename edge_type>
void add_undirected_graph(const vector<node_type> &node_labels, const vector<pair<uint, uint>> &edge_connectivity, const vector<edge_type> &edge_labels);


#endif /* GRAPH */
