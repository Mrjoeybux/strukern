//#include "/home/mrjoeybux/programs/gedlib/src/env/ged_env.hpp"

#ifndef GEDLIB_TEST
#define GEDLIB_TEST
#include "../domain/graph.h"
#include "src/env/ged_env.hpp"
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

void test() {
  ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;

  unordered_map<uint, vector<string>> n_attr;
  n_attr[0] = {"O"};
  n_attr[1] = {"H"};
  n_attr[2] = {"H"};
  unordered_map<pair<uint, uint>, vector<string>, pair_hash> e_attr;
  e_attr[make_pair(0, 1)] = {"bond"};
  e_attr[make_pair(0, 2)] = {"bond"};
  unordered_map<uint, vector<uint>> adjlist;
  adjlist[0] = {1, 2};
  adjlist[1] = {0};
  adjlist[2] = {0};
  bool directed = false;
  Graph<string, string> g0(adjlist, n_attr, e_attr, directed, "h20");

  n_attr[0] = {"O"};
  n_attr[1] = {"H"};
  n_attr[2] = {"H"};
  n_attr[3] = {"H"};

  e_attr[make_pair(0, 1)] = {"bond"};
  e_attr[make_pair(0, 2)] = {"bond"};
  e_attr[make_pair(0, 3)] = {"bond"};

  adjlist[0] = {1, 2, 3};
  adjlist[1] = {0};
  adjlist[2] = {0};
  adjlist[3] = {0};

  vector<string> nodenames = {"Atom"};
  vector<string> edgenames = {"Bond"};

  Graph<string, string> g1(adjlist, n_attr, e_attr, directed, "h30");
  ged::GEDGraph::GraphID id0 = g0.AddToGEDLIB(env, nodenames, edgenames);
  ged::GEDGraph::GraphID id1 = g1.AddToGEDLIB(env, nodenames, edgenames);
  env.save_as_gxl_graph(id0, "h20.gxl");
  env.save_as_gxl_graph(id1, "h30.gxl");

  /*ged::GEDGraph::GraphID id2 = env.add_graph("h30", "class2");
  ged::GXLLabel _20 = {{"Atom", "O"}};
  ged::GXLLabel _21 = {{"Atom", "H"}};
  ged::GXLLabel _22 = {{"Atom", "H"}};
  ged::GXLLabel _23 = {{"Atom", "H"}};
  env.add_node(id2, "0", _20);
  env.add_node(id2, "1", _21);
  env.add_node(id2, "2", _22);
  env.add_node(id2, "3", _23);

  env.add_edge(id2, "0", "1", bond);
  env.add_edge(id2, "0", "2", bond);
  env.add_edge(id2, "0", "3", bond);
  env.save_as_gxl_graph(id1, "h20.gxl");
  env.save_as_gxl_graph(id2, "h30.gxl");

  env.set_edit_costs(ged::Options::EditCosts::CONSTANT);
  env.init();

  env.set_method(ged::Options::GEDMethod::IPFP);
  env.init_method();
  env.run_method(id1, id2);
  double lb, ub;
  lb = env.get_lower_bound(id1, id2);
  ub = env.get_upper_bound(id1, id2);
  cout << "Lower bound: " << lb << "\nUpper bound: " << ub << endl;*/

  // ged::GEDGraph::GraphID id2 = env.add_graph("h30", "class2");
}
#endif /* GEDLIB_TEST */
