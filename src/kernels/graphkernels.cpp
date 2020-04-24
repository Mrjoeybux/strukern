#include "graphkernels.h"
#include "../domain/edit_costs.cpp"
#include "src/edit_costs/edit_costs.hpp"
#include <Dense>
#include <cstddef>
using namespace std;
using namespace ged;

template <typename bandwidthType, typename gedNodeLabels, typename gedEdgeLabels>
MatrixXd gedKernel<typename bandwidthType, typename gedNodeLabels, typename gedEdgeLabels>::computeKernelMatrix(
    const vector<string> &data, const unordered_map<string, bandwidthType> &params) const {
  GEDEnv<int, gedNodeLabels, gedEdgeLabels> env;
  vector<GEDGraph::GraphID> graph_ids = env.load_gxl_graphs_from_stream(data);
  EditCosts<gedNodeLabels, gedEdgeLabels> *edit_costs = this->init_edit_costs(env, params);
  this->init_methods(env, params);
  MatrixXd K(env.num_graphs(), env.num_graphs());
  uint i = 0, j;

  for (uint i = 0; i < env.num_graphs(); i++) {
    for (uint j = i; j < env.num_graphs(); j++) {
      env.run_method(graph_ids[i], graph_ids[j]);
      K(i, j) = env.get_upper_bound(graph_ids[i], graph_ids[j]);
      if (i != j) {
        K(j, i) = K(i, j);
      }
    }
  }
  delete edit_costs;
  return K;
};

template <typename bandwidthType, typename gedNodeLabels, typename gedEdgeLabels>
double gedKernel<typename bandwidthType, typename gedNodeLabels, typename gedEdgeLabels>::dot(
    const string &x1, const string &x2, const unordered_map<string, bandwidthType> &params) const {
  GEDEnv<int, gedNodeLabels, gedEdgeLabels> env;
  vector<string> data = {x1, x2};
  vector<GEDGraph::GraphID> graph_ids = env.load_gxl_graphs_from_stream(data);
  EditCosts<gedNodeLabels, gedEdgeLabels> *edit_costs = this->init_edit_costs(env, params);
  this->init_methods(env, params);
  env.run_method(graph_ids[0], graph_ids[1]);
  delete edit_costs;
  return env.get_upper_bound(graph_ids[0], graph_ids[1]);
}

EditCosts<string, string> *gedMABTS::init_edit_costs(GEDEnv<int, string, string> &env,
                                                     const unordered_map<string, LabelPairMap> &params) const {
  MoleculeAtomBondTypeSigmoid *edit_costs = new MoleculeAtomBondTypeSigmoid(params.at("AtomMap"), params.at("BondMap"));
  env.set_edit_costs(edit_costs);
  env.init(Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
  return edit_costs;
}

void gedMABTS_BIPARTITE::init_methods(GEDEnv<int, string, string> &env, const unordered_map<string, LabelPairMap> &params) const {
  env.set_method(Options::GEDMethod::BIPARTITE);
  env.init_methods();
}