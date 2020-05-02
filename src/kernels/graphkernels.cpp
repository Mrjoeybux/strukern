#include "graphkernels.h"
//#include "src/edit_costs/chem_1.hpp"
//#include "src/edit_costs/edit_costs.hpp"
#include <Dense>
using namespace std;
using namespace ged;

GEDKernel::GEDKernel(const GEDEditCosts &edit_costs, const GEDMethods &method) {
  switch (edit_costs) {
  case GEDEditCosts::MABTS:
    /*this->edit_cost_init = [this](GEDEnv<int, string, string> &env, const KernelParams &params) {
      return this->init_EditCosts_MABTS(env, params);
    };*/
    break;

  default:
    cout << "Unknown edit costs!" << endl;
  }

  switch (method) {
    /*case GEDMethods::BIPARTITE:
      this->method_init = [this](GEDEnv<int, string, string> &env, const KernelParams &params) { this->init_Method_BIPARTITE(env, params);
      }; break;*/

  default:
    cout << "Unknown method!" << endl;
  }
}

/*MatrixXd GEDKernel::computeKernelMatrix(const vector<string> &data, const KernelParams &params) const {
  GEDEnv<int, string, string> env;
  vector<GEDGraph::GraphID> graph_ids = env.load_gxl_graphs_from_stream(data);
  // EditCosts<string, string> *edit_costs = this->init_edit_costs(env, params);
  // this->init_methods(env, params);
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
  // delete edit_costs;
  return K;
};*/

double GEDKernel::dot(const string &x1, const string &x2, const KernelParams &params) const {
  GEDEnv<int, string, string> env;
  vector<string> data = {x1, x2};
  vector<GEDGraph::GraphID> graph_ids = env.load_gxl_graphs_from_stream(data);
  // EditCosts<string, string> *edit_costs = this->init_edit_costs(env, params);
  // this->init_methods(env, params);
  env.run_method(graph_ids[0], graph_ids[1]);
  // delete edit_costs;
  return env.get_upper_bound(graph_ids[0], graph_ids[1]);
}

/*EditCosts<string, string> *GEDKernel::init_edit_costs(GEDEnv<int, string, string> &env, const KernelParams &params) const {
  return this->edit_cost_init(env, params);
}*/

// void GEDKernel::init_methods(GEDEnv<int, string, string> &env, const KernelParams &params) const { this->method_init(env, params); }

/*EditCosts<string, string> *GEDKernel::init_EditCosts_MABTS(GEDEnv<int, string, string> &env, const KernelParams &params) const {
  MoleculeAtomBondTypeSigmoid *edit_costs = new MoleculeAtomBondTypeSigmoid(params.MABTS.at("AtomMap"), params.MABTS.at("BondMap"));
  // CHEM1<string, string> *edit_costs = new CHEM1<string, string>();
  env.set_edit_costs(edit_costs);
  env.init(Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
  return edit_costs;
}

void GEDKernel::init_Method_BIPARTITE(GEDEnv<int, string, string> &env, const KernelParams &params) const {
  env.set_method(Options::GEDMethod::BIPARTITE);
  env.init_method();
}*/
