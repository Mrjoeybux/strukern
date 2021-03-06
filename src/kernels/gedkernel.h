#ifndef GEDKERNEL
#define GEDKERNEL

#include "../ged_edit_costs/my_edit_costs.h"
#include "../ged_edit_costs/molecule.h"
#include "basekernel.h"
#include "src/edit_costs/edit_costs.hpp"
#include "src/env/ged_env.hpp"
#include "src/env/ged_graph.hpp"
#include <cstdio>
#include <fstream>
#include <functional>
#include <vector>

using namespace std;
using namespace ged;

struct CollectionResult{
    string Fname;
    vector<int> permute_indicator;
};

enum class GEDEditCosts {
    CatalystEUCLID, CatalystGAUSS, Constant, Dirac, ChemicallyInformed, Mutagenicity
};

enum class GEDMethods {
    BIPARTITE, IPFP
};

class GEDKernel : public Kernel<string> {
private:
    function<EditCosts<ged::GXLLabel, ged::GXLLabel> *(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &,
                                                       const KernelParams &)>
        edit_cost_init;

    function<void(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &, const KernelParams &)> method_init;

    string graph_dir;

    string uniqueID;

    string options;

public:

    GEDKernel(const GEDEditCosts &edit_costs, const GEDMethods &method, const string &graph_directory, const int &num_threads = 1);

    MatrixXd computeKernelMatrix(const vector<string> &data, const KernelParams &params) const;

    MatrixXd computeRectangularKernelMatrix(const vector<string> &data1, const vector<string> &data2,
                                            const KernelParams &params) const;

    // use computeKernelMatrix when applicable
    double dot(const string &x1, const string &x2, const KernelParams &params) const;

    // Adds edit costs and initiliases environment
    EditCosts<ged::GXLLabel, ged::GXLLabel> *init_edit_costs(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                                             const KernelParams &params) const;

    // Adds methods and initialises method.
    void init_methods(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env, const KernelParams &params) const;

    EditCosts<ged::GXLLabel, ged::GXLLabel> *
    init_EditCosts_CatalystEUCLID(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                  const KernelParams &params) const;

    EditCosts<ged::GXLLabel, ged::GXLLabel> *
    init_EditCosts_CatalystGAUSS(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                 const KernelParams &params) const;

    EditCosts<ged::GXLLabel, ged::GXLLabel> *
    init_EditCosts_Constant(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                            const KernelParams &params) const;

    EditCosts<ged::GXLLabel, ged::GXLLabel> *
    init_EditCosts_ChemicallyInformed(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                            const KernelParams &params) const;

    EditCosts<ged::GXLLabel, ged::GXLLabel> *
    init_EditCosts_Mutagenicity(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                      const KernelParams &params) const;

    EditCosts<ged::GXLLabel, ged::GXLLabel> *
    init_EditCosts_Dirac(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                            const KernelParams &params) const;

    void
    init_Method_BIPARTITE(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env, const KernelParams &params) const;

    void
    init_Method_IPFP(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env, const KernelParams &params) const;

    string write_collection_file(const vector<string> &data) const;

    string write_collection_file(const vector<string> &data1, const vector<string> &data2) const;

    void set_ID(const string &ID);
};

GEDKernel::GEDKernel(const GEDEditCosts &edit_costs, const GEDMethods &method, const string &graph_directory, const int &num_threads)
    : Kernel<string>(0.0) {
  switch (edit_costs) {
    case GEDEditCosts::CatalystEUCLID:
      this->edit_cost_init = [this](GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                    const KernelParams &params) {
          return this->init_EditCosts_CatalystEUCLID(env, params);
      };
      break;

    case GEDEditCosts::CatalystGAUSS:
      this->edit_cost_init = [this](GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                    const KernelParams &params) {
          return this->init_EditCosts_CatalystGAUSS(env, params);
      };
      break;

    case GEDEditCosts::ChemicallyInformed:
      this->edit_cost_init = [this](GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                    const KernelParams &params) {
          return this->init_EditCosts_ChemicallyInformed(env, params);
      };
      break;

    case GEDEditCosts::Dirac:
      this->edit_cost_init = [this](GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                    const KernelParams &params) {
          return this->init_EditCosts_Dirac(env, params);
      };
      break;
    case GEDEditCosts::Constant:
      this->edit_cost_init = [this](GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                    const KernelParams &params) {
          return this->init_EditCosts_Constant(env, params);
      };
      break;
    case GEDEditCosts::Mutagenicity:
      this->edit_cost_init = [this](GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                    const KernelParams &params) {
          return this->init_EditCosts_Mutagenicity(env, params);
      };
      break;

    default:
      cout << "Unknown edit costs!" << endl;
  }

  switch (method) {
    case GEDMethods::BIPARTITE:
      this->method_init = [this](GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                 const KernelParams &params) {
          this->init_Method_BIPARTITE(env, params);
      };
      break;

    case GEDMethods::IPFP:
      this->method_init = [this](GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                 const KernelParams &params) {
          this->init_Method_IPFP(env, params);
      };
      break;

    default:
      cout << "Unknown method!" << endl;
  }
  this->options = "--threads " + to_string(num_threads);
  if ('/' != graph_directory.back()) {
    this->graph_dir = graph_directory + "/";
  } else {
    this->graph_dir = graph_directory;
  }
}

EditCosts<ged::GXLLabel, ged::GXLLabel> *
GEDKernel::init_edit_costs(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                           const KernelParams &params) const {
  return this->edit_cost_init(env, params);
}

void
GEDKernel::init_methods(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env, const KernelParams &params) const {
  this->method_init(env, params);
}

double GEDKernel::dot(const string &x1, const string &x2, const KernelParams &params) const {
  // x1 and x2 are the filenames of the respective gxl files containing x1 and x2.
  GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
  vector<string> data = {x1, x2};
  vector<GEDGraph::GraphID> graph_ids;
  // GEDGraph::GraphID id1, id2;
  GEDGraph::GraphID id1 = env.load_gxl_graph(this->graph_dir, x1, Options::GXLNodeEdgeType::LABELED,
                                             Options::GXLNodeEdgeType::LABELED,
                                             std::unordered_set<std::string>(), std::unordered_set<std::string>());
  GEDGraph::GraphID id2 = env.load_gxl_graph(this->graph_dir, x2, Options::GXLNodeEdgeType::LABELED,
                                             Options::GXLNodeEdgeType::LABELED,
                                             std::unordered_set<std::string>(), std::unordered_set<std::string>());
  EditCosts<ged::GXLLabel, ged::GXLLabel> *edit_costs = this->init_edit_costs(env, params);
  this->init_methods(env, params);
  env.run_method(id1, id2);
  delete edit_costs;
  return env.get_upper_bound(id1, id2);
}

MatrixXd GEDKernel::computeKernelMatrix(const vector<string> &data, const KernelParams &params) const {

  /* data contains the filenames of the gxl files containing the graphs used in the computation.*/

  GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
  int n = data.size();
  string collection_file = this->write_collection_file(data);
  vector<GEDGraph::GraphID> graph_ids = env.load_gxl_graphs(this->graph_dir, collection_file);
  MatrixXd K = MatrixXd::Zero(env.num_graphs(), env.num_graphs());
  EditCosts<ged::GXLLabel, ged::GXLLabel> *edit_costs = this->init_edit_costs(env, params);
  this->init_methods(env, params);
  for (uint i = 0; i < env.num_graphs(); i++) {
    for (uint j = i + 1; j < env.num_graphs(); j++) {
      env.run_method(graph_ids[i], graph_ids[j]);
      K(i, j) = env.get_upper_bound(graph_ids[i], graph_ids[j]);
      K(j, i) = K(i, j);
    }
  }
  delete edit_costs;
  remove(collection_file.c_str());
  return K;
};

bool isEmptyString(const string &x) { return x.empty(); };

MatrixXd GEDKernel::computeRectangularKernelMatrix(const vector<string> &data1, const vector<string> &data2,
                                                   const KernelParams &params) const {

  /* data contains the filenames of the gxl files containing the graphs used in the computation.*/
  vector<string> compact_data1(data1), compact_data2(data2);
  compact_data1.erase(remove_if(compact_data1.begin(), compact_data1.end(), isEmptyString), compact_data1.end());
  compact_data2.erase(remove_if(compact_data2.begin(), compact_data2.end(), isEmptyString), compact_data2.end());
  int n = compact_data1.size(), m = compact_data2.size();
  GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
  string collection_file = this->write_collection_file(compact_data1, compact_data2);
  vector<GEDGraph::GraphID> graph_ids = env.load_gxl_graphs(this->graph_dir, collection_file);
  EditCosts<ged::GXLLabel, ged::GXLLabel> *edit_costs = this->init_edit_costs(env, params);
  this->init_methods(env, params);
  MatrixXd K = MatrixXd::Zero(n, m);
  uint i = 0, j = 0;
  for (uint i = 0; i < n; i++) {
    for (uint j = 0; j < m; j++) {
      env.run_method(graph_ids[i], graph_ids[n + j]);
      K(i, j) = env.get_upper_bound(graph_ids[i], graph_ids[n + j]);
    }
  }
  delete edit_costs;
  remove(collection_file.c_str());
  return K;
};

void GEDKernel::init_Method_BIPARTITE(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                      const KernelParams &params) const {
  env.set_method(Options::GEDMethod::BIPARTITE, this->options);
  env.init_method();
}

void GEDKernel::init_Method_IPFP(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                      const KernelParams &params) const {
  env.set_method(Options::GEDMethod::IPFP, this->options);
  env.init_method();
}

string GEDKernel::write_collection_file(const vector<string> &data) const {
  string collection_file = this->graph_dir + this->uniqueID + "_collection_file.xml", graphclass("no_class");
  ofstream file(collection_file);
  file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE GraphCollection SYSTEM>\n<GraphCollection>\n";
  for (uint i = 0; i < data.size(); i++) {
      file << "\t<graph file=\"" << data[i] << "\" class=\"" << graphclass << "\"/>\n";
    }
  file << "</GraphCollection>\n";
  file.close();
  return collection_file;
}

string GEDKernel::write_collection_file(const vector<string> &data1, const vector<string> &data2) const {
  string collection_file = this->graph_dir + this->uniqueID + "_collection_file.xml", graphclass("no_class");

  ofstream file(collection_file);
  file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE GraphCollection SYSTEM>\n<GraphCollection>\n";
  for (auto x : data1) {
    file << "\t<graph file=\"" << x << "\" class=\"" << graphclass << "\"/>\n";
  };
  for (auto x : data2) {
    file << "\t<graph file=\"" << x << "\" class=\"" << graphclass << "\"/>\n";
  };
  file << "</GraphCollection>\n";
  file.close();
  return collection_file;
}

void GEDKernel::set_ID(const string &ID) { this->uniqueID = ID; }

EditCosts<ged::GXLLabel, ged::GXLLabel> *
GEDKernel::init_EditCosts_CatalystEUCLID(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                         const KernelParams &params) const {
  
  MolecularDataset dataset = MolecularDataset::CATALYST;
  Euclidean *edit_costs = new Euclidean(params.MolecularEditCosts.at("base_node_ins"),
                                        params.MolecularEditCosts.at("base_node_del"),
                                        params.MolecularEditCosts.at("base_node_rel"),
                                        params.MolecularEditCosts.at("base_edge_ins"),
                                        params.MolecularEditCosts.at("base_edge_del"),
                                        params.MolecularEditCosts.at("base_edge_rel"),
                                        dataset);
  env.set_edit_costs(edit_costs);
  env.init(Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
  return edit_costs;
}

EditCosts<ged::GXLLabel, ged::GXLLabel> *
GEDKernel::init_EditCosts_CatalystGAUSS(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                        const KernelParams &params) const {
  
  MolecularDataset dataset = MolecularDataset::CATALYST;
  Gaussian *edit_costs = new Gaussian(params.MolecularEditCosts.at("base_node_ins"),
                                      params.MolecularEditCosts.at("base_node_del"),
                                      params.MolecularEditCosts.at("base_node_rel"),
                                      params.MolecularEditCosts.at("base_edge_ins"),
                                      params.MolecularEditCosts.at("base_edge_del"),
                                      params.MolecularEditCosts.at("base_edge_rel"),
                                      dataset);
  env.set_edit_costs(edit_costs);
  env.init(Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
  return edit_costs;
}

EditCosts<ged::GXLLabel, ged::GXLLabel> *
GEDKernel::init_EditCosts_Constant(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                   const KernelParams &params) const {
  env.set_edit_costs(ged::Options::EditCosts::CONSTANT);
  env.init(Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
  return nullptr;
}

EditCosts<ged::GXLLabel, ged::GXLLabel> *
GEDKernel::init_EditCosts_ChemicallyInformed(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                             const KernelParams &params) const {
  MolecularDataset dataset = MolecularDataset::INFORMED;
  
  ChemicallyInformed *edit_costs = new ChemicallyInformed(params.MolecularEditCosts.at("base_node_ins"),
                                                          params.MolecularEditCosts.at("base_node_del"),
                                                          params.MolecularEditCosts.at("base_node_rel"),
                                                          params.MolecularEditCosts.at("base_edge_ins"),
                                                          params.MolecularEditCosts.at("base_edge_del"),
                                                          params.MolecularEditCosts.at("base_edge_rel"),
                                                          dataset);
  env.set_edit_costs(edit_costs);
  env.init(Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
  return edit_costs;
}

EditCosts<ged::GXLLabel, ged::GXLLabel> *
GEDKernel::init_EditCosts_Dirac(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                             const KernelParams &params) const {
  
  MolecularDataset dataset = MolecularDataset::DIRAC;
  Dirac *edit_costs = new Dirac(params.MolecularEditCosts.at("base_node_ins"),
                                params.MolecularEditCosts.at("base_node_del"),
                                params.MolecularEditCosts.at("base_node_rel"),
                                params.MolecularEditCosts.at("base_edge_ins"),
                                params.MolecularEditCosts.at("base_edge_del"),
                                params.MolecularEditCosts.at("base_edge_rel"),
                                dataset);
  env.set_edit_costs(edit_costs);
  env.init(Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
  return edit_costs;
}

EditCosts<ged::GXLLabel, ged::GXLLabel> *
GEDKernel::init_EditCosts_Mutagenicity(GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> &env,
                                       const KernelParams &params) const {
  
  Mutagenicity *edit_costs = new Mutagenicity();
  env.set_edit_costs(edit_costs);
  env.init(Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
  return edit_costs;
}

#endif /* GEDKERNEL */
