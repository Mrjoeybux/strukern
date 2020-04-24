#include "src/edit_costs/edit_costs.hpp"
#include "src/env/common_types.hpp"
#include <Dense>
#include <unordered_map>
using namespace ged;
using namespace Eigen;
using namespace std;

struct pair_hash {
  template <typename T1, typename T2> size_t operator()(const pair<T1, T2> &pair) const {
    return hash<T1>()(pair.first) ^ hash<T2>()(pair.second);
  }
};

typedef unordered_map<pair<string, string>, double, pair_hash> LabelPairMap;

/*
Node properties generated from molecules:
atomic_num : uint
is_aromatic : bool
hybridization : enum
chiral_tag : enum
num_radical_electrons : int
formal_charge : uint
*/

/*
Edge properties generated from molecules:
bond_type_as_double : double
is_aromatic : bool
*/
class MoleculeAtomBondType : public EditCosts<string, string> {

private:
  LabelPairMap AtomMap;

  LabelPairMap BondMap;

public:
  MoleculeAtomBondType(const LabelPairMap &AtomMap, const LabelPairMap &BondMap);

  ~MoleculeAtomBondType();

  double node_ins_cost_fun(const string &node_label) const;

  double node_del_cost_fun(const string &node_label) const;

  double node_rel_cost_fun(const string &node_label_1, const string &node_label_2) const;

  double edge_ins_cost_fun(const string &edge_label) const;

  double edge_del_cost_fun(const string &edge_label) const;

  double edge_rel_cost_fun(const string &edge_label_1, const string &edge_label_2) const;

  double node_ins_deriv(const string &node_label) const;

  double node_del_deriv(const string &node_label) const;

  double node_rel_deriv(const string &node_label_1, const string &node_label_2) const;

  double edge_ins_deriv(const string &edge_label) const;

  double edge_del_deriv(const string &edge_label) const;

  double edge_rel_deriv(const string &edge_label_1, const string &edge_label_2) const;

  virtual double edit_param_fun(const double &edit_param) const;

  virtual double edit_param_deriv(const double &edit_param) const;
};

class MoleculeAtomBondTypeSigmoid : public MoleculeAtomBondType {

public:
  MoleculeAtomBondTypeSigmoid(const LabelPairMap &AtomMap, const LabelPairMap &BondMap) : MoleculeAtomBondType(AtomMap, BondMap);

  double edit_param_fun(const double &edit_param) const;

  double edit_param_deriv(const double &edit_param) const;
};
/*
class EuclideanDistance : public EditCosts<GXLLabel, GXLLabel> {
  ~EuclideanDistance();

  EuclideanDistance();

  double node_ins_cost_fun(const GXLLabel &node_label) const;

  double node_del_cost_fun(const GXLLabel &node_label) const;

  double node_rel_cost_fun(const GXLLabel &node_label_1, const GXLLabel &node_label_2) const;

  double edge_ins_cost_fun(const GXLLabel &edge_label) const;

  double edge_del_cost_fun(const GXLLabel &edge_label) const;

  double edge_rel_cost_fun(const GXLLabel &edge_label_1, const GXLLabel &edge_label_2) const;
};

class GaussianDistance : public EditCosts<GXLLabel, GXLLabel> {
private:
  double node_sigma;

  double edge_sigma;

public:
  ~GaussianDistance();

  GaussianDistance(const double &node_sigma, const double &edge_sigma);

  double node_ins_cost_fun(const GXLLabel &node_label) const;

  double node_del_cost_fun(const GXLLabel &node_label) const;

  double node_rel_cost_fun(const GXLLabel &node_label_1, const GXLLabel &node_label_2) const;

  double edge_ins_cost_fun(const GXLLabel &edge_label) const;

  double edge_del_cost_fun(const GXLLabel &edge_label) const;

  double edge_rel_cost_fun(const GXLLabel &edge_label_1, const GXLLabel &edge_label_2) const;
};*/