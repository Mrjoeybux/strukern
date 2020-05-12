#ifndef EDIT_COSTS
#define EDIT_COSTS
#include "../domain/common.h"
#include "src/edit_costs/edit_costs.hpp"
#include "src/env/common_types.hpp"
#include <Dense>
#include <unordered_map>
using namespace ged;
using namespace Eigen;
using namespace std;

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



class MoleculeAtomBondType : public EditCosts<ged::GXLLabel, ged::GXLLabel> {

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

  double edit_param_fun(const double &edit_param) const;

  // virtual double edit_param_deriv(const double &edit_param) const = 0;
};

double MoleculeAtomBondType::edit_param_fun(const double &edit_param) const { return edit_param; };

MoleculeAtomBondType::MoleculeAtomBondType(const LabelPairMap &AtomMap, const LabelPairMap &BondMap) {
  this->AtomMap = AtomMap;
  this->BondMap = BondMap;
}

MoleculeAtomBondType::~MoleculeAtomBondType(){};

double MoleculeAtomBondType::node_ins_cost_fun(const string &node_label) const {
  return this->edit_param_fun(this->AtomMap.at(make_pair(node_label, "INSERT")));
}

double MoleculeAtomBondType::node_del_cost_fun(const string &node_label) const {
  return this->edit_param_fun(this->AtomMap.at(make_pair(node_label, "DELETE")));
}

double MoleculeAtomBondType::node_rel_cost_fun(const string &node_label_1, const string &node_label_2) const {
  if (node_label_1 == node_label_2) {
    return 0.0;
  }
  if (node_label_1 < node_label_2) {
    return this->edit_param_fun(this->AtomMap.at(make_pair(node_label_1, node_label_2)));
  } else {
    return this->edit_param_fun(this->AtomMap.at(make_pair(node_label_1, node_label_2)));
  };
}

double MoleculeAtomBondType::edge_ins_cost_fun(const string &edge_label) const {
  return this->edit_param_fun(this->BondMap.at(make_pair(edge_label, "INSERT")));
}

double MoleculeAtomBondType::edge_del_cost_fun(const string &edge_label) const {
  return this->edit_param_fun(this->BondMap.at(make_pair(edge_label, "DELETE")));
}

double MoleculeAtomBondType::edge_rel_cost_fun(const string &edge_label_1, const string &edge_label_2) const {
  if (edge_label_1 == edge_label_2) {
    return 0.0;
  }
  if (edge_label_1 < edge_label_2) {
    return this->edit_param_fun(this->BondMap.at(make_pair(edge_label_1, edge_label_2)));
  } else {
    return this->edit_param_fun(this->BondMap.at(make_pair(edge_label_1, edge_label_2)));
  };
}

double MoleculeAtomBondType::node_ins_deriv(const string &node_label) const { return 0.0; };

double MoleculeAtomBondType::node_del_deriv(const string &node_label) const { return 0.0; };

double MoleculeAtomBondType::node_rel_deriv(const string &node_label_1, const string &node_label_2) const { return 0.0; };

double MoleculeAtomBondType::edge_ins_deriv(const string &edge_label) const { return 0.0; };

double MoleculeAtomBondType::edge_del_deriv(const string &edge_label) const { return 0.0; };

double MoleculeAtomBondType::edge_rel_deriv(const string &edge_label_1, const string &edge_label_2) const { return 0.0; };
/*class MoleculeAtomBondTypeSigmoid : public MoleculeAtomBondType {

public:
  MoleculeAtomBondTypeSigmoid(const LabelPairMap &AtomMap, const LabelPairMap &BondMap) : MoleculeAtomBondType(AtomMap, BondMap){};

  double edit_param_fun(const double &edit_param) const;

  double edit_param_deriv(const double &edit_param) const;
};*/

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
#endif /* EDIT_COSTS */
