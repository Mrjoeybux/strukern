#include "edit_costs.hpp"
#include <cmath>
#include <utility>

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

double MoleculeAtomBondTypeSigmoid::edit_param_fun(const double &edit_param) const { return 1.0 / (1.0 + exp(-1 * edit_param)); };

double MoleculeAtomBondTypeSigmoid::edit_param_deriv(const double &edit_param) const {
  double f = this->edit_param_fun(edit_param);
  return f * (1 - f);
};