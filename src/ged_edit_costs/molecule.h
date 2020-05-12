//
// Created by Mrjoeybux on 08/05/20.
//
#include "src/edit_costs/edit_costs.hpp"
#include "src/env/common_types.hpp"

using namespace ged;
using namespace std;

#ifndef STRUKERN_MOLECULE_H
#define STRUKERN_MOLECULE_H

enum class MolecularDataset {
    AHNEMAN, CATALYST
};

class Molecule : public EditCosts<ged::GXLLabel, ged::GXLLabel> {
private:
    double base_node_ins;
    double base_node_del;
    double base_node_rel;
    double base_edge_ins;
    double base_edge_del;
    double base_edge_rel;
    function<vector<double>(const ged::GXLLabel &)> node_label_extractor;

    function<vector<double>(const ged::GXLLabel &)> edge_label_extractor;

public:
    Molecule(const double &base_node_ins, const double &base_node_del, const double &base_node_rel,
             const double &base_edge_ins, const double &base_edge_del, const double &base_edge_rel,
             const MolecularDataset &dataset);

    double node_ins_cost_fun(const ged::GXLLabel &node_label) const;

    double node_del_cost_fun(const ged::GXLLabel &node_label) const;

    double node_rel_cost_fun(const ged::GXLLabel &node_label_1, const ged::GXLLabel &node_label_2) const;

    double edge_ins_cost_fun(const ged::GXLLabel &edge_label) const;

    double edge_del_cost_fun(const ged::GXLLabel &edge_label) const;

    double edge_rel_cost_fun(const ged::GXLLabel &edge_label_1, const ged::GXLLabel &edge_label_2) const;

    vector<double> ahneman_node_label(const ged::GXLLabel &node_label) const;

    vector<double> ahneman_edge_label(const ged::GXLLabel &edge_label) const;

    vector<double> catalyst_node_label(const ged::GXLLabel &node_label) const;

    vector<double> catalyst_edge_label(const ged::GXLLabel &edge_label) const;

    virtual double
    labelkernel(const vector<double> &label_1, const vector<double> &label_2, const double &parameter) const = 0;
};

Molecule::Molecule(const double &base_node_ins, const double &base_node_del, const double &base_node_rel,
                   const double &base_edge_ins, const double &base_edge_del, const double &base_edge_rel,
                   const MolecularDataset &dataset) {
  this->base_node_ins = base_node_ins;
  this->base_node_del = base_node_del;
  this->base_node_rel = base_node_rel;
  this->base_edge_ins = base_edge_ins;
  this->base_edge_del = base_edge_del;
  this->base_edge_rel = base_edge_rel;
  switch (dataset) {
    case MolecularDataset::AHNEMAN:
      this->node_label_extractor = [this](const GXLLabel &node_label) {
          return this->ahneman_node_label(node_label);
      };
      this->edge_label_extractor = [this](const GXLLabel &edge_label) {
          return this->ahneman_edge_label(edge_label);
      };
      break;
    case MolecularDataset::CATALYST:
      this->node_label_extractor = [this](const GXLLabel &node_label) {
          return this->catalyst_node_label(node_label);
      };
      this->edge_label_extractor = [this](const GXLLabel &edge_label) {
          return this->catalyst_edge_label(edge_label);
      };
      break;
  }
}

double Molecule::node_ins_cost_fun(const GXLLabel &node_label) const {
  return this->base_node_ins;
}

double Molecule::node_del_cost_fun(const GXLLabel &node_label) const {
  return this->base_node_del;
}

double Molecule::node_rel_cost_fun(const GXLLabel &node_label_1, const GXLLabel &node_label_2) const {
  if (node_label_1.at("AtomicNum") == node_label_2.at("AtomicNum")) {
    return 0;
  }
  return this->labelkernel(this->node_label_extractor(node_label_1), this->node_label_extractor(node_label_2),
                           this->base_node_rel);
}

double Molecule::edge_ins_cost_fun(const GXLLabel &edge_label) const {
  return this->base_edge_ins;
}

double Molecule::edge_del_cost_fun(const GXLLabel &edge_label) const {
  return this->base_edge_del;
}

double Molecule::edge_rel_cost_fun(const GXLLabel &edge_label_1, const GXLLabel &edge_label_2) const {
  if (edge_label_1.at("BondType") == edge_label_2.at("BondType")) {
    return 0.0;
  }
  return this->labelkernel(this->edge_label_extractor(edge_label_1), this->edge_label_extractor(edge_label_2),
                           this->base_edge_rel);
}

vector<double> Molecule::ahneman_node_label(const GXLLabel &node_label) const {
  return vector<double>();
}

vector<double> Molecule::ahneman_edge_label(const GXLLabel &edge_label) const {
  return vector<double>();
}

vector<double> Molecule::catalyst_node_label(const GXLLabel &node_label) const {
  vector<double> labels(4);
  labels[0] = stod(node_label.at("x"));
  labels[1] = stod(node_label.at("y"));
  labels[2] = stod(node_label.at("z"));
  labels[3] = stod(node_label.at("charge"));
  return labels;
}

vector<double> Molecule::catalyst_edge_label(const GXLLabel &edge_label) const {
  return vector<double>();
}


class Euclidean : public Molecule {
public:
    Euclidean(const double &base_node_ins, const double &base_node_del, const double &base_node_rel,
              const double &base_edge_ins, const double &base_edge_del, const double &base_edge_rel,
              const MolecularDataset &dataset) :
        Molecule(base_node_ins, base_node_del, base_node_rel, base_edge_ins, base_edge_del, base_edge_rel, dataset) {};

    double labelkernel(const vector<double> &label_1, const vector<double> &label_2, const double &parameter) const;
};

double
Euclidean::labelkernel(const vector<double> &label_1, const vector<double> &label_2, const double &parameter) const {
  if (label_1.empty() && label_2.empty()) {
    return parameter;
  }
  double value = 0.0;
  for (uint i = 0; i < label_1.size(); i++) {
    value += pow(label_1[i] - label_2[i], 2);
  }
  return parameter * sqrt(value);
}

class Gaussian : public Molecule {
public:
    Gaussian(const double &base_node_ins, const double &base_node_del, const double &base_node_rel,
             const double &base_edge_ins, const double &base_edge_del, const double &base_edge_rel,
             const MolecularDataset &dataset) :
        Molecule(base_node_ins, base_node_del, base_node_rel, base_edge_ins, base_edge_del, base_edge_rel, dataset) {};

    double labelkernel(const vector<double> &label_1, const vector<double> &label_2, const double &parameter) const;
};

double
Gaussian::labelkernel(const vector<double> &label_1, const vector<double> &label_2, const double &parameter) const {
  if (label_1.empty() && label_2.empty()) {
    return parameter;
  }
  double value = 0.0;
  for (uint i = 0; i < label_1.size(); i++) {
    value += pow(label_1[i] - label_2[i], 2);
  }
  return 1 - exp(-1 * (value / 2 * pow(parameter, 2)));
}


#endif //STRUKERN_MOLECULE_H
