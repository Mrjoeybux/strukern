#include "edit_costs.hpp"
double GaussianDistance::node_ins_cost_fun(const GXLLabel &node_label) const { return 1.0; }

double GaussianDistance::node_del_cost_fun(const GXLLabel &node_label) const { return 1.0; }

double GaussianDistance::node_rel_cost_fun(const GXLLabel &node_label_1, const GXLLabel &node_label_2) const { return 0.0; }

double GaussianDistance::edge_ins_cost_fun(const GXLLabel &edge_label) const { return 1.0; }

double GaussianDistance::edge_del_cost_fun(const GXLLabel &edge_label) const { return 1.0; }

double GaussianDistance::edge_rel_cost_fun(const GXLLabel &edge_label_1, const GXLLabel &edge_label_2) const { return 0.0; }