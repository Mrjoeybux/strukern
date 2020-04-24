#include "edit_costs.hpp"

double EuclideanDistance::node_ins_cost_fun(const GXLLabel &node_label) const { return 1.0; }

double EuclideanDistance::node_del_cost_fun(const GXLLabel &node_label) const { return 1.0; }

double EuclideanDistance::node_rel_cost_fun(const GXLLabel &node_label_1, const GXLLabel &node_label_2) const { return 0.0; }

double EuclideanDistance::edge_ins_cost_fun(const GXLLabel &edge_label) const { return 1.0; }

double EuclideanDistance::edge_del_cost_fun(const GXLLabel &edge_label) const { return 1.0; }

double EuclideanDistance::edge_rel_cost_fun(const GXLLabel &edge_label_1, const GXLLabel &edge_label_2) const { return 0.0; }