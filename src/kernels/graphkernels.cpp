#include "graphkernels.h"
//#include "src/edit_costs/chem_1.hpp"
//#include "src/edit_costs/edit_costs.hpp"
#include <Dense>
//#include "../libs/hungarian-algorithm-cpp/Hungarian.h"
#include "Hungarian.cpp"

using namespace std;


template<typename node_label_type, typename edge_label_type>
double OptimalMolecularAssignmentKernel<node_label_type, edge_label_type>::opt_assign_mol(
    const UndirectedGraph<node_label_type, edge_label_type> &x1,
    const UndirectedGraph<node_label_type, edge_label_type> &x2, const KernelParams &params) const{
  uint N = x1.number_of_nodes();
  uint M = x2.number_of_nodes();
  vector<double> costmatrixrow;
  vector<vector<double>> costmatrix;
  double val, max_element = 0.0;
  for (uint i = 0; i < N; i++) {
    for (uint j = 0; j < M; j++) {
      val = this->neighbour_kernel(x1, x2, i, j, params);
      costmatrixrow.push_back(val);
      if (val > max_element) {
        max_element = val;
      }
    }
    costmatrix.push_back(costmatrixrow);
    costmatrixrow.clear();
  }
  for (uint i = 0; i < N; i++) {
    for (uint j = 0; j < M; j++) {
      costmatrix[i][j] = max_element - costmatrix[i][j];
    }
  }
  vector<int> assignment;
  HungarianAlgorithm solver;
  double cost = solver.Solve(costmatrix, assignment);
  return N * max_element - cost;
}


template<typename node_label_type, typename edge_label_type>
double OptimalMolecularAssignmentKernel<node_label_type, edge_label_type>::dot(
    const UndirectedGraph<node_label_type, edge_label_type> &x1,
    const UndirectedGraph<node_label_type, edge_label_type> &x2,
    const KernelParams &params) const{
  if (x1.number_of_nodes() <= x2.number_of_nodes()) {
    return this->opt_assign_mol(x1, x2, params);
  } else {
    return this->opt_assign_mol(x2, x1, params);
  }
}

template<typename node_label_type, typename edge_label_type>
double OptimalMolecularAssignmentKernel<node_label_type, edge_label_type>::neighbour_kernel(
    const UndirectedGraph<node_label_type, edge_label_type> &x1,
    const UndirectedGraph<node_label_type, edge_label_type> &x2,
    const uint &node_1_id, const uint &node_2_id, const KernelParams &params) const{
  double Rl = 0.0;
  for (uint i = 1; i < params.MolecularRadius + 1; i++) {
    Rl += this->Rl(x1, x2, node_1_id, node_2_id, params, i) * (pow(params.AssignmentDecay, i));
  }
  return this->atom_kernel(x1, x2, node_1_id, node_2_id, params) + this->R0(x1, x2, node_1_id, node_2_id, params) + Rl;
}

template<typename node_label_type, typename edge_label_type>
double OptimalMolecularAssignmentKernel<node_label_type, edge_label_type>::atom_kernel(
    const UndirectedGraph<node_label_type, edge_label_type> &x1,
    const UndirectedGraph<node_label_type, edge_label_type> &x2,
    const uint &node_1_id, const uint &node_2_id, const KernelParams &params) const{
  return atom_kernel_function<node_label_type>(x1.get_node(node_1_id).attribute, x2.get_node(node_2_id).attribute, params);
}

template<typename node_label_type, typename edge_label_type>
double OptimalMolecularAssignmentKernel<node_label_type, edge_label_type>::bond_kernel(
    const UndirectedGraph<node_label_type, edge_label_type> &x1,
    const UndirectedGraph<node_label_type, edge_label_type> &x2,
    const uint &node_1_id, const uint &node_2_id, const uint &neighbour_id_1, const uint &neighbour_id_2,
    const KernelParams &params) const{
  return bond_kernel_function<edge_label_type>(x1.get_edge(node_1_id, neighbour_id_1).attribute,
                               x2.get_edge(node_2_id, neighbour_id_2).attribute, params);
}

template<typename node_label_type, typename edge_label_type>
double OptimalMolecularAssignmentKernel<node_label_type, edge_label_type>::opt_assign_R0(
    const UndirectedGraph<node_label_type, edge_label_type> &x1,
    const UndirectedGraph<node_label_type, edge_label_type> &x2,
    const uint &node_1_id, const uint &node_2_id, const KernelParams &params) const{
  uint N = x1.num_neighbours(node_1_id);
  uint M = x2.num_neighbours(node_2_id);
  double val, max_element = 0.0;
  vector<double> costmatrixrow;
  vector<vector<double>> costmatrix;
  for (const auto &neighbour1 : x1.get_neighbours(node_1_id)) {

    for (const auto &neighbour2 : x2.get_neighbours(node_2_id)) {
      val = this->atom_kernel(x1, x2, neighbour1, neighbour2, params) *
            this->bond_kernel(x1, x2, node_1_id, node_2_id, neighbour1, neighbour2, params);
      costmatrixrow.push_back(val);
      if (val > max_element) {
        max_element = val;
      }
    }
    costmatrix.push_back(costmatrixrow);
    costmatrixrow.clear();
  }
  for (uint i = 0; i < N; i++) {
    for (uint j = 0; j < M; j++) {
      costmatrix[i][j] = max_element - costmatrix[i][j];

    }
  }
  vector<int> assignment;
  HungarianAlgorithm solver;
  double cost = solver.Solve(costmatrix, assignment);
  return (N * max_element - cost) / M;
}

template<typename node_label_type, typename edge_label_type>
double OptimalMolecularAssignmentKernel<node_label_type, edge_label_type>::R0(
    const UndirectedGraph<node_label_type, edge_label_type> &x1,
    const UndirectedGraph<node_label_type, edge_label_type> &x2,
    const uint &node_1_id, const uint &node_2_id, const KernelParams &params) const{
  unordered_set<uint> node_1_neighbours, node_2_neighbours;
  if (x1.num_neighbours(node_1_id) < x2.num_neighbours(node_2_id)) {
    return this->opt_assign_R0(x1, x2, node_1_id, node_2_id, params);
  } else {
    return this->opt_assign_R0(x2, x1, node_2_id, node_1_id, params);
  }
}

template<typename node_label_type, typename edge_label_type>
double OptimalMolecularAssignmentKernel<node_label_type, edge_label_type>::Rl(
    const UndirectedGraph<node_label_type, edge_label_type> &x1,
    const UndirectedGraph<node_label_type, edge_label_type> &x2,
    const uint &node_1_id, const uint &node_2_id, const KernelParams &params, uint &it_num) const{
  if (it_num == 1) {
    double R1 = 0.0, count = 0;
    for (const auto &neighbour1 : x1.get_neighbours(node_1_id)) {
      for (const auto &neighbour2 : x2.get_neighbours(node_2_id)) {
        R1 += this->R0(x1, x2, neighbour1, neighbour2, params);
        count++;
      }
    }
    return R1 / count;
  } else {
    uint new_it = it_num - 1;
    double R_minus1 = 0.0, count = 0;
    for (const auto &neighbour1 : x1.get_neighbours(node_1_id)) {
      for (const auto &neighbour2 : x2.get_neighbours(node_2_id)) {
        R_minus1 += this->Rl(x1, x2, neighbour1, neighbour2, params, new_it);
        count++;
      }
    }
    return R_minus1 / count;
  }
}