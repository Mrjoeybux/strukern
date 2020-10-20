#ifndef GRAPHKERNELS
#define GRAPHKERNELS

#include "basekernel.h"
#include <functional>
#include <vector>
#include "../domain/graph.h"
#include "../kernels/vectorkernels.h"

using namespace std;

template<typename node_label_type>
double atom_kernel_function(const node_label_type &x1, const node_label_type &x2, const KernelParams &params){
  // If atom element types are equal
  if(x1[0] == x2[0]){
    double val = 0.0;
    double double_squared_sigma = 2*pow(params.RBFSigma, 2);
    for(uint i = 1; i < x1.size(); i++){
      val += pow(x1[i] - x2[i], 2);
    }
    return exp((-1.0 * val) / double_squared_sigma);
  }
  else{
    return 0.0;
  }
}

template<typename edge_label_type>
double bond_kernel_function(const edge_label_type &x1, const edge_label_type &x2, const KernelParams &params){
  // If bond orders are equal
  if(x1[0] == x2[0]) {
    // TODO: CHECK WHAT CAN BE COMPUTED
    // If both in aromatic system
    if (x1[1] == x2[1]) {
      double val = 0.0;
      double double_squared_sigma = 2 * pow(params.RBFSigma, 2);
      for (uint i = 2; i < x1.size(); i++) {
        val += pow(x1[i] - x2[i], 2);
      }
      return exp((-1.0 * val) / double_squared_sigma);
    }
    else{
      return 0.0;
    }
  }
  else{
    return 0.0;
  }
}




template<typename node_label_type, typename edge_label_type>
class OptimalMolecularAssignmentKernel : public Kernel<UndirectedGraph<node_label_type, edge_label_type>>{
public:

    double dot(const UndirectedGraph<node_label_type, edge_label_type> &x1, const UndirectedGraph<node_label_type, edge_label_type> &x2,
               const KernelParams &params) const;

    double opt_assign_mol(const UndirectedGraph<node_label_type, edge_label_type> &x1, const UndirectedGraph<node_label_type, edge_label_type> &x2,
                          const KernelParams &params) const;

    double opt_assign_R0(const UndirectedGraph<node_label_type, edge_label_type> &x1, const UndirectedGraph<node_label_type, edge_label_type> &x2,
        const uint &node_1_id, const uint &node_2_id, const KernelParams &params) const;

    double
    neighbour_kernel(const UndirectedGraph<node_label_type, edge_label_type> &x1, const UndirectedGraph<node_label_type, edge_label_type> &x2,
             const uint &node_1_id, const uint &node_2_id, const KernelParams &params) const;

    double
    atom_kernel(const UndirectedGraph<node_label_type, edge_label_type> &x1, const UndirectedGraph<node_label_type, edge_label_type> &x2,
                const uint &node_1_id, const uint &node_2_id, const KernelParams &params) const;

    double
    bond_kernel(const UndirectedGraph<node_label_type, edge_label_type> &x1, const UndirectedGraph<node_label_type, edge_label_type> &x2,
                const uint &node_1_id, const uint &node_2_id, const uint &neighbour_id_1, const uint &neighbour_id_2,
                const KernelParams &params) const;

    double R0(const UndirectedGraph<node_label_type, edge_label_type> &x1, const UndirectedGraph<node_label_type, edge_label_type> &x2,
              const uint &node_1_id, const uint &node_2_id, const KernelParams &params) const;

    double Rl(const UndirectedGraph<node_label_type, edge_label_type> &x1, const UndirectedGraph<node_label_type, edge_label_type> &x2,
              const uint &node_1_id, const uint &node_2_id, const KernelParams &params, uint &it_num) const;
};


#endif /* GRAPHKERNELS */
