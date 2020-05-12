#include "graphkernels.h"
//#include "src/edit_costs/chem_1.hpp"
//#include "src/edit_costs/edit_costs.hpp"
#include <Dense>
using namespace std;


double WeisfeilerLehman::dot(const string &x1, const string &x2, const KernelParams &params) const {
  return 0;
}

double WeisfeilerLehman::computeKernelMatrix(const vector<string> &data, const KernelParams &params) const {
  return Kernel::computeKernelMatrix(data, params);
}
