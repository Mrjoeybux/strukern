#ifndef GRAPHKERNELS
#define GRAPHKERNELS
#include "../domain/edit_costs.h"
#include "basekernel.h"
#include "src/edit_costs/edit_costs.hpp"
#include "src/env/ged_env.hpp"
#include <functional>
#include <vector>
using namespace std;
using namespace ged;

enum class GEDEditCosts { MABTS };

enum class GEDMethods { BIPARTITE };

class GEDKernel : public Kernel<string> {
private:
  // function<EditCosts<string, string> *(GEDEnv<int, string, string> &, const KernelParams &)> edit_cost_init;

  // function<void(GEDEnv<int, string, string> &, const KernelParams &)> method_init;

public:
  GEDKernel(const GEDEditCosts &edit_costs = GEDEditCosts::MABTS, const GEDMethods &method = GEDMethods::BIPARTITE);

  // MatrixXd computeKernelMatrix(const vector<string> &data, const KernelParams &params) const;

  // use computeKernelMatrix when applicable
  double dot(const string &x1, const string &x2, const KernelParams &params) const;
  // Adds edit costs and initiliases environment
  // EditCosts<string, string> *init_edit_costs(GEDEnv<int, string, string> &env, const KernelParams &params) const;
  // Adds methods and initialises method.
  // void init_methods(GEDEnv<int, string, string> &env, const KernelParams &params) const;

  // EditCosts<string, string> *init_EditCosts_MABTS(GEDEnv<int, string, string> &env, const KernelParams &params) const;

  // EditCosts<string, string> *init_EditCosts_CONSTANT(GEDEnv<int, string, string> &env, const KernelParams &params) const;

  // void init_Method_BIPARTITE(GEDEnv<int, string, string> &env, const KernelParams &params) const;
};

#endif /* GRAPHKERNELS */
