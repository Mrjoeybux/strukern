#ifndef GRAPHKERNELS
#define GRAPHKERNELS
#include "../domain/edit_costs.hpp"
#include "basekernel.h"
#include "src/edit_costs/edit_costs.hpp"
#include "src/env/ged_env.hpp"
#include <vector>
using namespace std;
using namespace ged;

template <typename bandwidthType, typename gedNodeLabels, typename gedEdgeLabels> class gedKernel : public Kernel<string, bandwidthType> {

public:
  MatrixXd computeKernelMatrix(const vector<string> &data, const unordered_map<string, bandwidthType> &params) const;

  /* use computeKernelMatrix when applicable */
  double dot(const string &x1, const string &x2, const unordered_map<string, bandwidthType> &params) const;
  /*Adds edit costs and initiliases environment*/
  virtual EditCosts<gedNodeLabels, gedEdgeLabels> *init_edit_costs(GEDEnv<int, gedNodeLabels, gedEdgeLabels> &env,
                                                                   const unordered_map<string, bandwidthType> &params) const = 0;
  /*Adds methods and initialises method. */
  virtual void init_methods(GEDEnv<int, gedNodeLabels, gedEdgeLabels> &env, const unordered_map<string, bandwidthType> &params) const = 0;
};

class gedMABTS : public gedKernel<LabelPairMap, string, string> {
public:
  EditCosts<string, string> *init_edit_costs(GEDEnv<int, string, string> &env, const unordered_map<string, LabelPairMap> &params) const;

  virtual void init_methods(GEDEnv<int, string, string> &env, const unordered_map<string, LabelPairMap> &params) const = 0;
};

class gedMABTS_LSAPE : public gedMABTS {

  void init_methods(GEDEnv<int, string, string> &env, const unordered_map<string, LabelPairMap> &params) const;
};
#endif /* GRAPHKERNELS */
