#ifndef GRAPHKERNELS
#define GRAPHKERNELS
#include "basekernel.h"
#include <functional>
#include <vector>
using namespace std;

class WeisfeilerLehman: public Kernel<string>{

    double dot(const string &x1, const string &x2, const KernelParams &params) const;

    double computeKernelMatrix(const vector<string> &data, const KernelParams &params) const;
};

#endif /* GRAPHKERNELS */
