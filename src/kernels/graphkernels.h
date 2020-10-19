#ifndef GRAPHKERNELS
#define GRAPHKERNELS
#include "basekernel.h"
#include <functional>
#include <vector>
#include <dlib/optimization/max_cost_assignment.h>
using namespace std;

/*class WeisfeilerLehman: public Kernel<string>{

    double dot(const string &x1, const string &x2, const KernelParams &params) const;

    double computeKernelMatrix(const vector<string> &data, const KernelParams &params) const;
};*/
template<typename T, typename base_T>
class OptimalAssignmentKernel : Kernel<T>{
private:
    const Kernel<base_T> *basekernel;

public:

    OptimalAssignmentKernel(Kernel<base_T> *basekernel){this->basekernel = basekernel;}

    matrix virtual build_cost_matrix(const T &x1, const T &x2, const KernelParams &params);


};


#endif /* GRAPHKERNELS */
