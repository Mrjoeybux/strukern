//
// Created by Mrjoeybux on 12/05/2020.
//
#ifndef STRUKERN_VECTORKERNELS_H
#define STRUKERN_VECTORKERNELS_H
#include <Dense>
#include "basekernel.h"
using namespace Eigen;

class VectorKernel : public Kernel<VectorXd>{
public:

    virtual double dot(const VectorXd &x1, const VectorXd &x2, const KernelParams &params) const = 0;

    MatrixXd computeKernelMatrix(const MatrixXd &data, const KernelParams &params) const;
};

class RBF : public VectorKernel{
public:

    double dot(const VectorXd &x1, const VectorXd &x2, const KernelParams &params) const;
};

class Polynomial : public VectorKernel{
public:

    double dot(const VectorXd &x1, const VectorXd &x2, const KernelParams &params) const;
};

#endif //STRUKERN_VECTORKERNELS_H
