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
    VectorKernel(const double &default_val = numeric_limits<double>::quiet_NaN());

    virtual double dot(const VectorXd &x1, const VectorXd &x2, const KernelParams &params) const = 0;

    MatrixXd computeKernelMatrix(const MatrixXd &data, const KernelParams &params) const;

    //MatrixXd computeRectangularKernelMatrix(const MatrixXd &data1, const MatrixXd &data2, const KernelParams &params) const ;
};

class RBFKernel : public VectorKernel{
public:
    RBFKernel(): VectorKernel(1.0){};

    double dot(const VectorXd &x1, const VectorXd &x2, const KernelParams &params) const;

    MatrixXd
    computeRectangularKernelMatrix(const vector<VectorXd> &data1, const vector<VectorXd> &data2, const KernelParams &params) const;
};

class Polynomial : public VectorKernel{
public:

    double dot(const VectorXd &x1, const VectorXd &x2, const KernelParams &params) const;
};

#endif //STRUKERN_VECTORKERNELS_H
