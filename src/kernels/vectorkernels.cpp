//
// Created by Mrjoeybux on 12/05/2020.
//
#include "./vectorkernels.h"

MatrixXd VectorKernel::computeKernelMatrix(const MatrixXd &data, const KernelParams &params) const {
// data shape = (n by m)
// n = number of examples, m = length of each example
uint n = data.rows(), m = data.cols();
MatrixXd K = MatrixXd::Zero(n, n);
for(uint i = 0; i < n; i++){
  for(uint j = i; j < n; j++){
    K(i, j) = this->dot(data.row(i), data.row(j), params);
  }
}
return K;
}

double RBF::dot(const VectorXd &x1, const VectorXd &x2, const KernelParams &params) const {
  VectorXd diff = x1 - x2;
  return exp(-0.5*diff.squaredNorm() / params.RBFSigma);
}

double Polynomial::dot(const VectorXd &x1, const VectorXd &x2, const KernelParams &params) const {
  return pow(x1.dot(x2), params.PolynomialPower);
}