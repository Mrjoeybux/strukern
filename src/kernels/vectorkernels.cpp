//
// Created by Mrjoeybux on 12/05/2020.
//
#include "./vectorkernels.h"
#include <math.h>

MatrixXd VectorKernel::computeKernelMatrix(const MatrixXd &data, const KernelParams &params) const {
// data shape = (n by m)
// n = number of examples, m = length of each example
  uint n = data.rows();
  MatrixXd K = MatrixXd::Zero(n, n);
  for (uint i = 0; i < n; i++) {
    for (uint j = i; j < n; j++) {
      K(i, j) = this->dot(data.row(i), data.row(j), params);
    }
  }
  return K;
}

bool isNanVector(const VectorXd &x) { return isnan(x[0]); };

VectorKernel::VectorKernel(const double &default_val) : Kernel<VectorXd>(default_val) {};

double RBFKernel::dot(const VectorXd &x1, const VectorXd &x2, const KernelParams &params) const {
  VectorXd diff = x1 - x2;
  return std::exp(-1.0 * (diff.squaredNorm() / (2.0 * pow(params.RBFSigma, 2))));
}

MatrixXd
RBFKernel::computeRectangularKernelMatrix(const vector<VectorXd> &data1, const vector<VectorXd> &data2,
                                          const KernelParams &params) const {
vector<VectorXd> compact_data1(data1), compact_data2(data2);
  compact_data1.erase(remove_if(compact_data1.begin(), compact_data1.end(), isNanVector), compact_data1.end());
  compact_data2.erase(remove_if(compact_data2.begin(), compact_data2.end(), isNanVector), compact_data2.end());
  int n = compact_data1.size(), m = compact_data2.size();
  MatrixXd K = MatrixXd::Zero(n, m);
  for(uint i = 0; i < n; i++){
    for(uint j = 0; j < m; j++){
      K(i, j) = this->dot(compact_data1[i], compact_data2[j], params);
    }
  }
  return K;

}

double Polynomial::dot(const VectorXd &x1, const VectorXd &x2, const KernelParams &params) const {
  return pow(x1.dot(x2), params.PolynomialPower);
}