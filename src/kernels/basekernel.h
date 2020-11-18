#ifndef BASEKERNEL
#define BASEKERNEL

#include "../domain/common.h"
#include <iostream>
#include <Core>
//#include <eigen3/Eigen/Core>
#include <memory>
#include <Dense>
#include <Eigenvalues>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>
#include <functional>

#ifndef LSAPE_IndexType
#define LSAPE_IndexType std::size_t
#endif
#include "lsap.h"
//#include "../lib/gedlib/ext/lsape.5/include/lsap.h"

using namespace std;
using namespace Eigen;

struct KernelParams {
    int StringCompressionLevel;
    int JPEGCompressionQuality;
    int TensorCompressionLevel;
    unordered_map<string, LabelPairMap> MABTS;
    unordered_map<string, int> LocalityImproved;
    unordered_map<string, double> MolecularEditCosts;
    double RBFSigma;
    int PolynomialPower;
    int TensorConcatDim;
    int MolecularRadius;
    double AssignmentDecay;
};

#define NUMERICAL_PRECISION 1e-12
#define VERYLARGECOST 1e12

enum class KernelType {
    Hilbert, Krein
};


enum class MultiInstanceMethod {
    SUM, OA, MAX, MIN
};


enum class CompressionDistanceMeasure {
    NCD, DIFF
};

template<typename T>
class Kernel {

public:
    Kernel(const double &default_val = numeric_limits<double>::quiet_NaN()){this->default_val = default_val;};

    double default_val;

    virtual double dot(const T &x1, const T &x2, const KernelParams &params) const = 0;

    virtual MatrixXd computeKernelMatrix(const vector<T> &data, const KernelParams &params) const;

    virtual MatrixXd
    computeRectangularKernelMatrix(const vector<T> &data1, const vector<T> &data2, const KernelParams &params) const;

    static MatrixXd normaliseHilbert(const MatrixXd &K);

    static MatrixXd normaliseKrein(const MatrixXd &K);

    /*virtual VectorXd quadratic_form_kmat_gradient(const VectorXd &u, const MatrixXd &kmat, const VectorXd &v, const vector<T> &X,
                                                  const KernelParams &params) const;

    virtual VectorXd quadratic_form_kmat_block_gradient(const VectorXd &u, const MatrixXd &kmat_block, const VectorXd &v, const vector<T> &X,
                                                        const vector<int> &support_args, const vector<int> &block_args, const MatrixXd &kmat,
                                                        const KernelParams &params) const;
    KernelType ktype;*/
};

template<typename T>
inline MatrixXd Kernel<T>::computeKernelMatrix(const vector<T> &data, const KernelParams &params) const {
  int n = data.size();
  MatrixXd K = MatrixXd::Zero(n, n);

  for (uint i = 0; i < n; i++) {
    for (uint j = i; j < n; j++) {
      K(i, j) = this->dot(data[i], data[j], params);
    }
  }
  return K.selfadjointView<Upper>();
}

template<typename T>
inline MatrixXd Kernel<T>::computeRectangularKernelMatrix(const vector<T> &data1, const vector<T> &data2,
                                                          const KernelParams &params) const {
  int n = data1.size(), m = data2.size();
  MatrixXd K = MatrixXd::Zero(n, m);

  for (uint i = 0; i < n; i++) {
    for (uint j = 0; j < m; j++) {
      K(i, j) = this->dot(data1[i], data2[j], params);
    }
  }
  return K;
}

template<typename T>
inline MatrixXd Kernel<T>::normaliseHilbert(const MatrixXd &K) {
  VectorXd kvec = K.diagonal().cwiseInverse().cwiseSqrt();
  return K.cwiseProduct(kvec * kvec.transpose());
}

template<typename T>
inline MatrixXd Kernel<T>::normaliseKrein(const MatrixXd &K) {
  SelfAdjointEigenSolver<MatrixXd> solver(K);
  MatrixXd Sigmap = solver.eigenvalues().cwiseMax(0).asDiagonal();
  MatrixXd Sigmam = solver.eigenvalues().cwiseMin(0).cwiseAbs().asDiagonal();
  return Kernel<T>::normaliseHilbert(solver.eigenvectors() * Sigmap * solver.eigenvectors().transpose()) -
                   Kernel<T>::normaliseHilbert(solver.eigenvectors() * Sigmam * solver.eigenvectors().transpose());
}

template<typename T>
class MultiInstanceKernel : public Kernel<vector<T>> {
private:
    const Kernel<T> *basekernel;

    function<double(const vector<T> &, const vector<T> &, const KernelParams &)> dot_func;

public:
    MultiInstanceKernel(const Kernel<T> *basekernel, const MultiInstanceMethod &method);

    double dot(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const;

    double sum(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const;

    double optimalassignment(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const;

    double max(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const;

    double min(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const;

    MatrixXd computeKernelMatrix(const vector<vector<T>> &data, const KernelParams &params) const;
};

template<typename T>
MultiInstanceKernel<T>::MultiInstanceKernel(const Kernel<T> *basekernel, const MultiInstanceMethod &method) {
  this->basekernel = basekernel;
  switch (method) {
    case MultiInstanceMethod::SUM:
      this->dot_func = [this](const vector<T> &x1, const vector<T> &x2, const KernelParams &params) {
          return this->sum(x1, x2, params);
      };
      break;
    case MultiInstanceMethod::OA:
      this->dot_func = [this](const vector<T> &x1, const vector<T> &x2, const KernelParams &params) {
          return this->optimalassignment(x1, x2, params);
      };
      break;
    case MultiInstanceMethod::MAX:
      this->dot_func = [this](const vector<T> &x1, const vector<T> &x2, const KernelParams &params) {
          return this->max(x1, x2, params);
      };
      break;
    case MultiInstanceMethod::MIN:
      this->dot_func = [this](const vector<T> &x1, const vector<T> &x2, const KernelParams &params) {
          return this->min(x1, x2, params);
      };
      break;
    default:
      cout << "Unknown Method!" << endl;
      break;
  }
}

template<typename T>
double MultiInstanceKernel<T>::dot(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const {
  return this->dot_func(x1, x2, params);
}

template<typename T>
double MultiInstanceKernel<T>::sum(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const {
  if (!isnan(this->basekernel->default_val)) {
    if (x1 == x2) {
      return this->basekernel->default_val;
    }
  }
  MatrixXd K = this->basekernel->computeRectangularKernelMatrix(x1, x2, params);
  return K.sum();
}

template<typename T>
double MultiInstanceKernel<T>::max(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const {
  if(!isnan(this->basekernel->default_val)){
    if(x1 == x2){
      return this->basekernel->default_val;
    }
  }
  MatrixXd K = this->basekernel->computeRectangularKernelMatrix(x1, x2, params);
  return K.maxCoeff();
}

template<typename T>
double MultiInstanceKernel<T>::min(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const {
  if(!isnan(this->basekernel->default_val)){
    if(x1 == x2){
      return this->basekernel->default_val;
    }
  }
  MatrixXd K = this->basekernel->computeRectangularKernelMatrix(x1, x2, params);
  return K.minCoeff();
}

template<typename T>
double
MultiInstanceKernel<T>::optimalassignment(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const {
  size_t nRows = x1.size(), nCols = x2.size();
  double *C = new double[nRows * nCols];
  size_t *rho = new size_t[nRows];
  double *u = new double[nRows], *v = new double[nCols];
  for (uint col = 0; col < nCols; col++) {
    for (uint row = 0; row < nRows; row++) {
      C[nRows * col + row] = this->basekernel->dot(x1[row], x2[col], params);
    }
  }
  lsape::hungarianLSAP<double>(C, nRows, nCols, rho, u, v);
  delete[] C, rho;
  double cost = lsape::permCost(u, nRows, v, nCols);
  delete[] u, v;
  return cost;
}

template<typename T>
MatrixXd MultiInstanceKernel<T>::computeKernelMatrix(const vector<vector<T>> &data, const KernelParams &params) const {
  int n = data.size();
  MatrixXd K = MatrixXd::Zero(n, n);

  for (uint i = 0; i < n; i++) {
    for (uint j = i; j < n; j++) {
      K(i, j) = this->dot(data[i], data[j], params);
      if (i != j) {
        K(j, i) = K(i, j);
      }
    }
  }
  return K;
}

template<typename T>
class AbstractCompressionKernel : public Kernel<T> {

public:
     AbstractCompressionKernel(const CompressionDistanceMeasure &measure){
       switch (measure) {
         case CompressionDistanceMeasure::NCD:
           this->distance_measure = [this](const double &comp_x1, const double &comp_x2, const double &comp_x1_x2, const double &comp_x2_x1){
               return this->NCD(comp_x1, comp_x2, comp_x1_x2, comp_x2_x1);
           };
           break;
         case CompressionDistanceMeasure::DIFF:
           this->distance_measure = [this](const double &comp_x1, const double &comp_x2, const double &comp_x1_x2, const double &comp_x2_x1){
               return this->difference(comp_x1, comp_x2, comp_x1_x2, comp_x2_x1);
           };
           break;
         default:
           cout << "Unknown distance measure!" << endl;
       }
     }

    function<double(const double &comp_x1, const double &comp_x2, const double &comp_x1_x2, const double &comp_x2_x1)> distance_measure;

    double dot(const T &x1, const T &x2, const KernelParams &params) const;

    MatrixXd computeKernelMatrix(const vector<T> &data, const KernelParams &params) const;

    MatrixXd computeRectangularKernelMatrix(const vector<T> &data1, const vector<T> &data2,
                                            const KernelParams &params) const;

    double NCD(const double &comp_x1, const double &comp_x2, const double &comp_x1_x2, const double &comp_x2_x1) const;

    double difference(const double &comp_x1, const double &comp_x2, const double &comp_x1_x2, const double &comp_x2_x1) const;

    virtual double compress(const T &x, const KernelParams &params) const = 0;

    virtual T concat(const T &x1, const T & x2, const KernelParams &params) const = 0;
};

template<typename T>
double AbstractCompressionKernel<T>::dot(const T &x1, const T &x2, const KernelParams &params) const {
  double comp_x1, comp_x2, comp_x1_x2, comp_x2_x1;
  comp_x1 = this->compress(x1, params);
  comp_x2 = this->compress(x2, params);
  comp_x1_x2 = this->compress(this->concat(x1, x2, params), params);
  comp_x2_x1 = this->compress(this->concat(x2, x1, params), params);
  return this->distance_measure(comp_x1, comp_x2, comp_x1_x2, comp_x2_x1);
}

template<typename T>
MatrixXd AbstractCompressionKernel<T>::computeKernelMatrix(const vector<T> &data, const KernelParams &params) const {
  int n = data.size();
  MatrixXd K = MatrixXd::Zero(n, n);
  double comp_i_j, comp_j_i;
  vector<double> comp_single;
  for(uint i = 0; i < n; i++){
    comp_single.push_back(this->compress(data[i], params));
  }
  for (uint i = 0; i < n; i++) {
    for (uint j = i; j < n; j++) {
      if(i == j){
        comp_i_j = this->compress(this->concat(data[i], data[j], params), params);
        K(i, i) = this->distance_measure(comp_single[i], comp_single[i], comp_i_j, comp_i_j);
      }else {
        comp_i_j = this->compress(this->concat(data[i], data[j], params), params);
        comp_j_i = this->compress(this->concat(data[j], data[i], params), params);
        K(i, j) = this->distance_measure(comp_single[i], comp_single[j], comp_i_j, comp_j_i);
      }
    }
  }
  return K.selfadjointView<Upper>();
}

template<typename T>
MatrixXd AbstractCompressionKernel<T>::computeRectangularKernelMatrix(const vector<T> &data1, const vector<T> &data2,
                                                                      const KernelParams &params) const {
  int n1 = data1.size();
  int n2 = data2.size();
  MatrixXd K = MatrixXd::Zero(n1, n2);
  double comp_i_j, comp_j_i;
  vector<double> comp_single_1, comp_single_2;
  for(uint i = 0; i < n1; i++){
    comp_single_1.push_back(this->compress(data1[i], params));
  }
  for(uint i = 0; i < n2; i++){
    comp_single_2.push_back(this->compress(data2[i], params));
  }
  for (uint i = 0; i < n1; i++) {
    for (uint j = 0; j < n2; j++) {
        comp_i_j = this->compress(this->concat(data1[i], data2[j], params), params);
        comp_j_i = this->compress(this->concat(data2[j], data1[i], params), params);
        K(i, j) = this->distance_measure(comp_single_1[i], comp_single_2[j], comp_i_j, comp_j_i);
      }
    }
  return K;
}

template<typename T>
double AbstractCompressionKernel<T>::NCD(const double &comp_x1, const double &comp_x2, const double &comp_x1_x2, const double &comp_x2_x1) const{
  return (((comp_x1_x2 + comp_x2_x1) / 2) - min(comp_x1, comp_x2)) / max(comp_x1, comp_x2);
}

template<typename T>
double AbstractCompressionKernel<T>::difference(const double &comp_x1, const double &comp_x2, const double &comp_x1_x2, const double &comp_x2_x1) const{
  return comp_x1 + comp_x2 - comp_x1_x2 - comp_x2_x1;
}


#endif /* BASEKERNEL */
