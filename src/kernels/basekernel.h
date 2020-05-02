#ifndef BASEKERNEL
#define BASEKERNEL
#include "../domain/common.h"
#include <iostream>
#include <Core>
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
#include <lsap.h>
using namespace std;
using namespace Eigen;

struct KernelParams {
  int ZlibCompressionLevel;
  int JPEGCompressionQuality;
  unordered_map<string, LabelPairMap> MABTS;
};

#define NUMERICAL_PRECISION 1e-12
#define VERYLARGECOST 1e12

enum class KernelType { Hilbert, Krein };

enum class MultiInstanceMethod {SUM, MAX};

template <typename T> class Kernel {

public:
  virtual double dot(const T &x1, const T &x2, const KernelParams &params) const = 0;

  virtual MatrixXd computeKernelMatrix(const vector<T> &data, const KernelParams &params) const;

  static MatrixXd normaliseHilbert(const MatrixXd &K);

  static MatrixXd normaliseKrein(const MatrixXd &K);

  /*virtual VectorXd quadratic_form_kmat_gradient(const VectorXd &u, const MatrixXd &kmat, const VectorXd &v, const vector<T> &X,
                                                const KernelParams &params) const;

  virtual VectorXd quadratic_form_kmat_block_gradient(const VectorXd &u, const MatrixXd &kmat_block, const VectorXd &v, const vector<T> &X,
                                                      const vector<int> &support_args, const vector<int> &block_args, const MatrixXd &kmat,
                                                      const KernelParams &params) const;
  KernelType ktype;*/
};

template <typename T> inline MatrixXd Kernel<T>::computeKernelMatrix(const vector<T> &data, const KernelParams &params) const {
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

template <typename T> inline MatrixXd Kernel<T>::normaliseHilbert(const MatrixXd &K) {
    VectorXd kvec = K.diagonal().cwiseInverse().cwiseSqrt();
    return K.cwiseProduct(kvec * kvec.transpose());
}

template <typename T> inline MatrixXd Kernel<T>::normaliseKrein(const MatrixXd &K) {
    SelfAdjointEigenSolver<MatrixXd> solver(K);
    MatrixXd Sigmap = solver.eigenvalues().cwiseMax(0).asDiagonal();
    MatrixXd Sigmam = solver.eigenvalues().cwiseMin(0).cwiseAbs().asDiagonal();
    return Kernel<T>::normaliseHilbert(solver.eigenvectors() * Sigmap * solver.eigenvectors().transpose()) -
           Kernel<T>::normaliseHilbert(solver.eigenvectors() * Sigmam * solver.eigenvectors().transpose());
}

template <typename T> class MultiInstanceKernel : public Kernel<vector<T>> {
private:
    const Kernel<T> *basekernel;

    function<double(const vector<T> &, const vector<T> &, const KernelParams &)> dot_func;

public:
    MultiInstanceKernel(const Kernel<T> *basekernel, const MultiInstanceMethod &method);

    double dot(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const;

    double dot_sum(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const;

    double max_sum(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const;

    MatrixXd computeKernelMatrix(const vector<vector<T>> &data, const KernelParams &params) const;
};

template<typename T>
MultiInstanceKernel<T>::MultiInstanceKernel(const Kernel<T> *basekernel, const MultiInstanceMethod &method) {
    this->basekernel = basekernel;
    switch (method) {
        case MultiInstanceMethod::SUM:
            this->dot_func = [this](const vector<T> &x1, const vector<T> &x2, const KernelParams &params) {
                return this->dot_sum(x1, x2, params);
            };
            break;
        case MultiInstanceMethod::MAX:
            this->dot_func = [this](const vector<T> &x1, const vector<T> &x2, const KernelParams &params) {
                return this->max_sum(x1, x2, params);
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
double MultiInstanceKernel<T>::dot_sum(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const {
    double sum = 0.0;
    for(uint i = 0; i < x1.size(); i++){
        for(uint j = 0; j < x2.size(); j++){
            sum += this->basekernel->dot(x1[i], x2[j], params);
        }
    }
    return sum;
}

template<typename T>
double MultiInstanceKernel<T>::max_sum(const vector<T> &x1, const vector<T> &x2, const KernelParams &params) const {
    size_t nRows = x1.size(), nCols = x2.size();
    double *C = new double[nRows*nCols];
    size_t *rho = new size_t[nRows];
    double *u = new double[nRows], *v = new double[nCols];
    for(uint col = 0; col < nCols; col++){
        for(uint row = 0; row < nRows; row++){
            C[nRows*col + row] = this->basekernel->dot(x1[row], x2[col], params);
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



template <typename T> class AbstractCompressionKernel : public Kernel<T> {

public:
    //KernelType ktype = KernelType::Krein;

    virtual double compress(const T &x, const int compressionlevel) const = 0;
};


#endif /* BASEKERNEL */
