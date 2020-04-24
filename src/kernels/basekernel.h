#ifndef BASEKERNEL
#define BASEKERNEL
#include <Core>
#include <Dense>
#include <Eigenvalues>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
using namespace Eigen;

#define NUMERICAL_PRECISION 1e-12

enum class KernelType { Hilbert, Krein };

template <typename T, typename bandwidthType> class Kernel {

public:
  virtual double dot(const T &x1, const T &x2, const unordered_map<string, bandwidthType> &params) const = 0;

  virtual MatrixXd computeKernelMatrix(const vector<T> &data, const unordered_map<string, bandwidthType> &params) const;

  static MatrixXd normaliseHilbert(const MatrixXd &K);

  static MatrixXd normaliseKrein(const MatrixXd &K);

  /*virtual VectorXd quadratic_form_kmat_gradient(const VectorXd &u, const MatrixXd &kmat, const VectorXd &v, const vector<T> &X,
                                                const unordered_map<string, bandwidthType> &params) const;

  virtual VectorXd quadratic_form_kmat_block_gradient(const VectorXd &u, const MatrixXd &kmat_block, const VectorXd &v, const vector<T> &X,
                                                      const vector<int> &support_args, const vector<int> &block_args, const MatrixXd &kmat,
                                                      const unordered_map<string, bandwidthType> &params) const;
  KernelType ktype;*/
};

template <typename T, typename bandwidthType> class AbstractCompressionKernel : public Kernel<T, bandwidthType> {

public:
  KernelType ktype = KernelType::Krein;

  virtual double compress(const T &x, const int compressionlevel) const = 0;
};

template <typename T, typename bandwidthType>
inline MatrixXd Kernel<T, bandwidthType>::computeKernelMatrix(const vector<T> &data,
                                                              const unordered_map<string, bandwidthType> &params) const {
  int n = data.size();
  MatrixXd K = MatrixXd::Zero(n, n);

  for (uint i = 0; i < n; i++) {
    for (uint j = 0; j < n; j++) {
      K(i, j) = this->dot(data[i], data[j], params);
      if (i != j) {
        K(j, i) = K(i, j);
      }
    }
  }
  return K;
}

template <typename T, typename bandwidthType> inline MatrixXd Kernel<T, bandwidthType>::normaliseHilbert(const MatrixXd &K) {
  VectorXd kvec = K.diagonal().cwiseInverse().cwiseSqrt();
  return K.cwiseProduct(kvec * kvec.transpose());
}

template <typename T, typename bandwidthType> inline MatrixXd Kernel<T, bandwidthType>::normaliseKrein(const MatrixXd &K) {
  SelfAdjointEigenSolver<MatrixXd> solver(K);
  MatrixXd Sigmap = solver.eigenvalues().cwiseMax(0).asDiagonal();
  MatrixXd Sigmam = solver.eigenvalues().cwiseMin(0).cwiseAbs().asDiagonal();
  return Kernel<T, bandwidthType>::normaliseHilbert(solver.eigenvectors() * Sigmap * solver.eigenvectors().transpose()) -
         Kernel<T, bandwidthType>::normaliseHilbert(solver.eigenvectors() * Sigmam * solver.eigenvectors().transpose());
}
#endif /* BASEKERNEL */
