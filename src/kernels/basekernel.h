#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include<vector>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace py = pybind11;
using namespace std;
using namespace Eigen; 

template<typename T>
class Kernel{

public:

    virtual double dot(const T &x1, const T &x2, py::dict &params) const = 0;

    virtual MatrixXd computeKernelMatrix(const vector<T> &data, py::dict &params) const;

    static MatrixXd normaliseHilbert(const MatrixXd &K);

    static MatrixXd normaliseKrein(const MatrixXd &K);

};

template<typename T>
inline MatrixXd Kernel<T>::computeKernelMatrix(const vector<T> &data, py::dict &params) const {
    int n = data.size();
    MatrixXd K = MatrixXd::Zero(n, n);

    for (uint i = 0; i < n; i++){
        for (uint j = 0; j < n; j++){
            K(i, j) = this->dot(data[i], data[j], params);
            if (i != j){
                K(j, i) = K(i, j);
            }
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
    uint n = K.rows();
    MatrixXd Sigmap = MatrixXd::Zero(n, n);
    MatrixXd Sigmam = MatrixXd::Zero(n, n);
    SelfAdjointEigenSolver<MatrixXd> solver(K);
    cout << solver.eigenvalues() << endl;
    for(uint i = 0; i < n; i++){
        if(solver.eigenvalues()[i] >= 0.0){
            Sigmap(i, i) = solver.eigenvalues()[i];
        }
        else
        {
            Sigmam(i, i) = std::abs(solver.eigenvalues()[i]);
        }
    }
    return Kernel<T>::normaliseHilbert(solver.eigenvectors()* Sigmap * solver.eigenvectors().transpose()) - 
           Kernel<T>::normaliseHilbert(solver.eigenvectors()* Sigmam * solver.eigenvectors().transpose());
}

template<typename T>
class DiracKernel : public Kernel<T>{

public:

    virtual double dot(const T &x1, const T &x2, py::dict &params) const;

};

template<typename T>
inline double DiracKernel<T>::dot(const T &x1, const T &x2, py::dict &params) const{
    if (x1 == x2) {
        return 1.0;
    }
    return 0.0;
};