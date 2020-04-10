#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include <vector>
#include <string>

using namespace Eigen;
namespace py = pybind11;
using namespace std;
using namespace pybind11::literals;

template<typename T>
class KernelMatrix{
    public:
    MatrixXd kmat;
    MatrixXd &compute_kernel_matrix1(){
        MatrixXd m(2,2);
        m(0,0) = 3;
        m(1,0) = 2.5;
        m(0,1) = -1;
        m(1,1) = m(1,0) + m(0,1);
        kmat = m;
        return kmat;
    }

    MatrixXd &compute_kernel_matrix2(){
        MatrixXd m(2,2);
        m(0,0) = -3;
        m(1,0) = -2.5;
        m(0,1) = 1;
        m(1,1) = m(1,0) + m(0,1);
        ///kmat = m;
        return m;
    }

    void print_dict(py::dict dict){
        for (auto item : dict)
            cout << "key: " << item.first << ", value = " << item.second << endl;
    }

    void print_list(vector<T> &list){
        for (auto item : list)
            cout << item << endl;
    }
};


int main(){
    return 0;
}

PYBIND11_MODULE(kernelmat, m) {
    py::class_<KernelMatrix<string>>(m, "KernelMatrix")
        .def(py::init<>())
        .def("compute_kernel_matrix1", &KernelMatrix<string>::compute_kernel_matrix1)
        .def("print_dict", &KernelMatrix<string>::print_dict)
        .def("print_list", &KernelMatrix<string>::print_list);
};