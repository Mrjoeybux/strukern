#include "./tensorkernel.h"
//
// Created by Mrjoeybux on 23/09/2020.
//

Tensor CubicTensorCompressionKernel::concat_dim1(const Tensor &x1, const Tensor &x2) const{
  int rows = x1.shape(0);
  int cols = x1.shape(1);
  int depth = x1.shape(2);
  int numels = rows * cols * depth;
  Tensor output({2 * rows, cols, depth});
  py::buffer_info x1buf = x1.request(), x2buf = x2.request(), outbuf = output.request();
  double *x1ptr = (double *) x1buf.ptr, *x2ptr = (double *) x2buf.ptr, *outptr = (double *) outbuf.ptr;
  for (uint i = 0; i < rows; i++) {
    for (uint j = 0; j < cols; j++) {
      for (uint k = 0; k < depth; k++) {
        outptr[k + depth * (j + cols * i)] = x1ptr[k + depth * (j + cols * i)];
        outptr[numels + k + depth * (j + cols * i)] = x2ptr[k + depth * (j + cols * i)];
      }
    }
  }
  return output;
}

Tensor CubicTensorCompressionKernel::concat_dim2(const Tensor &x1, const Tensor &x2) const{
  int rows = x1.shape(0);
  int cols = x1.shape(1);
  int depth = x1.shape(2);
  Tensor output({rows, 2 * cols, depth});
  py::buffer_info x1buf = x1.request(), x2buf = x2.request(), outbuf = output.request();
  double *x1ptr = (double *) x1buf.ptr, *x2ptr = (double *) x2buf.ptr, *outptr = (double *) outbuf.ptr;
  for (uint i = 0; i < rows; i++) {
    for (uint j = 0; j < 2*cols; j++) {
      for (uint k = 0; k < depth; k++) {
        if (j < cols) {
          outptr[k + depth * (j + 2 * cols * i)] = x1ptr[k + depth * (j + cols * i)];
        } else {
          outptr[k + depth * (j + 2 * cols * i)] = x2ptr[k + depth * ((j - cols) + cols * i)];
        }
      }
    }
  }
  return output;
}

Tensor CubicTensorCompressionKernel::concat_dim3(const Tensor &x1, const Tensor &x2) const{
  int rows = x1.shape(0);
  int cols = x1.shape(1);
  int depth = x1.shape(2);
  Tensor output({rows, cols, 2 * depth});
  py::buffer_info x1buf = x1.request(), x2buf = x2.request(), outbuf = output.request();
  double *x1ptr = (double *) x1buf.ptr, *x2ptr = (double *) x2buf.ptr, *outptr = (double *) outbuf.ptr;
  for (uint i = 0; i < rows; i++) {
    for (uint j = 0; j < cols; j++) {
      for (uint k = 0; k < 2 * depth; k++) {
        if (k < depth) {
          outptr[k + 2 * depth * (j + cols * i)] = x1ptr[k + depth * (j + cols * i)];
        } else {
          outptr[k + 2 * depth * (j + cols * i)] = x2ptr[(k - depth) + depth * (j + cols * i)];
        }
      }
    }
  }
  return output;
}

Tensor CubicTensorCompressionKernel::concat(const Tensor &x1, const Tensor &x2, const KernelParams &params) const {
  if(params.TensorConcatDim == 1){
    return this->concat_dim1(x1, x2);
  } else if(params.TensorConcatDim == 2){
    return this->concat_dim2(x1, x2);
  } else if(params.TensorConcatDim == 3){
    return this->concat_dim3(x1, x2);
  } else {
    throw domain_error("dimension = " + to_string(params.TensorConcatDim) + " is not a valid concatenation dimension!");
  }
}

double TThreshCompressionKernel::compress(const Tensor &x, const KernelParams &params) const{
  return tthresh_compress(x, params.TensorCompressionLevel);
}

/*double TThreshCompressionKernel::difference(const Tensor &x1, const Tensor &x2, const KernelParams &params) const {
  double comp_x1, comp_x2, comp_x1_x2, comp_x2_x1;
  comp_x1 = this->compress(x1, params.TensorCompressionLevel);
  comp_x2 = this->compress(x2, params.TensorCompressionLevel);
  comp_x1_x2 = this->compress(this->concat(x1, x2, params.TensorConcatDim), params.TensorCompressionLevel);
  comp_x2_x1 = this->compress(this->concat(x2, x1, params.TensorConcatDim), params.TensorCompressionLevel);

  return comp_x1 + comp_x2 - comp_x1_x2 - comp_x2_x1;
}

double TThreshCompressionKernel::NCD(const Tensor &x1, const Tensor &x2, const KernelParams &params) const {
  double comp_x1, comp_x2, comp_x1_x2, comp_x2_x1;
  comp_x1 = this->compress(x1, params.TensorCompressionLevel);
  comp_x2 = this->compress(x2, params.TensorCompressionLevel);
  comp_x1_x2 = this->compress(this->concat(x1, x2, params.TensorConcatDim), params.TensorCompressionLevel);
  comp_x2_x1 = this->compress(this->concat(x2, x1, params.TensorConcatDim), params.TensorCompressionLevel);
  return max(comp_x1_x2 - comp_x1, comp_x2_x1 - comp_x2) / max(comp_x1, comp_x2);
}*/


/*double TThreshCompressionKernel::dot(const Tensor &x1, const Tensor &x2, const KernelParams &params) const {
  return this->distance_measure(x1, x2, params);
}*/





