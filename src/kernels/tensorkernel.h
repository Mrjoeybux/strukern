//
// Created by Mrjoeybux on 23/09/2020.
//
#include "../tthresh_usage/my_compress.hpp"
#include "../kernels/basekernel.h"
#ifndef STRUKERN_TENSORKERNEL_H
#define STRUKERN_TENSORKERNEL_H


class CubicTensorCompressionKernel : public AbstractCompressionKernel<Tensor> {

public:
    using AbstractCompressionKernel<Tensor>::AbstractCompressionKernel;
    // dim1 = rows
    Tensor concat_dim1(const Tensor &x1, const Tensor &x2) const;
    // dim2 = cols
    Tensor concat_dim2(const Tensor &x1, const Tensor &x2) const;
    // dim3 = cols
    Tensor concat_dim3(const Tensor &x1, const Tensor &x2) const;

    Tensor concat(const Tensor &x1, const Tensor &x2, const KernelParams &params) const;

    virtual double compress(const Tensor &x1, const KernelParams &params) const = 0;

};

class TThreshCompressionKernel : public CubicTensorCompressionKernel{

public:

    using CubicTensorCompressionKernel::CubicTensorCompressionKernel;

    double compress(const Tensor &x, const KernelParams &params) const;
};

#endif //STRUKERN_TENSORKERNEL_H
