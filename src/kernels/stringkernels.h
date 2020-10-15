#ifndef STRINGKERNELS
#define STRINGKERNELS

#include "./basekernel.h"
#include <string>
#include "../domain/ppm_models.h"

using namespace std;

typedef my_compress_stream::order_1 o1;
typedef my_compress_stream::order_2 o2;
typedef my_compress_stream::order_3 o3;
typedef my_compress_stream::order_4 o4;
typedef my_compress_stream::order_5 o5;


class StringCompressionKernel : public AbstractCompressionKernel<string> {

public:
    using AbstractCompressionKernel<string>::AbstractCompressionKernel;

    string concat(const string &x1, const string &x2, const KernelParams &params) const override;

    virtual double compress(const string &x, const KernelParams &params) const = 0;

};


class ZlibCompressionKernel : public StringCompressionKernel {

public:
    using StringCompressionKernel::StringCompressionKernel;

    double compress(const string &x, const KernelParams &params) const override;
};

class PPMCompressionKernel : public StringCompressionKernel {

public:
    using StringCompressionKernel::StringCompressionKernel;

    double compress(const string &x, const KernelParams &params) const override;

};

class LocalityImprovedKernel : public Kernel<string> {
public:
    double sub_window(const string &x1_substr, const string &x2_substr, const int &d1) const;

    double dot(const string &x1, const string &x2, const KernelParams &params) const override;
};

#endif /* STRINGKERNELS */