#ifndef STRINGKERNELS
#define STRINGKERNELS
#include "./basekernel.h"
#include <string>
using namespace std;

class StringCompressionKernel : public AbstractCompressionKernel<string, int> {

public:
  double dot(const string &x1, const string &x2, const unordered_map<string, int> &params) const override;
};

class ZlibCompressionKernel : public StringCompressionKernel {

public:
  double compress(const string &x, const int compressionlevel) const override;
};

#endif /* STRINGKERNELS */