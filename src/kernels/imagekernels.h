#ifndef IMAGEKERNELS
#define IMAGEKERNELS
#include "basekernel.h"
#include <Dense>
using namespace Eigen;

typedef Matrix<unsigned char, Dynamic, Dynamic, RowMajor> JPEGImageMat;

class hImageCompressionKernel : public AbstractCompressionKernel<JPEGImageMat, int> {

public:
  double dot(const JPEGImageMat &x1, const JPEGImageMat &x2, const unordered_map<string, int> &params) const override;
};

class vImageCompressionKernel : public AbstractCompressionKernel<JPEGImageMat, int> {

public:
  double dot(const JPEGImageMat &x1, const JPEGImageMat &x2, const unordered_map<string, int> &params) const;
};

class hvImageCompressionKernel : public AbstractCompressionKernel<JPEGImageMat, int> {

public:
  double dot(const JPEGImageMat &x1, const JPEGImageMat &x2, const unordered_map<string, int> &params) const override;
};

class vJPEGCompressionKernel : public vImageCompressionKernel {
public:
  double compress(const JPEGImageMat &x, const int compressionlevel) const;
};

#endif /* IMAGEKERNELS */
