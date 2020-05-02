#ifndef IMAGEKERNELS
#define IMAGEKERNELS
#include "basekernel.h"
#include <Dense>
#include <functional>
using namespace Eigen;
using namespace std;

typedef Eigen::Matrix<unsigned char, Dynamic, Dynamic, RowMajor> JPEGImageMat;

enum class CompressionMethod { Vertical, Horizontal, Both };

class ImageCompressionKernel : public AbstractCompressionKernel<JPEGImageMat> {
private:
  function<double(const JPEGImageMat &, const JPEGImageMat &, const KernelParams &)> comp_method;

public:
  ImageCompressionKernel(const CompressionMethod &method);

  double dotHorizontal(const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) const;

  double dotVertical(const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) const;

  double dotBoth(const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) const;

  double dot(const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) const;
};

class JPEGCompressionKernel : public ImageCompressionKernel {
public:
  JPEGCompressionKernel(const CompressionMethod &method) : ImageCompressionKernel(method){};

  double compress(const JPEGImageMat &x, const int compressionlevel) const;
};

#endif /* IMAGEKERNELS */
