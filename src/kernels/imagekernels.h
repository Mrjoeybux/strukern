#ifndef IMAGEKERNELS
#define IMAGEKERNELS
#include "basekernel.h"
#include <Dense>
#include <functional>
using namespace Eigen;
using namespace std;

typedef py::array_t<unsigned char, py::array::c_style> ImageMat;

//typedef Eigen::Matrix<unsigned char, Dynamic, Dynamic, RowMajor> ImageMat;

enum class CompressionMethod { Vertical, Horizontal, Both };

class ImageCompressionKernel : public AbstractCompressionKernel<ImageMat> {
private:
  function<double(const ImageMat &, const ImageMat &, const KernelParams &)> comp_method;

public:
  ImageCompressionKernel(const CompressionMethod &method);

  double dotHorizontal(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const;

  double dotVertical(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const;

  double dotBoth(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const;

  double dot(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const;

  ImageMat concatVertical(const ImageMat &x, const ImageMat &y) const;

  ImageMat concatHorizontal(const ImageMat &x, const ImageMat &y) const;
};

class JPEGCompressionKernel : public ImageCompressionKernel {
public:
  JPEGCompressionKernel(const CompressionMethod &method) : ImageCompressionKernel(method){};

  double compress(const ImageMat &x, const int compressionlevel) const;
};

#endif /* IMAGEKERNELS */
