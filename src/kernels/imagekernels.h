#ifndef IMAGEKERNELS
#define IMAGEKERNELS

#include "basekernel.h"
#include "jpeglib.h"
#include <Dense>
#include <functional>
#include <jpeglib.h>

using namespace Eigen;
using namespace std;

typedef py::array_t<unsigned char, py::array::c_style> ImageMat;

//typedef Eigen::Matrix<unsigned char, Dynamic, Dynamic, RowMajor> ImageMat;

enum class CompressionMethod {
    Vertical, Horizontal, Both
};

enum class ImageType {
    BW, RGB
};

class ImageCompressionKernel : public AbstractCompressionKernel<ImageMat> {
private:
    function<double(const ImageMat &, const ImageMat &, const KernelParams &)> comp_method;

    function<ImageMat(const ImageMat &, const ImageMat &)> concatVertical;

    function<ImageMat(const ImageMat &, const ImageMat &)> concatHorizontal;

public:
    ImageCompressionKernel(const CompressionMethod &method, const ImageType &image_type);

    double dotHorizontal(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const;

    double dotVertical(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const;

    double dotBoth(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const;

    double dot(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const;

    ImageMat concatVerticalBW(const ImageMat &x, const ImageMat &y) const;

    ImageMat concatHorizontalBW(const ImageMat &x, const ImageMat &y) const;

    ImageMat concatVerticalRGB(const ImageMat &x, const ImageMat &y) const;

    ImageMat concatHorizontalRGB(const ImageMat &x, const ImageMat &y) const;
};

class JPEGCompressionKernel : public ImageCompressionKernel {
private:
    int input_components;
    J_COLOR_SPACE colour;
public:
    JPEGCompressionKernel(const CompressionMethod &method, const ImageType &image_type) : ImageCompressionKernel(method,
                                                                                                                 image_type) {
      switch (image_type) {
        case ImageType::BW:
          this->input_components = 1;
          this->colour = JCS_GRAYSCALE;
          break;

        case ImageType::RGB:
          this->input_components = 3;
          this->colour = JCS_RGB;
          break;

        default:
          cout << "Unknown image type!" << endl;
          break;
      }
    };

    double compress(const ImageMat &x, const int compressionlevel) const;
};

#endif /* IMAGEKERNELS */
