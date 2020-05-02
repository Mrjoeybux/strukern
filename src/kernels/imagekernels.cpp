#include "imagekernels.h"
#include "jpeglib.h"
#include <iostream>
//#include <Dense>

ImageCompressionKernel::ImageCompressionKernel(const CompressionMethod &method) {
  switch (method) {
  case CompressionMethod::Vertical:
    this->comp_method = [this](const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) {
      return dotVertical(x1, x2, params);
    };
    break;

  case CompressionMethod::Horizontal:
    this->comp_method = [this](const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) {
      return dotHorizontal(x1, x2, params);
    };
    break;

  case CompressionMethod::Both:
    this->comp_method = [this](const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) {
      return dotBoth(x1, x2, params);
    };
    break;

  default:
    cout << "Unknown compression method!" << endl;
    break;
  }
}

double ImageCompressionKernel::dotVertical(const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) const {
  JPEGImageMat x1abovex2(2 * x1.rows(), x1.cols()), x2abovex1(2 * x1.rows(), x1.cols());
  x1abovex2 << x1, x2;
  x2abovex1 << x2, x1;
  return this->compress(x1, params.JPEGCompressionQuality) + this->compress(x2, params.JPEGCompressionQuality) -
         this->compress(x1abovex2, params.JPEGCompressionQuality) - this->compress(x2abovex1, params.JPEGCompressionQuality);
};

double ImageCompressionKernel::dotHorizontal(const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) const {
  JPEGImageMat x1nexttox2(x1.rows(), 2 * x1.cols()), x2nexttox1(x1.rows(), 2 * x1.cols());
  x1nexttox2 << x1, x2;
  x2nexttox1 << x2, x1;
  return this->compress(x1, params.JPEGCompressionQuality) + this->compress(x2, params.JPEGCompressionQuality) -
         this->compress(x1nexttox2, params.JPEGCompressionQuality) - this->compress(x2nexttox1, params.JPEGCompressionQuality);
};

double ImageCompressionKernel::dotBoth(const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) const {
  JPEGImageMat x1abovex2(2 * x1.rows(), x1.cols()), x2abovex1(2 * x1.rows(), x1.cols()), x1nexttox2(x1.rows(), 2 * x1.cols()),
      x2nexttox1(x1.rows(), 2 * x1.cols());
  x1abovex2 << x1, x2;
  x2abovex1 << x2, x1;
  x1nexttox2 << x1, x2;
  x2nexttox1 << x2, x1;
  return 2 * (this->compress(x1, params.JPEGCompressionQuality) + this->compress(x2, params.JPEGCompressionQuality)) -
         this->compress(x1nexttox2, params.JPEGCompressionQuality) - this->compress(x2nexttox1, params.JPEGCompressionQuality) -
         this->compress(x1abovex2, params.JPEGCompressionQuality) - this->compress(x2abovex1, params.JPEGCompressionQuality);
}

double ImageCompressionKernel::dot(const JPEGImageMat &x1, const JPEGImageMat &x2, const KernelParams &params) const {
  return this->comp_method(x1, x2, params);
}

double JPEGCompressionKernel::compress(const JPEGImageMat &x, const int compressionlevel) const {
  unsigned char *raw_image = (unsigned char *)x.data();
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW row_pointer[1];
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  unsigned char *buffer = NULL;
  unsigned long size;
  jpeg_mem_dest(&cinfo, &buffer, &size);
  cinfo.image_width = x.cols();
  cinfo.image_height = x.rows();
  // TODO: MAKE USABLE FOR COLOR IMAGES
  cinfo.input_components = 1;
  cinfo.in_color_space = JCS_GRAYSCALE;

  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, compressionlevel, TRUE);
  jpeg_start_compress(&cinfo, TRUE);

  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer[0] = &raw_image[cinfo.next_scanline * cinfo.image_width * cinfo.input_components];
    jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);
  free(buffer);
  return size;
}