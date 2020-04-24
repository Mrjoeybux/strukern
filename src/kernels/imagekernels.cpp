#include "imagekernels.h"
#include "jpeglib.h"
#include <iostream>
//#include <Dense>

double vImageCompressionKernel::dot(const JPEGImageMat &x1, const JPEGImageMat &x2, const unordered_map<string, int> &params) const {
  JPEGImageMat x1x2(x1.rows(), 2 * x1.cols()), x2x1(x1.rows(), 2 * x1.cols());
  x1x2 << x1, x2;
  x2x1 << x2, x1;
  return this->compress(x1, params.at("compressionlevel")) + this->compress(x2, params.at("compressionlevel")) -
         this->compress(x1x2, params.at("compressionlevel")) - this->compress(x2x1, params.at("compressionlevel"));
};

double vJPEGCompressionKernel::compress(const JPEGImageMat &x, const int compressionlevel) const {
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