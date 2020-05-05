#include "imagekernels.h"
#include "jpeglib.h"
#include <iostream>
#include <pybind11/pybind11.h>

namespace py = pybind11;
//#include <Dense>

ImageCompressionKernel::ImageCompressionKernel(const CompressionMethod &method) {
  switch (method) {
  case CompressionMethod::Vertical:
    this->comp_method = [this](const ImageMat &x1, const ImageMat &x2, const KernelParams &params) {
      return dotVertical(x1, x2, params);
    };
    break;

  case CompressionMethod::Horizontal:
    this->comp_method = [this](const ImageMat &x1, const ImageMat &x2, const KernelParams &params) {
      return dotHorizontal(x1, x2, params);
    };
    break;

  case CompressionMethod::Both:
    this->comp_method = [this](const ImageMat &x1, const ImageMat &x2, const KernelParams &params) {
      return dotBoth(x1, x2, params);
    };
    break;

  default:
    cout << "Unknown compression method!" << endl;
    break;
  }
}

double ImageCompressionKernel::dotVertical(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const {
    ImageMat x1abovex2 = this->concatVertical(x1, x2);
    ImageMat x2abovex1 = this->concatVertical(x2, x1);
  return this->compress(x1, params.JPEGCompressionQuality) + this->compress(x2, params.JPEGCompressionQuality) -
         this->compress(x1abovex2, params.JPEGCompressionQuality) - this->compress(x2abovex1, params.JPEGCompressionQuality);
};

double ImageCompressionKernel::dotHorizontal(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const {
    ImageMat x1nexttox2 = this->concatHorizontal(x1, x2);
    ImageMat x2nexttox1 = this->concatHorizontal(x2, x1);
  return this->compress(x1, params.JPEGCompressionQuality) + this->compress(x2, params.JPEGCompressionQuality) -
         this->compress(x1nexttox2, params.JPEGCompressionQuality) - this->compress(x2nexttox1, params.JPEGCompressionQuality);
};

double ImageCompressionKernel::dotBoth(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const {
    ImageMat x1abovex2 = this->concatVertical(x1, x2);
    ImageMat x2abovex1 = this->concatVertical(x2, x1);
    ImageMat x1nexttox2 = this->concatHorizontal(x1, x2);
    ImageMat x2nexttox1 = this->concatHorizontal(x2, x1);
  return 2 * (this->compress(x1, params.JPEGCompressionQuality) + this->compress(x2, params.JPEGCompressionQuality)) -
         this->compress(x1nexttox2, params.JPEGCompressionQuality) - this->compress(x2nexttox1, params.JPEGCompressionQuality) -
         this->compress(x1abovex2, params.JPEGCompressionQuality) - this->compress(x2abovex1, params.JPEGCompressionQuality);
}

double ImageCompressionKernel::dot(const ImageMat &x1, const ImageMat &x2, const KernelParams &params) const {
  return this->comp_method(x1, x2, params);
}

ImageMat ImageCompressionKernel::concatVertical(const ImageMat &x, const ImageMat &y) const {
    int rows = x.shape(0);
    int cols = x.shape(1);
    int depth =  x.shape(2);
    int numels = rows*cols*depth;
    py::array_t<unsigned char> z({2*rows, cols, depth}, {cols*depth, depth, 1});
    py::buffer_info xbuf = x.request(), ybuf = y.request(), zbuf = z.request();
    unsigned char *xptr = (unsigned char *) xbuf.ptr, *yptr = (unsigned char *) ybuf.ptr, *zptr = (unsigned char *) zbuf.ptr;
    for(uint i = 0; i < rows; i++){
        for(uint j = 0; j < cols; j++){
            for(uint k = 0; k < depth; k++){
                zptr[k + depth*(j + cols*i)] = xptr[k + depth*(j + cols*i)];
                zptr[numels + k + depth*(j + cols*i)] = yptr[k + depth*(j + cols*i)];
            }
        }
    }
    return z;
}

ImageMat ImageCompressionKernel::concatHorizontal(const ImageMat &x, const ImageMat &y) const {
    int rows = x.shape(0);
    int cols = x.shape(1);
    int depth =  x.shape(2);
    py::array_t<unsigned char> z({rows, 2*cols, depth}, {2*cols*depth, depth, 1});
    py::buffer_info xbuf = x.request(), ybuf = y.request(), zbuf = z.request();
    unsigned char *xptr = (unsigned char *) xbuf.ptr, *yptr = (unsigned char *) ybuf.ptr, *zptr = (unsigned char *) zbuf.ptr;
    for(uint i = 0; i < rows; i++){
        for(uint j = 0; j < 2*cols; j++){
            for(uint k = 0; k < depth; k++){
                if(j < cols){
                    zptr[k + depth*(j + 2*cols*i)] = xptr[k + depth*(j + cols*i)];
                }
                else{
                    zptr[k + depth*(j + 2*cols*i)] = yptr[k + depth*((j - cols) + cols*i)];
                }

            }
        }
    }
    return z;
}

double JPEGCompressionKernel::compress(const ImageMat &x, const int compressionlevel) const {
  unsigned char *raw_image = (unsigned char *)x.data();
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW row_pointer[1];
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  unsigned char *buffer = NULL;
  unsigned long size;
  jpeg_mem_dest(&cinfo, &buffer, &size);
    cinfo.image_height = x.shape(0);
    cinfo.image_width = x.shape(1);
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;

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

/*double compressColour(const py::array_t<unsigned char> &x, const int compressionlevel){
    unsigned char *raw_image = (unsigned char *)x.data();
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    JSAMPROW row_pointer[1];
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    unsigned char *buffer = NULL;
    unsigned long size;
    jpeg_mem_dest(&cinfo, &buffer, &size);
    cinfo.image_height = x.shape(0);
    cinfo.image_width = x.shape(1);
    // TODO: MAKE USABLE FOR COLOR IMAGES
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;

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

py::array_t<unsigned char, py::array::c_style> concatVertical(const py::array_t<unsigned char, py::array::c_style> &x, const py::array_t<unsigned char, py::array::c_style> &y) {
    int rows = x.shape(0);
    int cols = x.shape(1);
    int depth =  x.shape(2);
    int numels = rows*cols*depth;
    py::array_t<unsigned char> z({2*rows, cols, depth}, {cols*depth, depth, 1});
    py::buffer_info xbuf = x.request(), ybuf = y.request(), zbuf = z.request();
    unsigned char *xptr = (unsigned char *) xbuf.ptr, *yptr = (unsigned char *) ybuf.ptr, *zptr = (unsigned char *) zbuf.ptr;
    for(uint i = 0; i < rows; i++){
        for(uint j = 0; j < cols; j++){
            for(uint k = 0; k < depth; k++){
                zptr[k + depth*(j + cols*i)] = xptr[k + depth*(j + cols*i)];
                zptr[numels + k + depth*(j + cols*i)] = yptr[k + depth*(j + cols*i)];
            }
        }
    }
    return z;
}

py::array_t<unsigned char, py::array::c_style> concatHorizontal(const py::array_t<unsigned char, py::array::c_style> &x, const
py::array_t<unsigned char, py::array::c_style> &y) {
    int rows = x.shape(0);
    int cols = x.shape(1);
    int depth =  x.shape(2);
    py::array_t<unsigned char> z({rows, 2*cols, depth}, {2*cols*depth, depth, 1});
    py::buffer_info xbuf = x.request(), ybuf = y.request(), zbuf = z.request();
    unsigned char *xptr = (unsigned char *) xbuf.ptr, *yptr = (unsigned char *) ybuf.ptr, *zptr = (unsigned char *) zbuf.ptr;
    int xcolcounter = 0, ycolcounter = 0;
    for(uint i = 0; i < rows; i++){
        for(uint j = 0; j < 2*cols; j++){
            for(uint k = 0; k < depth; k++){
                if(j < cols){
                    zptr[k + depth*(j + 2*cols*i)] = xptr[k + depth*(j + cols*i)];
                }
                else{
                    zptr[k + depth*(j + 2*cols*i)] = yptr[k + depth*((j - cols) + cols*i)];
                }

            }
        }
    }
    return z;
}*/