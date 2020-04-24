
#include "kernels/gedlib_test.h"
#include "kernels/imagekernels.h"
#include "kernels/jpeg_test.h"
#include <Dense>
#include <iostream>

int main() {
  test();
  // char infilename[] = "/home/mrjoeybux/Firefox_wallpaper.jpeg", outfilename[] = "out.jpeg";
  // char *in = infilename, *out = outfilename;
  // unsigned char *image = read_jpeg_file(in);
  // double *double_array = (double *)image;
  // Map<MatrixXd> x(double_array, 2160, 3840);
  // cout << x(0, 0);
  std::ifstream infile("/home/mrjoeybux/coding/strukern/src/mnist_image_converted.dat");
  std::string line, newline;
  MatrixXi x(28, 28);
  int i = 0, j = 0;
  while (std::getline(infile, line)) {
    if (!line.empty() && line[line.length() - 1] == '\n') {
      line.erase(line.length() - 1);
    }
    x(i, j) = std::stoi(line);
    if (j == 27) {
      j = 0;
      i++;
    } else {
      j++;
    }
  }
  vJPEGCompressionKernel kf;
  JPEGImageMat y = x.cast<unsigned char>();
  int size = kf.compress(y, 1);
  cout << "Compression size: " << size << "\n";
  // write_jpeg_file(out);
  return 0;
}