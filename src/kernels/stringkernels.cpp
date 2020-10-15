#include "./stringkernels.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <zlib.h>

string StringCompressionKernel::concat(const string &x1, const string &x2, const KernelParams &params) const{
  return x1 + x2;
}

double ZlibCompressionKernel::compress(const string &x, const KernelParams &params) const {
  z_stream zs; // z_stream is zlib's control structure
  memset(&zs, 0, sizeof(zs));

  if (deflateInit(&zs, params.StringCompressionLevel) != Z_OK)
    throw(std::runtime_error("deflateInit failed while compressing."));

  zs.next_in = (Bytef *)x.data();
  zs.avail_in = x.size(); // set the z_stream's input

  int ret;
  char outbuffer[32768];
  std::string outstring;

  // retrieve the compressed bytes blockwise
  do {
    zs.next_out = reinterpret_cast<Bytef *>(outbuffer);
    zs.avail_out = sizeof(outbuffer);

    ret = deflate(&zs, Z_FINISH);

    if (outstring.size() < zs.total_out) {
      // append the block to the output string
      outstring.append(outbuffer, zs.total_out - outstring.size());
    }
  } while (ret == Z_OK);

  deflateEnd(&zs);

  if (ret != Z_STREAM_END) { // an error occurred that was not EOF
    std::ostringstream oss;
    oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
    throw(std::runtime_error(oss.str()));
  }

  return outstring.size();
}

double PPMCompressionKernel::compress(const string &x, const KernelParams &params) const {
  istringstream is(x);
  ostringstream os;
  switch (params.StringCompressionLevel) {
    case 1: {
      o1 compressor;
      compressor.compress(is, os);
    }break;
    case 2:{
      o2 compressor;
      compressor.compress(is, os);
    }break;
    case 3:{
      o3 compressor;
      compressor.compress(is, os);
    }break;
    case 4:{
      o4 compressor;
      compressor.compress(is, os);
    }break;
    case 5:{
      o5 compressor;
      compressor.compress(is, os);
    }break;
  }
  return os.tellp();
}

double LocalityImprovedKernel::dot(const string &x1, const string &x2, const KernelParams &params) const {
    int L = params.LocalityImproved.at("sub_window_length");
    int d1 = params.LocalityImproved.at("d1");
    int d2 = params.LocalityImproved.at("d2");
    double sum = 0;
    for(uint i = L; i < x1.size() - L; i++){
        sum += this->sub_window(x1.substr(i - L, 2*L + 1), x2.substr(i - L, 2*L + 1), d1);
    }
    return pow(sum, d2);
}

double LocalityImprovedKernel::sub_window(const string &x1_substr, const string &x2_substr, const int &d1) const {
    double sum = 0;
    for(uint i = 0; i < x1_substr.size(); i++){
        if(x1_substr[i] == x2_substr[i]){
            sum += 1;
        }
    }
    return pow(sum, d1);
}
