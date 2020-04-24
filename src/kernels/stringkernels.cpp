#include "./stringkernels.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <zlib.h>

double StringCompressionKernel::dot(const string &x1, const string &x2, const unordered_map<string, int> &params) const {
  return this->compress(x1, params.at("compressionlevel")) + this->compress(x2, params.at("compressionlevel")) -
         this->compress(x1 + x2, params.at("compressionlevel")) - this->compress(x2 + x1, params.at("compressionlevel"));
}

double ZlibCompressionKernel::compress(const string &x, int compressionlevel) const {
  z_stream zs; // z_stream is zlib's control structure
  memset(&zs, 0, sizeof(zs));

  if (deflateInit(&zs, compressionlevel) != Z_OK)
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