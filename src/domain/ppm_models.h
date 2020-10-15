//
// Created by Mrjoeybux on 29/07/2020.
//

#ifndef STRUKERN_PPM_MODELS_H
#define STRUKERN_PPM_MODELS_H
#include <dlib/entropy_decoder.h>
#include <dlib/entropy_decoder_model.h>
#include <dlib/entropy_encoder.h>
#include <dlib/entropy_encoder_model/entropy_encoder_model_kernel_5.h>
#include <dlib/compress_stream/compress_stream_kernel_1.h>
#include <dlib/crc32/crc32_kernel_1.h>
using namespace dlib;

/*enum class PPM_ORDER{
    ONE, TWO, THREE, FOUR, FIVE
};*/

template <unsigned long alphabet_size, typename entropy_encoder>
class my_entropy_encoder_model
{
    my_entropy_encoder_model() = default;

public:
    // kernel_5
    typedef     entropy_encoder_model_kernel_5<alphabet_size,entropy_encoder,25000,1>
        kernel_5_1;
    typedef     entropy_encoder_model_kernel_5<alphabet_size,entropy_encoder,50000,2>
        kernel_5_2;
    typedef     entropy_encoder_model_kernel_5<alphabet_size,entropy_encoder,100000,3>
        kernel_5_3;
    typedef     entropy_encoder_model_kernel_5<alphabet_size,entropy_encoder,200000,4>
        kernel_5_4;
    typedef     entropy_encoder_model_kernel_5<alphabet_size,entropy_encoder,1000000,5>
        kernel_5_5;


};

class my_compress_stream{

    typedef entropy_decoder_model<257,entropy_decoder::kernel_2a>::kernel_1b fcd1;

public:

    typedef compress_stream_kernel_1<
    my_entropy_encoder_model<257, entropy_encoder::kernel_1a>::kernel_5_1,
    fcd1,
    crc32::kernel_1a> order_1;

    typedef compress_stream_kernel_1<
    my_entropy_encoder_model<257, entropy_encoder::kernel_1a>::kernel_5_2,
    fcd1,
    crc32::kernel_1a> order_2;

    typedef compress_stream_kernel_1<
    my_entropy_encoder_model<257, entropy_encoder::kernel_1a>::kernel_5_3,
    fcd1,
    crc32::kernel_1a> order_3;

    typedef compress_stream_kernel_1<
    my_entropy_encoder_model<257, entropy_encoder::kernel_1a>::kernel_5_4,
    fcd1,
    crc32::kernel_1a> order_4;

    typedef compress_stream_kernel_1<
    my_entropy_encoder_model<257, entropy_encoder::kernel_1a>::kernel_5_5,
    fcd1,
    crc32::kernel_1a> order_5;
};



#endif //STRUKERN_PPM_MODELS_H
