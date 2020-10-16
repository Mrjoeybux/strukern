/*
 * Copyright (c) 2016-2017, Rafael Ballester-Ripoll
 *                          (Visualization and MultiMedia Lab, University of Zurich),
 *                          rballester@ifi.uzh.ch
 *
 * Licensed under the LGPLv3.0 (https://github.com/rballester/tthresh/blob/master/LICENSE)
 */

#ifndef __IO_HPP__
#define __IO_HPP__

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "tthresh.hpp"
//#include "zlib.h"

// Avoids corruption of the input and output data on Windows/MS-DOS systems
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

struct {
    uint64_t wbytes;
    int8_t wbit;
    FILE *file; // File handle to read/write from/to
    uint8_t inout[CHUNK]; // Buffer to write the results of inflation/deflation
    uint8_t buf[CHUNK]; // Buffer used for the read/write operations
    int32_t bufstart = 0;
    int32_t bufend = 0;
    size_t total_written_bytes = 0; // Used to compute the final file size
} zs; // Read/write state for zlib interfacing

/*********/
// Writing
/*********/

// Call open_wbit() before write_bits()
// If write_bits() has been called, call close_wbit() before write_stream()

void write_stream(unsigned char *buf, size_t bytes_to_write)
{
    zs.total_written_bytes += bytes_to_write;
}

void open_wbit() {
    zs.wbytes = 0;
    zs.wbit = 63;
}

// Assumption: to_write <= 64
void write_bits(uint64_t bits, char to_write) {
    if (to_write <= zs.wbit+1) {
        zs.wbytes |= bits << (zs.wbit+1-to_write);
        zs.wbit -= to_write;
    }
    else {
        if (zs.wbit > -1)
            zs.wbytes |= bits >> (to_write-(zs.wbit+1));
        write_stream(reinterpret_cast<unsigned char *> (&zs.wbytes), sizeof(zs.wbytes));
        to_write -= zs.wbit+1;
        zs.wbytes = 0;
        zs.wbytes |= bits << (64-to_write);
        zs.wbit = 63-to_write;
    }
}

void close_wbit() {
    // Write any reamining bits
    if (zs.wbit < 63)
        write_stream(reinterpret_cast < unsigned char *> (&zs.wbytes), sizeof(zs.wbytes));
}


#endif // IO_HPP
