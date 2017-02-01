#ifndef _TYPE_DEF_H_
#define _TYPE_DEF_H_

#include <armadillo>
#include <stdint.h>

namespace mat_t {
    typedef unsigned char bit;
    typedef arma::uword   uint;
    typedef arma::sword   sint;
    typedef float         val;
    typedef bool          logic;
    const   mat_t::val    eps   = arma::datum::eps;
    const   mat_t::bit    bit_0 = 0b0;
    const   mat_t::bit    bit_1 = 0b1;
    const   mat_t::logic  no    = false;
    const   mat_t::logic  yes   = true;
}

namespace aud_t {
    typedef int8_t        int8;
    typedef int16_t       int16;
    typedef int32_t       int32;
    typedef unsigned long uint;
    typedef long          sint;
    typedef double        fpt;
    typedef bool          logic;
    const mat_t::logic    no  = false;
    const mat_t::logic    yes = true;
}

#endif
