#ifndef _TYPE_DEF_H_
#define _TYPE_DEF_H_

#include <armadillo>

namespace type {
    typedef unsigned char bitval;
    typedef arma::uword   uint_val;
    typedef arma::sword   int_val;
    typedef float         val;
    typedef bool          logic;
}

namespace cnst {
    const type::val    eps   = arma::datum::eps;
    const type::bitval bit_0 = 0b0;
    const type::bitval bit_1 = 0b1;
    const type::logic  no    = false;
    const type::logic  yes   = true;
}

#endif
