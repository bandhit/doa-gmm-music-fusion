#ifndef _TYPE_DEF_H_
#define _TYPE_DEF_H_

#include <armadillo>

//#define PORT_ADIO_USING_STD_MEMCPY

namespace type {
    using sint8    = int8_t;
    using sint16   = int16_t;
    using sint32   = int32_t;
    using sint64   = int64_t;
    using uint8    = uint8_t;
    using uint16   = uint16_t;
    using uint32   = uint32_t;
    using uint64   = uint64_t;
    using bit      = uint8;
    using sint     = sint32;
    using uint     = uint32;
    using val      = float;
    using cx_val   = std::complex<val>;
    using bit_vec  = arma::Col<bit>;
    using sint_vec = arma::Col<sint>;
    using uint_vec = arma::Col<uint>;
    using vec      = arma::Col<val>;
    using cx_vec   = arma::Col<cx_val>;
    using bit_mat  = arma::Mat<bit>;
    using sint_mat = arma::Mat<sint>;
    using uint_mat = arma::Mat<uint>;
    using mat      = arma::Mat<val>;
    using cx_mat   = arma::Mat<cx_val>;
    template <typename t>
    using abs_vec  = arma::Col<t>;
    template <typename t>
    using abs_mat  = arma::Mat<t>;
}

namespace cnst {
    const type::bit b_0 = 0b0;
    const type::bit b_1 = 0b1;
    const type::val eps = 1e-8;
}

#endif
