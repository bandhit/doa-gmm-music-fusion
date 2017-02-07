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
    using bit_cub  = arma::Cube<bit>;
    using sint_cub = arma::Cube<sint>;
    using uint_cub = arma::Cube<uint>;
    using cub      = arma::Cube<val>;
    using cx_cub   = arma::Cube<cx_val>;
    template <typename t>
    using abs_vec  = arma::Col<t>;
    template <typename t>
    using abs_mat  = arma::Mat<t>;
    template <typename t>
    using abs_cub  = arma::Cube<t>;
}

namespace cnst {
    const type::bit b_0 = 0b0;
    const type::bit b_1 = 0b1;
    using mat = arma::Datum<type::val>;
}

namespace optr {
    static type::val approx_equal (type::val a, type::val b, type::val eps) {
        return (std::abs(a - b) < eps);
    }
}

#endif
