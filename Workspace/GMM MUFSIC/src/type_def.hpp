#ifndef _TYPE_DEF_H_
#define _TYPE_DEF_H_

#include <armadillo>

#define ALL_USING_DOUBLE_TYPE
//#define PORT_ADIO_USING_STD_MEMCPY
#define MATH_USING_FFTW3
#define FFTW3_FLAG FFTW_ESTIMATE // FFTW_MEASURE or FFTW_ESTIMATE
#define FFTW3_USING_THREADS
#define FFTW3_NUMBER_OF_THREADS 2

namespace type {
    using sint8        = int8_t;
    using sint16       = int16_t;
    using sint32       = int32_t;
    using sint64       = int64_t;
    using uint8        = uint8_t;
    using uint16       = uint16_t;
    using uint32       = uint32_t;
    using uint64       = uint64_t;
    using sint         = sint32;
    using uint         = uint32;
    #if not defined(ALL_USING_DOUBLE_TYPE)
    using val          = float;
    #else
    using val          = double;
    #endif
    using cx_val       = std::complex<val>;
    using sint_vec     = arma::Col<sint>;
    using uint_vec     = arma::Col<uint>;
    using vec          = arma::Col<val>;
    using cx_vec       = arma::Col<cx_val>;
    using sint_row_vec = arma::Row<sint>;
    using uint_row_vec = arma::Row<uint>;
    using row_vec      = arma::Row<val>;
    using cx_row_vec   = arma::Row<cx_val>;
    using sint_mat     = arma::Mat<sint>;
    using uint_mat     = arma::Mat<uint>;
    using mat          = arma::Mat<val>;
    using cx_mat       = arma::Mat<cx_val>;
    using sint_cub     = arma::Cube<sint>;
    using uint_cub     = arma::Cube<uint>;
    using cub          = arma::Cube<val>;
    using cx_cub       = arma::Cube<cx_val>;
    template <typename t>
    using abs_vec      = arma::Col<t>;
    template <typename t>
    using abs_row_vec  = arma::Row<t>;
    template <typename t>
    using abs_mat      = arma::Mat<t>;
    template <typename t>
    using abs_cub      = arma::Cube<t>;
    using idx_vec      = arma::uvec;
}

namespace cnst {
    using mat = arma::Datum<type::val>;
    const type::cx_val i = type::cx_val(0, 1);
    // reduce precision
    #if not defined(ALL_USING_DOUBLE_TYPE)
    const type::val tol = mat::eps * (type::val) 1E+2;
    #else
    const type::val tol = mat::eps * (type::val) 1E+2;
    #endif
}

#endif
