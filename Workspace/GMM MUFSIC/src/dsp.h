#ifndef _MUFSIC_H_
#define _MUFSIC_H_

#include <type_def.hpp>

#if defined(MATH_USING_FFTW3)
namespace fftw3 {
    void fft (type::cx_vec& out, type::cx_vec& in, const type::uint n);
    void ifft (type::cx_vec& out, type::cx_vec& in, const type::uint n);
    #if defined(FFTW3_USING_THREADS)
    void init (void);
    void de_init (void);
    #endif
}
#endif

void blak_man_win_fcn (type::vec& out, const type::uint n_pnt);

void hamg_man_win_fcn (type::vec& out, const type::uint n_pnt);

void stft_slice (type::cx_cub&    out,
                 type::vec&       out_frq,
                 type::vec&       out_tim,
                 const type::mat& in,
                 const type::vec& win,
                 const type::val  sam_frq,
                 const type::uint n_fft,
                 const type::uint hop_size);

void istft_slice (type::mat&          out,
                  type::vec&          out_tim,
                  const type::cx_cub& in,
                  const type::vec&    win,
                  const type::val     sam_frq,
                  const type::uint    n_fft,
                  const type::uint    hop_size);

#endif
