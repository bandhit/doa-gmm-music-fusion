#ifndef _MUFSIC_H_
#define _MUFSIC_H_

#include <type_def.hpp>

namespace fftw3 {
    void fft_1d (type::cx_vec& out, const type::cx_vec& in, const type::uint n);
}

void blak_man_win_fcn (type::vec& out, const type::uint n_pnt);

void hamg_man_win_fcn (type::vec& out, const type::uint n_pnt);

void stft_cub (type::cx_cub&    out_stft,
               type::vec&       out_frq,
               type::vec&       out_tim,
               const type::mat& in_istft,
               const type::vec& win,
               const type::val  sam_frq,
               const type::uint n_fft,
               const type::uint hop_size);

/*void istft_cub (type::mat&    out_istft,
                type::cx_cub& in_stft,
                type::vec&    win,
                type::val     sam_frq,
                type::uint    nfft,
                type::uint    hop_size);*/

#endif
