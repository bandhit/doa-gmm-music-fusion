#ifndef _MUFSIC_H_
#define _MUFSIC_H_

#include <type_def.hpp>

void blak_man_win_fcn (type::vec& out,
                       type::uint n_pnt);

void hamg_man_win_fcn (type::vec& out,
                       type::uint n_pnt);

void stft_cub (type::cx_cub& out_stft,
               type::vec&    out_frq,
               type::vec&    out_tim,
               type::mat&    in_istft,
               type::vec&    win,
               type::val     sam_frq,
               type::uint    n_fft,
               type::uint    hop_size);

/*void istft_cub (type::mat&    out_istft,
                type::cx_cub& in_stft,
                type::vec&    win,
                type::val     sam_frq,
                type::uint    nfft,
                type::uint    hop_size);*/

#endif
