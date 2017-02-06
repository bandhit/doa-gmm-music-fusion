#ifndef _MUFSIC_H_
#define _MUFSIC_H_

#include <type_def.hpp>

void blak_man_win_fcn (type::mat& out,
                       type::uint n);

void hamg_man_win_fcn (type::mat& out,
                       type::uint n);

void stft_cub (type::cx_cub& out,
               type::mat&    in,
               type::val     sam_frq,
               void*         win_fcn(type::mat& win,
                                     type::uint n_win),
               type::uint    win_size,
               type::uint    nfft,
               type::uint    hop_size);

void istft_cub (void);

#endif
