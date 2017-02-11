#ifndef _LPS_H_
#define _LPS_H_

#include <type_def.hpp>

void frwd_w_fcn (type::cx_mat&        out,
                 const type::row_vec& phi_rad,
                 const type::uint     n_sam,
                 const type::val      d,
                 const type::val      fun_frq,
                 const type::val      frq,
                 const type::val      lamb);

void frwd_w_fcn (type::cx_vec&    out,
                 const type::val  phi_rad,
                 const type::uint n_sam,
                 const type::val  d,
                 const type::val  fun_frq,
                 const type::val  frq,
                 const type::val  lamb);

void music_min_eig (type::vec&          abs_out,
                    const type::cx_mat& in,
                    const type::vec&    phi_rad,
                    const type::val     d,
                    const type::val     c,
                    const type::val     fun_frq,
                    const type::val     frq,
                    const type::val     min_eig);

#endif
