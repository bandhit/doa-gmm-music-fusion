#ifndef _LPS_H_
#define _LPS_H_

#include <type_def.hpp>

void music_min_eig (type::vec&       abs_out,
                    const type::mat& in,
                    const type::vec& phi_rad,
                    const type::val  d,
                    const type::val  c,
                    const type::val  frq,
                    const type::val  min_eig);

#endif
