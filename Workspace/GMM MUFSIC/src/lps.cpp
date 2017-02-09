#include <lps.h>

void music_min_eig (type::vec&       abs_out,
                    const type::mat& in,
                    const type::vec& phi_rad,
                    const type::val  d,
                    const type::val  c,
                    const type::val  frq,
                    const type::val  min_eig) {
    type::uint n_tim = in.n_rows;
    type::uint n_sam = in.n_cols;
    type::val  lamb  = c / frq;
}
