#include <mufsic.h>

void blak_man_win_fcn (type::vec& out,
                       type::uint n_pnt) {
    type::uint n_haft;
    type::uint idx;
    if ((n_pnt % 2) == 0) {
        n_haft = n_pnt / 2;
        idx    = n_haft - 1;
    }
    else {
        n_haft = (n_pnt + 1) / 2;
        idx    = n_haft - 2;
    }
    type::vec n = arma::regspace<type::vec>(0, 1, n_haft - 1);
    out.resize(n_pnt);
    out.rows(0, n_haft - 1) =  0.42
                            - (0.50 * arma::cos((2 * cnst::mat::pi * n) / (n_pnt - 1)))
                            + (0.08 * arma::cos((4 * cnst::mat::pi * n) / (n_pnt - 1)));
    out.rows(n_haft,
             n_haft + idx)  = flipud(out.rows(0, idx));
}

void hamg_man_win_fcn (type::vec& out,
                       type::uint n_pnt) {
    type::vec n = arma::regspace<type::vec>(0, 1, n_pnt - 1);
    out =  0.54
        - (0.46 * arma::cos(2 * cnst::mat::pi * (n / (n_pnt - 1))));
}

void stft_cub (type::cx_cub& out_stft,
               type::vec&    out_frq,
               type::vec&    out_tim,
               type::mat&    in_istft,
               type::vec&    win,
               type::val     sam_frq,
               type::uint    n_fft,
               type::uint    hop_size) {
    type::uint n_tim     = in_istft.n_rows;
    type::uint n_sam     = in_istft.n_cols;
    type::uint win_size  = win.n_elem;
    type::uint n_frq_fft = std::ceil((n_fft + 1) / 2);
    type::uint n_tim_fft = 1 + std::floor((n_tim - win_size) / hop_size);
    out_frq = (sam_frq / n_fft)
            * arma::regspace<type::vec>(0, 1, n_frq_fft - 1);
    out_tim = arma::regspace<type::vec>(win_size / 2,
                                        hop_size,
                                        (win_size / 2) + ((n_tim_fft - 1) * hop_size))
            / sam_frq;

    type::uint_vec head_hop(n_tim_fft);
    type::uint_vec tail_hop(n_tim_fft);
    for (type::uint idx = 0, idx_hop = 0; idx < n_tim_fft; idx++,
                                                           idx_hop = idx_hop + hop_size) {
        head_hop(idx) = idx_hop;
        tail_hop(idx) = idx_hop + win_size - 1;
    }

    type::cx_vec tmp_fft;
    type::uint tail_frq_fft = n_frq_fft - 1;
    out_stft.resize(n_frq_fft, n_tim_fft, n_sam);
    for (type::uint idx_sam = 0, idx_tim_fft; idx_sam < n_sam; idx_sam++) {
        for (idx_tim_fft = 0; idx_tim_fft < n_tim_fft; idx_tim_fft++) {
            tmp_fft = arma::fft(in_istft.col(idx_sam)
                                        .subvec(head_hop(idx_tim_fft),
                                                tail_hop(idx_tim_fft))
                               % win,
                                n_fft);
            out_stft.slice(idx_sam)
                    .col(idx_tim_fft) = tmp_fft.rows(0, tail_frq_fft);
        }
    }
}
