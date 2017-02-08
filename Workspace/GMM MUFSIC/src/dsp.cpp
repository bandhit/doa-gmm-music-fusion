#include <dsp.h>

#if defined(MATH_USING_FFTW3)
#include <fftw3.h>

void fftw3::fft (type::cx_vec& out, type::cx_vec& in, const type::uint n) {
    #if not defined(ALL_USING_DOUBLE_TYPE)
    fftwf_plan p = fftwf_plan_dft_1d(n,
                                    (type::val(*)[2])in.memptr(),
                                    (type::val(*)[2])out.memptr(),
                                     FFTW_FORWARD,
                                     FFTW3_FLAG);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    #else
    fftw_plan p = fftw_plan_dft_1d(n,
                                   (type::val(*)[2])&in(0),
                                   (type::val(*)[2])&out(0),
                                   FFTW_FORWARD,
                                   FFTW3_FLAG);
    fftw_execute(p);
    fftw_destroy_plan(p);
    #endif
}

#if defined(FFTW3_USING_THREADS)
void fftw3::init_fft (void) {
    #if not defined(ALL_USING_DOUBLE_TYPE)
    fftwf_init_threads();
    fftwf_plan_with_nthreads(FFTW3_NUMBER_OF_THREADS);
    #else
    fftw_init_threads();
    fftw_plan_with_nthreads(FFTW3_NUMBER_OF_THREADS);
    #endif
}

void fftw3::de_init_fft (void) {
    #if not defined(ALL_USING_DOUBLE_TYPE)
    fftwf_cleanup_threads();
    #else
    fftw_cleanup_threads();
    #endif
}

#endif
#endif

void blak_man_win_fcn (type::vec& out, const type::uint n_pnt) {
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

void hamg_man_win_fcn (type::vec& out, const type::uint n_pnt) {
    type::vec n = arma::regspace<type::vec>(0, 1, n_pnt - 1);
    out =  0.54
        - (0.46 * arma::cos(2 * cnst::mat::pi * (n / (n_pnt - 1))));
}

void stft_slice (type::cx_cub&    out_stft,
                 type::vec&       out_frq,
                 type::vec&       out_tim,
                 const type::mat& in,
                 const type::vec& win,
                 const type::val  sam_frq,
                 const type::uint n_fft,
                 const type::uint hop_size) {
    #if defined(MATH_USING_FFTW3)
    #if defined(FFTW3_USING_THREADS)
    fftw3::init_fft();
    #endif
    type::cx_mat cx_in = arma::conv_to<type::cx_mat>::from(in);
    type::cx_vec tmp_cx_in(n_fft);
    #endif
    type::uint     n_sam     = in.n_cols;
    type::uint     n_tim     = in.n_rows;
    type::uint     win_size  = win.n_elem;
    type::uint     n_frq_fft = std::ceil((n_fft + 1) / 2);
    type::uint     n_tim_fft = 1 + std::floor((n_tim - win_size) / hop_size);
    type::uint_vec head_hop(n_tim_fft);
    type::uint_vec tail_hop(n_tim_fft);
    for (type::uint idx = 0, idx_hop = 0; idx < n_tim_fft; idx++, idx_hop = idx_hop + hop_size) {
        head_hop(idx) = idx_hop;
        tail_hop(idx) = idx_hop + win_size - 1;
    }

    out_stft.resize(n_frq_fft, n_tim_fft, n_sam);
    type::cx_vec tmp_fft(n_fft);
    type::uint   tail_frq_fft = n_frq_fft - 1;
    for (type::uint idx_sam = 0, idx_tim_fft; idx_sam < n_sam; idx_sam++) {
        for (idx_tim_fft = 0; idx_tim_fft < n_tim_fft; idx_tim_fft++) {
    #if not defined(MATH_USING_FFTW3)
            tmp_fft = arma::fft(win % in.col(idx_sam)
                                        .subvec(head_hop(idx_tim_fft), tail_hop(idx_tim_fft)),
                                n_fft);
    #else
            tmp_cx_in = win % cx_in.col(idx_sam)
                                   .subvec(head_hop(idx_tim_fft), tail_hop(idx_tim_fft));
            fftw3::fft(tmp_fft, tmp_cx_in, n_fft);
            out_stft.slice(idx_sam)
                    .col(idx_tim_fft) = tmp_fft.subvec(0, tail_frq_fft);
        }
    #endif
    }

    out_frq = (sam_frq / n_fft) * arma::regspace<type::vec>(0, 1, n_frq_fft - 1);
    out_tim = arma::regspace<type::vec>(win_size / 2, hop_size,
                                        (win_size / 2) + ((n_tim_fft - 1) * hop_size)) / sam_frq;
    #if defined(MATH_USING_FFTW3) and defined(FFTW3_USING_THREADS)
    fftw3::de_init_fft();
    #endif
}
