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

void fftw3::ifft (type::cx_vec& out, type::cx_vec& in, const type::uint n) {
    #if not defined(ALL_USING_DOUBLE_TYPE)
    fftwf_plan p = fftwf_plan_dft_1d(n,
                                    (type::val(*)[2])in.memptr(),
                                    (type::val(*)[2])out.memptr(),
                                     FFTW_BACKWARD,
                                     FFTW3_FLAG);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    #else
    fftw_plan p = fftw_plan_dft_1d(n,
                                   (type::val(*)[2])&in(0),
                                   (type::val(*)[2])&out(0),
                                   FFTW_BACKWARD,
                                   FFTW3_FLAG);
    fftw_execute(p);
    fftw_destroy_plan(p);
    #endif
}

#if defined(FFTW3_USING_THREADS)
void fftw3::init (void) {
    #if not defined(ALL_USING_DOUBLE_TYPE)
    fftwf_init_threads();
    fftwf_plan_with_nthreads(FFTW3_NUMBER_OF_THREADS);
    #else
    fftw_init_threads();
    fftw_plan_with_nthreads(FFTW3_NUMBER_OF_THREADS);
    #endif
}

void fftw3::de_init (void) {
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
    out.set_size(n_pnt);
    out.rows(0, n_haft - 1) =  0.42
                            - (0.50 * arma::cos((2 * cnst::mat::pi * n) / (n_pnt - 1)))
                            + (0.08 * arma::cos((4 * cnst::mat::pi * n) / (n_pnt - 1)));
    out.rows(n_haft,
             n_haft + idx)  = arma::flipud(out.rows(0, idx));
}

void hamg_man_win_fcn (type::vec& out, const type::uint n_pnt) {
    type::vec n = arma::regspace<type::vec>(0, 1, n_pnt - 1);
    out =  0.54
        - (0.46 * arma::cos(2 * cnst::mat::pi * (n / (n_pnt - 1))));
}

void stft_slice (type::cx_cub&    out,
                 type::vec&       out_frq,
                 type::vec&       out_tim,
                 const type::mat& in,
                 const type::vec& win,
                 const type::val  sam_frq,
                 const type::uint n_fft,
                 const type::uint hop_size) {
    type::uint     n_sam     = in.n_cols;
    type::uint     n_tim     = in.n_rows;
    type::uint     win_size  = win.n_elem;
    type::uint     n_frq_fft = std::ceil(((type::val) n_fft + 1) / 2);
    type::uint     n_tim_fft = 1 + std::floor(((type::val) n_tim - win_size) / hop_size);

    type::uint_vec head_hop(n_tim_fft);
    type::uint_vec tail_hop(n_tim_fft);
    for (type::uint i = 0, i_hop = 0; i < n_tim_fft; i++, i_hop = i_hop + hop_size) {
        head_hop(i) = i_hop;
        tail_hop(i) = i_hop + win_size - 1;
    }

    type::uint h_frq_fft = 0;
    type::uint t_frq_fft = n_frq_fft - 1;

    #if defined(MATH_USING_FFTW3)
    #if defined(FFTW3_USING_THREADS)
    fftw3::init();
    #endif
    type::cx_mat cx_in = arma::conv_to<type::cx_mat>::from(in);
    type::cx_vec tmp_cx_in(n_fft);
    #endif

    out.set_size(n_frq_fft, n_tim_fft, n_sam);
    type::cx_vec tmp_fft(n_fft);
    for (type::uint i_sam = 0, i_tim_fft; i_sam < n_sam; i_sam++) {
        for (i_tim_fft = 0; i_tim_fft < n_tim_fft; i_tim_fft++) {
    #if not defined(MATH_USING_FFTW3)
            tmp_fft = arma::fft(win % in.col(i_sam)
                                        .subvec(head_hop(i_tim_fft),
                                                tail_hop(i_tim_fft)),
                                n_fft);
    #else
            tmp_cx_in = win % cx_in.col(i_sam)
                                   .subvec(head_hop(i_tim_fft),
                                           tail_hop(i_tim_fft));
            fftw3::fft(tmp_fft, tmp_cx_in, n_fft);
    #endif
            out.slice(i_sam)
               .col(i_tim_fft) = tmp_fft.subvec(h_frq_fft, t_frq_fft);
        }
    }

    #if defined(MATH_USING_FFTW3) and defined(FFTW3_USING_THREADS)
    fftw3::de_init();
    #endif

    out_frq  = arma::regspace<type::vec>(0, 1, n_frq_fft - 1);
    out_frq *= (sam_frq / n_fft);
    out_tim  = arma::regspace<type::vec>(win_size / 2, hop_size,
                                         (win_size / 2) + ((n_tim_fft - 1) * hop_size));
    out_tim /= sam_frq;
}

void istft_slice (type::mat&          out,
                  type::vec&          out_tim,
                  const type::cx_cub& in,
                  const type::vec&    win,
                  const type::val     sam_frq,
                  const type::uint    n_fft,
                  const type::uint    hop_size) {
    type::uint n_frq_fft  = in.n_rows;
    type::uint n_tim_fft  = in.n_cols;
    type::uint n_tim_ifft = n_fft + (n_tim_fft - 1) * hop_size;
    type::uint n_sam      = in.n_slices;

    type::uint_vec head_hop(n_tim_fft);
    type::uint_vec tail_hop(n_tim_fft);
    for (type::uint i = 0, i_hop = 0; i < n_tim_fft; i++, i_hop = i_hop + hop_size) {
        head_hop(i) = i_hop;
        tail_hop(i) = i_hop + n_fft - 1;
    }

    type::uint h_l_tmp_in = 0;
    type::uint t_l_tmp_in = n_frq_fft - 1;
    type::uint h_c_tmp_in = 1;
    type::uint t_c_tmp_in = n_frq_fft - 1;
    type::uint h_r_tmp_in = n_frq_fft;
    type::uint t_r_tmp_in = (2 * n_frq_fft) - 2;
    if ((n_fft % 2) == 0) {
        t_c_tmp_in--;
        t_r_tmp_in--;
    }
    type::uint n_tmp_in = t_r_tmp_in + 1;

    #if defined(MATH_USING_FFTW3) and defined(FFTW3_USING_THREADS)
    fftw3::init();
    #endif

    out.set_size(n_tim_ifft, n_sam);
    out.zeros();
    type::cx_vec tmp_in(n_tmp_in);
    type::cx_vec tmp_ifft(n_tmp_in);
    for (type::uint i_sam = 0, i_tim_fft; i_sam < n_sam; i_sam++) {
        for (i_tim_fft = 0; i_tim_fft < n_tim_fft; i_tim_fft++) {
            tmp_in.rows(h_l_tmp_in,
                        t_l_tmp_in) = in.slice(i_sam).col(i_tim_fft);
            tmp_in.rows(h_r_tmp_in,
                        t_r_tmp_in) = arma::conj(arma::flipud(tmp_in.rows(h_c_tmp_in,
                                                                          t_c_tmp_in)));
    #if not defined(MATH_USING_FFTW3)
            tmp_ifft = arma::ifft(tmp_in, n_tmp_in);
    #else
            fftw3::ifft(tmp_ifft, tmp_in, n_tmp_in);
    #endif
            out.col(i_sam)
               .rows(head_hop(i_tim_fft),
                     tail_hop(i_tim_fft)) += win % arma::real(tmp_ifft);
        }
    }

    #if defined(MATH_USING_FFTW3) and defined(FFTW3_USING_THREADS)
    fftw3::de_init();
    #endif

    #if defined(MATH_USING_FFTW3)
    out /= n_tmp_in;
    #endif
    out *= hop_size;
    out /= arma::as_scalar(arma::sum(arma::pow(win, 2)));

    out_tim  = arma::regspace<type::vec>(0, 1, n_tim_ifft - 1);
    out_tim /= sam_frq;
}

void sam_cov_zero_mean (type::mat& out, const type::mat& in) {
    type::vec avg_in = arma::mean(in, 0).t();
    out              = arma::cov(in, 1) + (avg_in * avg_in.t());
}

void sam_cov_zero_mean (type::cx_mat& out, const type::cx_mat& in) {
    type::cx_vec avg_in = arma::mean(in, 0).st();
    out                 = arma::cov(in, 1) + arma::conj(avg_in * avg_in.t());
}

void root (type::cx_vec& out, const type::vec& in) {
    type::mat a  = arma::diagmat(arma::ones<type::vec>(in.n_elem - 2), -1);
    a.row(0)     = in.subvec(1, in.n_elem - 1).t();
    a.row(0)    /= - in(0);
    arma::eig_gen(out, a);
}

void root (type::cx_vec& out, const type::cx_vec& in) {
    type::cx_vec d = - in.subvec(1, in.n_elem - 1) / in(0);
    if (arma::approx_equal(arma::imag(d), arma::zeros<type::vec>(in.n_elem - 1), "absdiff", cnst::tol)) {
        type::mat    a = arma::diagmat(arma::ones<type::vec>(in.n_elem - 2), -1);
        a.row(0)       = arma::real(d).t();
        arma::eig_gen(out, a);
    }
    else {
        type::cx_mat a = arma::diagmat(arma::ones<type::cx_vec>(in.n_elem - 2), -1);
        a.row(0)       = d.st();
        arma::eig_gen(out, a);
    }
}
