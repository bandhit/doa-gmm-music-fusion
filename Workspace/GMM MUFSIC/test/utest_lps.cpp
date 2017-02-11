#ifdef USING_UNITTEST

#include <gtest/gtest.h>
#include <lps.h>
#include <dsp.h>
#include <misc.h>

TEST(TEST_LPS, MUSIC_1) {
    type::mat in;
    cap_tim cap_tim_1_obj;
    cap_tim_1_obj.tic();
    in.load("test/ex_lps_1.csv", arma::csv_ascii);
    cap_tim_1_obj.toc();
    type::val    sam_frq       = 44.1E+3;
    type::uint   n_fft         = std::pow(2, (13 + 1));
    type::uint   win_size_stft = std::floor(100E-3 / (1 / sam_frq));
    type::uint   hop_size      = std::ceil(win_size_stft / std::pow(2, 4));
    type::vec    win_stft;
    blak_man_win_fcn(win_stft, win_size_stft);
    type::cx_cub out_stft;
    type::vec    out_frq;
    type::vec    out_tim_stft;
    cap_tim cap_tim_2_obj;
    cap_tim_2_obj.tic();
    stft_slice(out_stft, out_frq, out_tim_stft, in, win_stft, sam_frq, n_fft, hop_size);
    out_stft /= arma::sum(win_stft);
    cap_tim_2_obj.toc();

    type::uint t_tim = out_stft.n_cols - 1;

    type::val fun_frq     = 3.4E+3;;
    type::val c           = 340;
    type::val d           = (c / fun_frq) / 2;
    type::val min_eig     = 1E-3;
    type::val roi_min_frq = fun_frq / 4;
    type::val roi_max_frq = fun_frq;
    type::vec phi_rad = arma::regspace<type::vec>(-cnst::mat::pi / 2,
                                                  1E-2,
                                                  cnst::mat::pi / 2);

    cap_tim cap_tim_3_obj;
    cap_tim_3_obj.tic();
    type::idx_vec min_idx = arma::find(out_frq >= roi_min_frq, 1, "first");
    type::idx_vec max_idx = arma::find(out_frq <= roi_max_frq, 1, "last");
    cap_tim_3_obj.toc();
    type::uint    h_idx   = !min_idx.is_empty() ? min_idx(0) : 0;
    type::uint    t_idx   = !max_idx.is_empty() ? max_idx(0) : t_tim;
    type::uint    nfrq    = t_idx - h_idx + 1;

    type::uint   h_frq;
    type::cx_mat tmp_out_stft;
    type::vec    tmp_out;
    type::mat    out(phi_rad.n_elem, nfrq);
    cap_tim cap_tim_4_obj;
    cap_tim_4_obj.tic();
    for (type::uint i = 0; i < nfrq; i++) {
        h_frq        = h_idx + i;
        tmp_out_stft = arma::reshape(out_stft.tube(h_frq, 0, h_frq, t_tim),
                                     out_stft.n_cols, out_stft.n_slices, 1);
        music_min_eig(tmp_out, tmp_out_stft, phi_rad, d, c, fun_frq, out_frq(h_frq), min_eig);
        out.col(i)   = tmp_out;
    }
    out = out.t();
    out = arma::fliplr(out);
    out = arma::flipud(out);
    out = 10 * arma::log10(out + (type::val) 1E-2);
    out.elem(arma::find(out <= 0)).zeros();
    cap_tim_4_obj.toc();

    std::cout << GLOG() << "Open  around " << cap_tim_1_obj.get_tim() / 1000.0 << " msec." << std::endl;
    std::cout << GLOG() << "STFT  around " << cap_tim_2_obj.get_tim() / 1000.0 << " msec." << std::endl;
    std::cout << GLOG() << "Cut   around " << cap_tim_3_obj.get_tim() / 1000.0 << " msec." << std::endl;
    std::cout << GLOG() << "MUSIC around " << cap_tim_4_obj.get_tim() / 1000.0 << " msec." << std::endl;

    phi_rad.save("out/ex_music_1.csv", arma::csv_ascii);
    out_frq = out_frq.subvec(h_idx, t_idx);
    out_frq.save("out/ex_music_2.csv", arma::csv_ascii);
    out.save("out/ex_music_3.csv", arma::csv_ascii);

    SUCCEED();
}

#endif
