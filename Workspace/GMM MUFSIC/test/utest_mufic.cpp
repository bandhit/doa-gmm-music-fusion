#ifdef USING_UNITTEST

#include <gtest/gtest.h>
#include <mufsic.h>
#include <port_adio.h>
#include <misc.h>

TEST(TEST_MATH_WIN_FCN, BLACK_MAN_WIN_FCN_1) {
    type::vec  true_win = {0,
                           0.090453424354128,
                           0.459182957545964,
                           0.920363618099908,
                           0.920363618099908,
                           0.459182957545964,
                           0.090453424354128,
                           0};
    type::vec  win;
    type::uint n = true_win.n_rows;
    blak_man_win_fcn(win, n);
    EXPECT_TRUE(arma::approx_equal(win, true_win, "absdiff", cnst::mat::eps));
}

TEST(TEST_MATH_WIN_FCN, BLACK_MAN_WIN_FCN_2) {
    type::vec  true_win = {0,
                           0.130000000000000,
                           0.630000000000000,
                           1.000000000000000,
                           0.630000000000000,
                           0.130000000000000,
                           0};
    type::vec  win;
    type::uint n = true_win.n_rows;
    blak_man_win_fcn(win, n);
    EXPECT_TRUE(arma::approx_equal(win, true_win, "absdiff", cnst::mat::eps));
}

TEST(TEST_MATH_WIN_FCN, HAMMING_WIN_FCN_1) {
    type::vec  true_win = {0.080000000000000,
                           0.253194691144983,
                           0.642359629619905,
                           0.954445679235113,
                           0.954445679235113,
                           0.642359629619905,
                           0.253194691144983,
                           0.080000000000000};
    type::vec  win;
    type::uint n = true_win.n_rows;
    hamg_man_win_fcn(win, n);
    EXPECT_TRUE(arma::approx_equal(win, true_win, "absdiff", cnst::mat::eps));
}

TEST(TEST_MATH_WIN_FCN, HAMMING_WIN_FCN_2) {
    type::vec  true_win = {0.080000000000000,
                           0.310000000000000,
                           0.770000000000000,
                           1,
                           0.770000000000000,
                           0.310000000000000,
                           0.080000000000000};
    type::vec  win;
    type::uint n = true_win.n_rows;
    hamg_man_win_fcn(win, n);
    EXPECT_TRUE(arma::approx_equal(win, true_win, "absdiff", cnst::mat::eps));
}

TEST(TEST_MATH_FFT, STFT_1) {
    type::cx_cub out_stft;
    type::vec    out_frq;
    type::vec    out_tim;
    type::vec    sim_tim;
    type::mat    sim_in_istft;
    type::vec    win;
    type::val    sam_frq;
    type::val    sam_tim;
    type::val    tim;
    type::uint   n_fft;
    type::uint   hop_size;
    type::uint   win_size;

    sam_frq  = 192E+3;
    sam_tim  = 1 / sam_frq;
    tim      = 1;
    n_fft    = std::pow(2, (13 + 1)) - 1;
    win_size = std::floor(100E-3 / (1 / sam_frq));
    hop_size = std::ceil(win_size / std::pow(2, 4));
    blak_man_win_fcn(win, win_size);

    sim_tim = arma::regspace<type::vec>(0, sam_tim, tim);
    sim_in_istft.resize(sim_tim.n_elem, 2);
    sim_in_istft.col(0) = arma::sin(2 * cnst::mat::pi * 5E+3 * sim_tim);
    sim_in_istft.col(1) = arma::sin(2 * cnst::mat::pi * 2E+3 * sim_tim);
    //sim_in_istft.col(2) = arma::sin(2 * cnst::mat::pi * 1E+3 * sim_tim);

    cap_tim cap_tim_obj;
    cap_tim_obj.tic();
    stft_cub(out_stft, out_frq, out_tim, sim_in_istft, win, sam_frq, n_fft, hop_size);
    cap_tim_obj.toc();
    std::cout << GLOG() << "Around " << cap_tim_obj.get_tim() << " usec." << std::endl;
    FAIL() << "Not done yet!";
}

#endif
