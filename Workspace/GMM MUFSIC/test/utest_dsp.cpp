
#ifdef USING_UNITTEST

#include <gtest/gtest.h>
#include <dsp.h>
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
    EXPECT_TRUE(arma::approx_equal(win, true_win, "absdiff", cnst::tol));
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
    EXPECT_TRUE(arma::approx_equal(win, true_win, "absdiff", cnst::tol));
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
    EXPECT_TRUE(arma::approx_equal(win, true_win, "absdiff", cnst::tol));
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
    EXPECT_TRUE(arma::approx_equal(win, true_win, "absdiff", cnst::tol));
}

TEST(TEST_MATH_FFT, STFT_1) {
    type::cx_cub out_stft;
    type::vec    out_frq;
    type::vec    out_tim;
    type::vec    sim_tim;
    type::mat    sim_in;
    type::vec    win;
    type::val    sam_frq;
    type::val    sam_tim;
    type::val    tim;
    type::uint   n_fft;
    type::uint   hop_size;
    type::uint   win_size;

    for (type::uint i = 1; i <= 13; i++) {
        sam_frq  = 16E+3;
        sam_tim  = 1 / sam_frq;
        tim      = 1;
        n_fft    = std::pow(2, (i + 1)) - 1;
        win_size = std::floor(100E-3 / (1 / sam_frq));
        hop_size = std::ceil(win_size / std::pow(2, 4));
        blak_man_win_fcn(win, win_size);

        sim_tim = arma::regspace<type::vec>(0, sam_tim, tim);
        sim_in.resize(sim_tim.n_elem, 8);
        sim_in.col(0) = arma::sin(2 * cnst::mat::pi * 2E+2 * sim_tim);
        sim_in.col(1) = arma::sin(2 * cnst::mat::pi * 4E+2 * sim_tim);
        sim_in.col(2) = arma::sin(2 * cnst::mat::pi * 6E+2 * sim_tim);
        sim_in.col(3) = arma::sin(2 * cnst::mat::pi * 8E+2 * sim_tim);
        sim_in.col(4) = arma::sin(2 * cnst::mat::pi * 2E+3 * sim_tim);
        sim_in.col(5) = arma::sin(2 * cnst::mat::pi * 4E+3 * sim_tim);
        sim_in.col(6) = arma::sin(2 * cnst::mat::pi * 6E+3 * sim_tim);
        sim_in.col(7) = arma::sin(2 * cnst::mat::pi * 8E+3 * sim_tim);

        cap_tim cap_tim_obj;
        cap_tim_obj.tic();
        stft_slice(out_stft, out_frq, out_tim, sim_in, win, sam_frq, n_fft, hop_size);
        cap_tim_obj.toc();
        std::cout << GLOG() << "  16K around " << cap_tim_obj.get_tim() / 1000.0 << " msec,"
                            << " with next pow of 2 #fft "   << i << "." <<  std::endl;
        sim_tim.reset();
        sim_in.reset();
        out_stft.reset();
        out_frq.reset();
        out_tim.reset();
    }
    SUCCEED();
}

TEST(TEST_MATH_FFT, STFT_2) {
    type::cx_cub out_stft;
    type::vec    out_frq;
    type::vec    out_tim;
    type::vec    sim_tim;
    type::mat    sim_in;
    type::vec    win;
    type::val    sam_frq;
    type::val    sam_tim;
    type::val    tim;
    type::uint   n_fft;
    type::uint   hop_size;
    type::uint   win_size;

    for (type::uint i = 1; i <= 13; i++) {
        sam_frq  = 44.1E+3;
        sam_tim  = 1 / sam_frq;
        tim      = 1;
        n_fft    = std::pow(2, (i + 1)) - 1;
        win_size = std::floor(100E-3 / (1 / sam_frq));
        hop_size = std::ceil(win_size / std::pow(2, 4));
        blak_man_win_fcn(win, win_size);

        sim_tim = arma::regspace<type::vec>(0, sam_tim, tim);
        sim_in.resize(sim_tim.n_elem, 8);
        sim_in.col(0) = arma::sin(2 * cnst::mat::pi * 2E+2 * sim_tim);
        sim_in.col(1) = arma::sin(2 * cnst::mat::pi * 4E+2 * sim_tim);
        sim_in.col(2) = arma::sin(2 * cnst::mat::pi * 6E+2 * sim_tim);
        sim_in.col(3) = arma::sin(2 * cnst::mat::pi * 8E+2 * sim_tim);
        sim_in.col(4) = arma::sin(2 * cnst::mat::pi * 2E+3 * sim_tim);
        sim_in.col(5) = arma::sin(2 * cnst::mat::pi * 4E+3 * sim_tim);
        sim_in.col(6) = arma::sin(2 * cnst::mat::pi * 6E+3 * sim_tim);
        sim_in.col(7) = arma::sin(2 * cnst::mat::pi * 8E+3 * sim_tim);

        cap_tim cap_tim_obj;
        cap_tim_obj.tic();
        stft_slice(out_stft, out_frq, out_tim, sim_in, win, sam_frq, n_fft, hop_size);
        cap_tim_obj.toc();
        std::cout << GLOG() << "44.1K around " << cap_tim_obj.get_tim() / 1000.0 << " msec,"
                            << " with next pow of 2 #fft "   << i << "." << std::endl;
        sim_tim.reset();
        sim_in.reset();
        out_stft.reset();
        out_frq.reset();
        out_tim.reset();
    }
    SUCCEED();
}

TEST(TEST_MATH_FFT, STFT_3) {
    type::cx_cub out_stft;
    type::vec    out_frq;
    type::vec    out_tim;
    type::vec    sim_tim;
    type::mat    sim_in;
    type::vec    win;
    type::val    sam_frq;
    type::val    sam_tim;
    type::val    tim;
    type::uint   n_fft;
    type::uint   hop_size;
    type::uint   win_size;

    for (type::uint i = 1; i <= 13; i++) {
        sam_frq  = 192E+3;
        sam_tim  = 1 / sam_frq;
        tim      = 1;
        n_fft    = std::pow(2, (i + 1)) - 1;
        win_size = std::floor(100E-3 / (1 / sam_frq));
        hop_size = std::ceil(win_size / std::pow(2, 4));
        blak_man_win_fcn(win, win_size);

        sim_tim = arma::regspace<type::vec>(0, sam_tim, tim);
        sim_in.resize(sim_tim.n_elem, 8);
        sim_in.col(0) = arma::sin(2 * cnst::mat::pi * 2E+2 * sim_tim);
        sim_in.col(1) = arma::sin(2 * cnst::mat::pi * 4E+2 * sim_tim);
        sim_in.col(2) = arma::sin(2 * cnst::mat::pi * 6E+2 * sim_tim);
        sim_in.col(3) = arma::sin(2 * cnst::mat::pi * 8E+2 * sim_tim);
        sim_in.col(4) = arma::sin(2 * cnst::mat::pi * 2E+3 * sim_tim);
        sim_in.col(5) = arma::sin(2 * cnst::mat::pi * 4E+3 * sim_tim);
        sim_in.col(6) = arma::sin(2 * cnst::mat::pi * 6E+3 * sim_tim);
        sim_in.col(7) = arma::sin(2 * cnst::mat::pi * 8E+3 * sim_tim);

        cap_tim cap_tim_obj;
        cap_tim_obj.tic();
        stft_slice(out_stft, out_frq, out_tim, sim_in, win, sam_frq, n_fft, hop_size);
        cap_tim_obj.toc();
        std::cout << GLOG() << " 192K around " << cap_tim_obj.get_tim() / 1000.0 << " msec,"
                            << " with next pow of 2 #fft "   << i << "." <<  std::endl;
        sim_tim.reset();
        sim_in.reset();
        out_stft.reset();
        out_frq.reset();
        out_tim.reset();
    }
    SUCCEED();
}

TEST(TEST_MATH_FFT, STFT_4) {
    type::cx_cub out_stft;
    type::vec    out_frq;
    type::vec    out_tim;
    type::vec    sim_tim;
    type::mat    sim_in;
    type::vec    win;
    type::val    sam_frq;
    type::val    sam_tim;
    type::val    tim;
    type::uint   n_fft;
    type::uint   hop_size;
    type::uint   win_size;

    sam_frq  = 44.1E+3;
    sam_tim  = 1 / sam_frq;
    tim      = 1;
    n_fft    = std::pow(2, (11 + 1)) - 1;
    win_size = std::floor(100E-3 / (1 / sam_frq));
    hop_size = std::ceil(win_size / std::pow(2, 4));
    blak_man_win_fcn(win, win_size);

    sim_tim = arma::regspace<type::vec>(0, sam_tim, tim);
    sim_in.resize(sim_tim.n_elem, 8);
    sim_in.col(0) = arma::sin(2 * cnst::mat::pi * 2E+2 * sim_tim);
    sim_in.col(1) = arma::sin(2 * cnst::mat::pi * 4E+2 * sim_tim);
    sim_in.col(2) = arma::sin(2 * cnst::mat::pi * 6E+2 * sim_tim);
    sim_in.col(3) = arma::sin(2 * cnst::mat::pi * 8E+2 * sim_tim);
    sim_in.col(4) = arma::sin(2 * cnst::mat::pi * 2E+3 * sim_tim);
    sim_in.col(5) = arma::sin(2 * cnst::mat::pi * 4E+3 * sim_tim);
    sim_in.col(6) = arma::sin(2 * cnst::mat::pi * 6E+3 * sim_tim);
    sim_in.col(7) = arma::sin(2 * cnst::mat::pi * 8E+3 * sim_tim);

    cap_tim cap_tim_obj;
    cap_tim_obj.tic();
    stft_slice(out_stft, out_frq, out_tim, sim_in, win, sam_frq, n_fft, hop_size);
    cap_tim_obj.toc();
    std::cout << GLOG() << "Around " << cap_tim_obj.get_tim() / 1000.0 << " msec." << std::endl;
    type::mat out_stft_sam = arma::abs(out_stft.slice(7));
    out_stft_sam.save("out/ex_stft_1.csv", arma::csv_ascii);
    out_frq.save("out/ex_stft_2.csv", arma::csv_ascii);
    out_tim.save("out/ex_stft_3.csv", arma::csv_ascii);
    out_stft_sam.reset();
    sim_tim.reset();
    sim_in.reset();
    out_stft.reset();
    out_frq.reset();
    out_tim.reset();
    SUCCEED();
}

#endif
