#ifdef USING_UNITTEST
//#ifdef DUMMY

#include <gtest/gtest.h>
#include <dsp.h>
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
    type::cx_cub out;
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
        sim_in.set_size(sim_tim.n_elem, 8);
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
        stft_slice(out, out_frq, out_tim, sim_in, win, sam_frq, n_fft, hop_size);
        cap_tim_obj.toc();
        std::cout << GLOG() << "  16K around " << cap_tim_obj.get_tim() / 1000.0 << " msec,"
                            << " with next pow of 2 #fft "   << i << "." <<  std::endl;
        sim_tim.reset();
        sim_in.reset();
        out.reset();
        out_frq.reset();
        out_tim.reset();
    }
    SUCCEED();
}

TEST(TEST_MATH_FFT, STFT_2) {
    type::cx_cub out;
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
        sim_in.set_size(sim_tim.n_elem, 8);
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
        stft_slice(out, out_frq, out_tim, sim_in, win, sam_frq, n_fft, hop_size);
        cap_tim_obj.toc();
        std::cout << GLOG() << "44.1K around " << cap_tim_obj.get_tim() / 1000.0 << " msec,"
                            << " with next pow of 2 #fft "   << i << "." << std::endl;
        sim_tim.reset();
        sim_in.reset();
        out.reset();
        out_frq.reset();
        out_tim.reset();
    }
    SUCCEED();
}

TEST(TEST_MATH_FFT, STFT_3) {
    type::cx_cub out;
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
        sim_in.set_size(sim_tim.n_elem, 8);
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
        stft_slice(out, out_frq, out_tim, sim_in, win, sam_frq, n_fft, hop_size);
        cap_tim_obj.toc();
        std::cout << GLOG() << " 192K around " << cap_tim_obj.get_tim() / 1000.0 << " msec,"
                            << " with next pow of 2 #fft "   << i << "." <<  std::endl;
        sim_tim.reset();
        sim_in.reset();
        out.reset();
        out_frq.reset();
        out_tim.reset();
    }
    SUCCEED();
}

TEST(TEST_MATH_FFT, STFT_ODD) {
    type::cx_cub out;
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
    sim_in.set_size(sim_tim.n_elem, 8);
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
    stft_slice(out, out_frq, out_tim, sim_in, win, sam_frq, n_fft, hop_size);
    cap_tim_obj.toc();
    std::cout << GLOG() << "Around " << cap_tim_obj.get_tim() / 1000.0 << " msec." << std::endl;
    type::mat out_stft_sam = arma::abs(out.slice(7));
    out_stft_sam.save("out/ex_stft_1.csv", arma::csv_ascii);
    out_frq.save("out/ex_stft_2.csv", arma::csv_ascii);
    out_tim.save("out/ex_stft_3.csv", arma::csv_ascii);
    out_stft_sam.reset();
    sim_tim.reset();
    sim_in.reset();
    out.reset();
    out_frq.reset();
    out_tim.reset();
    SUCCEED();
}

TEST(TEST_MATH_FFT, STFT_EVEN) {
    type::cx_cub out;
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
    n_fft    = std::pow(2, (11 + 1));
    win_size = std::floor(100E-3 / (1 / sam_frq));
    hop_size = std::ceil(win_size / std::pow(2, 4));
    blak_man_win_fcn(win, win_size);

    sim_tim = arma::regspace<type::vec>(0, sam_tim, tim);
    sim_in.set_size(sim_tim.n_elem, 8);
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
    stft_slice(out, out_frq, out_tim, sim_in, win, sam_frq, n_fft, hop_size);
    cap_tim_obj.toc();
    std::cout << GLOG() << "Around " << cap_tim_obj.get_tim() / 1000.0 << " msec." << std::endl;
    type::mat out_stft_sam = arma::abs(out.slice(7));
    out_stft_sam.save("out/ex_stft_4.csv", arma::csv_ascii);
    out_frq.save("out/ex_stft_5.csv", arma::csv_ascii);
    out_tim.save("out/ex_stft_6.csv", arma::csv_ascii);
    out_stft_sam.reset();
    sim_tim.reset();
    sim_in.reset();
    out.reset();
    out_frq.reset();
    out_tim.reset();
    SUCCEED();
}

TEST(TEST_MATH_FFT, STFT_ISTFT_ODD) {
    type::cx_cub out_stft;
    type::mat    out_istft;
    type::vec    out_frq;
    type::vec    out_tim_stft;
    type::vec    out_tim;
    type::vec    sim_tim;
    type::mat    sim_in;
    type::vec    win_stft;
    type::vec    win_istft;
    type::val    sam_frq;
    type::val    sam_tim;
    type::val    tim;
    type::uint   n_fft;
    type::uint   hop_size;
    type::uint   win_size_stft;
    type::uint   win_size_istft;

    sam_frq        = 44.1E+3;
    sam_tim        = 1 / sam_frq;
    tim            = 1;
    n_fft          = std::pow(2, (11 + 1)) - 1;
    win_size_stft  = std::floor(100E-3 / (1 / sam_frq));
    win_size_istft = n_fft;
    hop_size       = std::ceil(win_size_stft / std::pow(2, 4));
    blak_man_win_fcn(win_stft, win_size_stft);
    blak_man_win_fcn(win_istft, win_size_istft);

    sim_tim = arma::regspace<type::vec>(0, sam_tim, tim);
    sim_in.set_size(sim_tim.n_elem, 8);
    sim_in.col(0) = arma::sin(2 * cnst::mat::pi * 2E+2 * sim_tim);
    sim_in.col(1) = arma::sin(2 * cnst::mat::pi * 4E+2 * sim_tim);
    sim_in.col(2) = arma::sin(2 * cnst::mat::pi * 6E+2 * sim_tim);
    sim_in.col(3) = arma::sin(2 * cnst::mat::pi * 8E+2 * sim_tim);
    sim_in.col(4) = arma::sin(2 * cnst::mat::pi * 2E+3 * sim_tim);
    sim_in.col(5) = arma::sin(2 * cnst::mat::pi * 4E+3 * sim_tim);
    sim_in.col(6) = arma::sin(2 * cnst::mat::pi * 6E+3 * sim_tim);
    sim_in.col(7) = arma::sin(2 * cnst::mat::pi * 8E+3 * sim_tim);

    stft_slice(out_stft, out_frq, out_tim_stft, sim_in, win_stft, sam_frq, n_fft, hop_size);

    cap_tim cap_tim_obj;
    cap_tim_obj.tic();
    istft_slice(out_istft, out_tim, out_stft, win_istft, sam_frq, n_fft, hop_size);
    cap_tim_obj.toc();
    std::cout << GLOG() << "Around " << cap_tim_obj.get_tim() / 1000.0 << " msec." << std::endl;
    out_istft.save("out/ex_istft_1.csv", arma::csv_ascii);
    out_tim.save("out/ex_istft_2.csv", arma::csv_ascii);

    out_istft.reset();
    out_tim.reset();
    sim_tim.reset();
    sim_in.reset();
    out_stft.reset();
    out_frq.reset();
    out_tim_stft.reset();
    SUCCEED();
}

TEST(TEST_MATH_FFT, STFT_ISTFT_EVEN) {
    type::cx_cub out_stft;
    type::mat    out_istft;
    type::vec    out_frq;
    type::vec    out_tim_stft;
    type::vec    out_tim;
    type::vec    sim_tim;
    type::mat    sim_in;
    type::vec    win_stft;
    type::vec    win_istft;
    type::val    sam_frq;
    type::val    sam_tim;
    type::val    tim;
    type::uint   n_fft;
    type::uint   hop_size;
    type::uint   win_size_stft;
    type::uint   win_size_istft;

    sam_frq        = 44.1E+3;
    sam_tim        = 1 / sam_frq;
    tim            = 1;
    n_fft          = std::pow(2, (11 + 1));
    win_size_stft  = std::floor(100E-3 / (1 / sam_frq));
    win_size_istft = n_fft;
    hop_size       = std::ceil(win_size_stft / std::pow(2, 4));
    blak_man_win_fcn(win_stft, win_size_stft);
    blak_man_win_fcn(win_istft, win_size_istft);

    sim_tim = arma::regspace<type::vec>(0, sam_tim, tim);
    sim_in.set_size(sim_tim.n_elem, 8);
    sim_in.col(0) = arma::sin(2 * cnst::mat::pi * 2E+2 * sim_tim);
    sim_in.col(1) = arma::sin(2 * cnst::mat::pi * 4E+2 * sim_tim);
    sim_in.col(2) = arma::sin(2 * cnst::mat::pi * 6E+2 * sim_tim);
    sim_in.col(3) = arma::sin(2 * cnst::mat::pi * 8E+2 * sim_tim);
    sim_in.col(4) = arma::sin(2 * cnst::mat::pi * 2E+3 * sim_tim);
    sim_in.col(5) = arma::sin(2 * cnst::mat::pi * 4E+3 * sim_tim);
    sim_in.col(6) = arma::sin(2 * cnst::mat::pi * 6E+3 * sim_tim);
    sim_in.col(7) = arma::sin(2 * cnst::mat::pi * 8E+3 * sim_tim);

    stft_slice(out_stft, out_frq, out_tim_stft, sim_in, win_stft, sam_frq, n_fft, hop_size);

    cap_tim cap_tim_obj;
    cap_tim_obj.tic();
    istft_slice(out_istft, out_tim, out_stft, win_istft, sam_frq, n_fft, hop_size);
    cap_tim_obj.toc();
    std::cout << GLOG() << "Around " << cap_tim_obj.get_tim() / 1000.0 << " msec." << std::endl;
    out_istft.save("out/ex_istft_3.csv", arma::csv_ascii);
    out_tim.save("out/ex_istft_4.csv", arma::csv_ascii);

    out_istft.reset();
    out_tim.reset();
    sim_tim.reset();
    sim_in.reset();
    out_stft.reset();
    out_frq.reset();
    out_tim_stft.reset();
    SUCCEED();
}

TEST(TEST_MATH_COV, COV_1) {
    type::mat  in =
        {{0.162182308193243, 0.6892145031400080, 0.53834243526005700, 0.8173032206534330},
         {0.794284540683907, 0.7481515928237100, 0.99613471662688600, 0.8686947053635100},
         {0.311215042044805, 0.4505415985024980, 0.07817552875318370, 0.0844358455109103},
         {0.528533135506213, 0.0838213779969326, 0.44267826977544600, 0.3997826490988970},
         {0.165648729499781, 0.2289769687168190, 0.10665277018058400, 0.2598704028506540},
         {0.601981941401637, 0.9133373615016700, 0.96189808085505400, 0.8000684802243080},
         {0.262971284540144, 0.1523780189692230, 0.00463422413406744, 0.4314138274635450},
         {0.654079098476782, 0.8258169774895470, 0.77491046471150200, 0.9106475944295230}};
    type::vec  avg_in;
    type::uint n_sam;
    type::uint n_var;
    type::val  tem_i;
    type::val  tem_j;
    type::mat  cov_1;
    type::mat  cov_2;
    cap_tim    cap_tim_1_obj;
    cap_tim    cap_tim_2_obj;


    cap_tim_1_obj.tic();
    n_sam  = in.n_rows;
    n_var  = in.n_cols;
    avg_in = arma::mean(in, 0).t();
    cov_1.set_size(n_var, n_var);
    for (type::uint i = 0, j; i < n_var; i++) {
        tem_i = arma::as_scalar(avg_in(i));
        for (j = 0; j < n_var; j++) {
            tem_j = arma::as_scalar(avg_in(j));
            cov_1(i, j) = arma::sum((in.col(i) - tem_i) % (in.col(j) - tem_j));
        }
    }
    cov_1 /= n_sam;
    cap_tim_1_obj.toc();

    cap_tim_2_obj.tic();
    cov_2 = arma::cov(in, 1);
    cap_tim_2_obj.toc();

    std::cout << GLOG() << "Native covariance matrix,    around "
                        << cap_tim_1_obj.get_tim() / 1000.0 << " msec." << std::endl;
    std::cout << GLOG() << "Armadillo covariance matrix, around "
                        << cap_tim_2_obj.get_tim() / 1000.0 << " msec." << std::endl;

    EXPECT_TRUE(arma::approx_equal(cov_1, cov_2, "absdiff", cnst::tol));
}

TEST(TEST_MATH_COV, COV_2) {
    type::mat re_in =
        {{0.162182308193243, 0.6892145031400080, 0.53834243526005700, 0.8173032206534330},
         {0.794284540683907, 0.7481515928237100, 0.99613471662688600, 0.8686947053635100},
         {0.311215042044805, 0.4505415985024980, 0.07817552875318370, 0.0844358455109103},
         {0.528533135506213, 0.0838213779969326, 0.44267826977544600, 0.3997826490988970},
         {0.165648729499781, 0.2289769687168190, 0.10665277018058400, 0.2598704028506540},
         {0.601981941401637, 0.9133373615016700, 0.96189808085505400, 0.8000684802243080},
         {0.262971284540144, 0.1523780189692230, 0.00463422413406744, 0.4314138274635450},
         {0.654079098476782, 0.8258169774895470, 0.77491046471150200, 0.9106475944295230}};
    type::mat im_in =
        {{0.181847028302853, 0.8530311177218940, 0.18390778828241700, 0.3377194098213770},
         {0.263802916521990, 0.6220551314850660, 0.23995252566490300, 0.9000538464176620},
         {0.145538980384717, 0.3509523808922710, 0.41726706908437000, 0.3692467811202150},
         {0.136068558708664, 0.5132495398670530, 0.04965443032574210, 0.1112027552937870},
         {0.869292207640089, 0.4018080337519420, 0.90271610991528100, 0.7802520683211380},
         {0.579704587365570, 0.0759666916908419, 0.94478718972164600, 0.3897388369612530},
         {0.549860201836332, 0.2399161535536580, 0.49086409246808000, 0.2416912859138330},
         {0.144954798223727, 0.1233189348351660, 0.48925263840001900, 0.4039121455881150}};
    type::cx_mat in(re_in.n_rows, re_in.n_cols);
    in.set_real(re_in);
    in.set_imag(im_in);
    type::cx_vec avg_in;
    type::uint   n_sam;
    type::uint   n_var;
    type::cx_val tem_i;
    type::cx_val tem_j;
    type::cx_mat cov_1;
    type::cx_mat cov_2;
    cap_tim      cap_tim_1_obj;
    cap_tim      cap_tim_2_obj;

    cap_tim_1_obj.tic();
    n_sam  = in.n_rows;
    n_var  = in.n_cols;
    avg_in = arma::mean(in, 0).st();
    cov_1.set_size(n_var, n_var);
    for (type::uint i = 0, j; i < n_var; i++) {
        tem_i = arma::as_scalar(avg_in(i));
        for (j = 0; j < n_var; j++) {
            tem_j = arma::as_scalar(avg_in(j));
            cov_1(i, j) = arma::sum((in.col(i) - tem_i) % arma::conj(in.col(j) - tem_j));
        }
    }
    cov_1 = arma::conj(cov_1) / n_sam;
    cap_tim_1_obj.toc();

    cap_tim_2_obj.tic();
    cov_2 = arma::cov(in, 1);
    cap_tim_2_obj.toc();

    std::cout << GLOG() << "Native complex covariance matrix,    around "
                        << cap_tim_1_obj.get_tim() / 1000.0 << " msec." << std::endl;
    std::cout << GLOG() << "Armadillo complex covariance matrix, around "
                        << cap_tim_2_obj.get_tim() / 1000.0 << " msec." << std::endl;

    EXPECT_TRUE(arma::approx_equal(cov_1, cov_2, "absdiff", cnst::tol));
}

TEST(TEST_MATH_COV, COV_3) {
    type::mat  in =
        {{0.162182308193243, 0.6892145031400080, 0.53834243526005700, 0.8173032206534330},
         {0.794284540683907, 0.7481515928237100, 0.99613471662688600, 0.8686947053635100},
         {0.311215042044805, 0.4505415985024980, 0.07817552875318370, 0.0844358455109103},
         {0.528533135506213, 0.0838213779969326, 0.44267826977544600, 0.3997826490988970},
         {0.165648729499781, 0.2289769687168190, 0.10665277018058400, 0.2598704028506540},
         {0.601981941401637, 0.9133373615016700, 0.96189808085505400, 0.8000684802243080},
         {0.262971284540144, 0.1523780189692230, 0.00463422413406744, 0.4314138274635450},
         {0.654079098476782, 0.8258169774895470, 0.77491046471150200, 0.9106475944295230}};
    type::vec  avg_in;
    type::uint n_sam;
    type::uint n_var;
    type::mat  cov_1;
    type::mat  cov_2;
    cap_tim    cap_tim_1_obj;
    cap_tim    cap_tim_2_obj;


    cap_tim_1_obj.tic();
    n_sam  = in.n_rows;
    n_var  = in.n_cols;
    avg_in = arma::mean(in, 0).t();
    cov_1.set_size(n_var, n_var);
    for (type::uint i = 0, j; i < n_var; i++) {
        for (j = 0; j < n_var; j++) {
            cov_1(i, j) = arma::sum(in.col(i) % in.col(j));
        }
    }
    cov_1 /= n_sam;
    cap_tim_1_obj.toc();

    cap_tim_2_obj.tic();
    sam_cov_zero_mean(cov_2, in);
    cap_tim_2_obj.toc();

    std::cout << GLOG() << "Native covariance matrix with zero mean,    around "
                        << cap_tim_1_obj.get_tim() / 1000.0 << " msec." << std::endl;
    std::cout << GLOG() << "Armadillo covariance matrix with zero mean, around "
                        << cap_tim_2_obj.get_tim() / 1000.0 << " msec." << std::endl;

    EXPECT_TRUE(arma::approx_equal(cov_1, cov_2, "absdiff", cnst::tol));
}

TEST(TEST_MATH_COV, COV_4) {
    type::mat re_in =
        {{0.162182308193243, 0.6892145031400080, 0.53834243526005700, 0.8173032206534330},
         {0.794284540683907, 0.7481515928237100, 0.99613471662688600, 0.8686947053635100},
         {0.311215042044805, 0.4505415985024980, 0.07817552875318370, 0.0844358455109103},
         {0.528533135506213, 0.0838213779969326, 0.44267826977544600, 0.3997826490988970},
         {0.165648729499781, 0.2289769687168190, 0.10665277018058400, 0.2598704028506540},
         {0.601981941401637, 0.9133373615016700, 0.96189808085505400, 0.8000684802243080},
         {0.262971284540144, 0.1523780189692230, 0.00463422413406744, 0.4314138274635450},
         {0.654079098476782, 0.8258169774895470, 0.77491046471150200, 0.9106475944295230}};
    type::mat im_in =
        {{0.181847028302853, 0.8530311177218940, 0.18390778828241700, 0.3377194098213770},
         {0.263802916521990, 0.6220551314850660, 0.23995252566490300, 0.9000538464176620},
         {0.145538980384717, 0.3509523808922710, 0.41726706908437000, 0.3692467811202150},
         {0.136068558708664, 0.5132495398670530, 0.04965443032574210, 0.1112027552937870},
         {0.869292207640089, 0.4018080337519420, 0.90271610991528100, 0.7802520683211380},
         {0.579704587365570, 0.0759666916908419, 0.94478718972164600, 0.3897388369612530},
         {0.549860201836332, 0.2399161535536580, 0.49086409246808000, 0.2416912859138330},
         {0.144954798223727, 0.1233189348351660, 0.48925263840001900, 0.4039121455881150}};
    type::cx_mat in(re_in.n_rows, re_in.n_cols);
    in.set_real(re_in);
    in.set_imag(im_in);
    type::cx_vec avg_in;
    type::uint   n_sam;
    type::uint   n_var;
    type::cx_mat cov_1;
    type::cx_mat cov_2;
    cap_tim      cap_tim_1_obj;
    cap_tim      cap_tim_2_obj;

    cap_tim_1_obj.tic();
    n_sam  = in.n_rows;
    n_var  = in.n_cols;
    avg_in = arma::mean(in, 0).st();
    cov_1.set_size(n_var, n_var);
    for (type::uint i = 0, j; i < n_var; i++) {
        for (j = 0; j < n_var; j++) {
            cov_1(i, j) = arma::sum(in.col(i) % arma::conj(in.col(j)));
        }
    }
    cov_1 = arma::conj(cov_1) / n_sam;
    cap_tim_1_obj.toc();

    cap_tim_2_obj.tic();
    sam_cov_zero_mean(cov_2, in);
    cap_tim_2_obj.toc();

    std::cout << GLOG() << "Native complex covariance matrix with zero mean,    around "
                        << cap_tim_1_obj.get_tim() / 1000.0 << " msec." << std::endl;
    std::cout << GLOG() << "Armadillo complex covariance matrix with zero mean, around "
                        << cap_tim_2_obj.get_tim() / 1000.0 << " msec." << std::endl;

    EXPECT_TRUE(arma::approx_equal(cov_1, cov_2, "absdiff", cnst::tol));
}

TEST(TEST_MATH_COV, COV_5) {
    type::mat  in = arma::randu<type::mat>(44101, 8);
    type::vec  avg_in;
    type::uint n_sam;
    type::uint n_var;
    type::mat  cov_1;
    type::mat  cov_2;
    cap_tim    cap_tim_1_obj;
    cap_tim    cap_tim_2_obj;


    cap_tim_1_obj.tic();
    n_sam  = in.n_rows;
    n_var  = in.n_cols;
    avg_in = arma::mean(in, 0).t();
    cov_1.set_size(n_var, n_var);
    for (type::uint i = 0, j; i < n_var; i++) {
        for (j = 0; j < n_var; j++) {
            cov_1(i, j) = arma::sum(in.col(i) % in.col(j));
        }
    }
    cov_1 /= n_sam;
    cap_tim_1_obj.toc();

    cap_tim_2_obj.tic();
    sam_cov_zero_mean(cov_2, in);
    cap_tim_2_obj.toc();

    std::cout << GLOG() << "Native covariance matrix with zero mean,    around "
                        << cap_tim_1_obj.get_tim() / 1000.0 << " msec." << std::endl;
    std::cout << GLOG() << "Armadillo covariance matrix with zero mean, around "
                        << cap_tim_2_obj.get_tim() / 1000.0 << " msec." << std::endl;

    EXPECT_TRUE(arma::approx_equal(cov_1, cov_2, "absdiff", cnst::tol));
}

TEST(TEST_MATH_COV, COV_6) {
    type::cx_mat in = arma::randu<type::cx_mat>(44101, 8);
    type::cx_vec avg_in;
    type::uint   n_sam;
    type::uint   n_var;
    type::cx_mat cov_1;
    type::cx_mat cov_2;
    cap_tim      cap_tim_1_obj;
    cap_tim      cap_tim_2_obj;

    cap_tim_1_obj.tic();
    n_sam  = in.n_rows;
    n_var  = in.n_cols;
    avg_in = arma::mean(in, 0).st();
    cov_1.set_size(n_var, n_var);
    for (type::uint i = 0, j; i < n_var; i++) {
        for (j = 0; j < n_var; j++) {
            cov_1(i, j) = arma::sum(in.col(i) % arma::conj(in.col(j)));
        }
    }
    cov_1 = arma::conj(cov_1) / n_sam;
    cap_tim_1_obj.toc();

    cap_tim_2_obj.tic();
    sam_cov_zero_mean(cov_2, in);
    cap_tim_2_obj.toc();

    std::cout << GLOG() << "Native complex covariance matrix with zero mean,    around "
                        << cap_tim_1_obj.get_tim() / 1000.0 << " msec." << std::endl;
    std::cout << GLOG() << "Armadillo complex covariance matrix with zero mean, around "
                        << cap_tim_2_obj.get_tim() / 1000.0 << " msec." << std::endl;

    EXPECT_TRUE(arma::approx_equal(cov_1, cov_2, "absdiff", cnst::tol));
}

TEST(TEST_MATH_ROOT, ROOT_1) {
    type::vec    in  = {1, 2, 1};
    type::cx_vec ans = {type::cx_val(-1, 0), type::cx_val(-1, 0)};
    type::cx_vec out;
    root(out, in);
    EXPECT_TRUE(arma::approx_equal(out, ans, "absdiff", cnst::tol));
}

TEST(TEST_MATH_ROOT, ROOT_2) {
    type::vec    in  = {1, 1, 1};
    type::cx_vec ans = {type::cx_val(-0.5, +0.866025403784439),
                        type::cx_val(-0.5, -0.866025403784439)};
    type::cx_vec out;
    root(out, in);
    EXPECT_TRUE(arma::approx_equal(out, ans, "absdiff", cnst::tol));
}

TEST(TEST_MATH_ROOT, ROOT_3) {
    type::cx_vec in  = {type::cx_val(1, 0),  type::cx_val(2, 0), type::cx_val(1, 0)};
    type::cx_vec ans = {type::cx_val(-1, 0), type::cx_val(-1, 0)};
    type::cx_vec out;
    root(out, in);
    EXPECT_TRUE(arma::approx_equal(out, ans, "absdiff", cnst::tol));
}

TEST(TEST_MATH_ROOT, ROOT_4) {
    type::cx_vec in  = {type::cx_val(1, 1),  type::cx_val(2, 1), type::cx_val(1, 1)};
    type::cx_vec ans = {type::cx_val(-1, 1), type::cx_val(-0.5, -0.5)};
    type::cx_vec out;
    root(out, in);
    EXPECT_TRUE(arma::approx_equal(out, ans, "absdiff", cnst::tol));
}

#endif
