#ifndef USING_UNITTEST

#include <iostream>
#include <signal.h>

#include <port_adio.h>
#include <lps.h>
#include <dsp.h>
#include <misc.h>

bool is_end = false;

int main (void) {
    // common configuration
    type::val   sam_frq       = 44.1E+3;
    // interface configuration
    std::string ch_str        = "UA-101: USB Audio";
    type::uint  n_ch          = 2;
    type::uint  frm_per_buf   = 16;
    type::uint  tim_s         = 1;
    using       adio_type     = type::sint32;
    // sfft configuration
    type::uint  n_fft         = std::pow(2, (12 + 1));
    type::uint  win_size_stft = std::floor(100E-3 / (1 / sam_frq));
    type::uint  hop_size      = std::ceil(win_size_stft / std::pow(2, 4));
    type::vec   win_stft;
    blak_man_win_fcn(win_stft, win_size_stft);
    type::val   sum_win_stft  = arma::sum(win_stft);
    // root music configuration
    type::val   fun_frq       = 3.4E+3;
    type::val   c             = 340;
    type::val   d             = (c / fun_frq) / 2; // 5 cm
    type::val   min_eig       = 1E-3;
    type::val   roi_min_frq   = fun_frq / 6.8; // 500 Hz
    type::val   roi_max_frq   = fun_frq / (3.4 / 3.0); // 3 kHz
    type::val   amp           = 100.0;
    // gmm configuration
    type::uint  n_gaus        = 1;
    type::uint  km_iter       = 1000;
    type::uint  em_iter       = 1000;
    type::val   var           = 1E-6;
    // initialization
    bool        is_err        = false;
    std::string err_str;
    cap_tim     cap_tim_obj;
    signal(SIGINT, [](int sig){is_end = true;});
    // main
    pa_err err_open = open_port();
    if (err_open == pa_no_err) {
        try {
            port_adio<adio_type>* obj = new port_adio<adio_type>(ch_str, n_ch, frm_per_buf, tim_s, sam_frq);
            obj->open_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            obj->strt_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            cap_tim_obj.tic();
            while (is_end != true) {
                cap_tim_obj.toc();
                type::val calc_tim_s = cap_tim_obj.get_tim() / 1000000.0;
                if (calc_tim_s > tim_s) {
                    std::cout << "warning, calculation time is more than capture time, "
                              << calc_tim_s << " sec." << std::endl;
                }
                while(obj->is_stem_acti());
                cap_tim_obj.tic();
                obj->clos_stem();
                if (obj->is_err()) {
                    throw get_err_str(obj->get_err());
                }
                obj->move_buf_recd();
                obj->open_stem();
                if (obj->is_err()) {
                    throw get_err_str(obj->get_err());
                }
                obj->strt_stem();
                if (obj->is_err()) {
                    throw get_err_str(obj->get_err());
                }
                type::mat in = arma::conv_to<type::mat>::from(obj->get_recd());
                in -= obj->get_post_pos_off();
                in /= obj->get_post_pos_amp();
                in *= amp;

                type::cx_cub out_stft;
                type::vec    out_frq;
                type::vec    out_tim_stft;
                stft_slice(out_stft, out_frq, out_tim_stft, in, win_stft, sam_frq, n_fft, hop_size);
                out_stft /= sum_win_stft;

                type::uint    t_tim   = out_stft.n_cols - 1;
                type::idx_vec min_idx = arma::find(out_frq >= roi_min_frq, 1, "first");
                type::idx_vec max_idx = arma::find(out_frq <= roi_max_frq, 1, "last");
                type::uint    h_idx   = !min_idx.is_empty() ? min_idx(0) : 0;
                type::uint    t_idx   = !max_idx.is_empty() ? max_idx(0) : t_tim;
                type::uint    nfrq    = t_idx - h_idx + 1;
                type::uint    h_frq;
                type::cx_mat  tmp_out_stft;
                type::vec     tmp_phi_rad;
                type::vec     tmp_mag_rot;
                type::vec     phi_rad;
                for (type::uint i = 0; i < nfrq; i++) {
                    h_frq        = h_idx + i;
                    tmp_out_stft = arma::reshape(out_stft.tube(h_frq, 0, h_frq, t_tim),
                                                 out_stft.n_cols, out_stft.n_slices, 1);
                    root_music_min_eig(tmp_phi_rad, tmp_mag_rot, tmp_out_stft, d, c, fun_frq, out_frq(h_frq), min_eig);
                    if (!tmp_phi_rad.is_empty()) {
                        phi_rad = arma::join_cols(phi_rad, tmp_phi_rad);
                    }
                }

                if (!phi_rad.is_empty()) {
                    type::vec phi_deg = phi_rad;
                    phi_deg *= (type::val) 180;
                    phi_deg /= cnst::mat::pi;
                    arma::gmm_priv::gmm_diag<type::val> model;
                    model.learn(phi_deg.t(), n_gaus, arma::maha_dist, arma::static_spread, km_iter, em_iter, var, false);
                    std::cout << model.means;
                }
                else {
                    std::cout << "[skipped]" << std::endl;
                }
            }
        }
        catch (std::exception& e) {
            is_err  = true;
            err_str = e.what();
        }
        catch (std::string& e) {
            is_err  = true;
            err_str = e;
        }
        pa_err err_clos = clos_port();
        if ((is_err == false) && (err_clos != pa_no_err)) {
            is_err  = true;
            err_str = get_err_str(err_clos);
        }
    }
    else {
        is_err  = true;
        err_str = get_err_str(err_open);
    }
    if (is_err) {
        std::cout << err_str << std::endl;
    }
    return 0;
}

#else

#include <gtest/gtest.h>

GTEST_API_ int main (int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#endif
