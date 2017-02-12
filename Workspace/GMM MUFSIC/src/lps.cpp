#include <lps.h>
#include <dsp.h>
#include <complex>
#include <type_traits>

void frwd_w_fcn (type::cx_mat&        out,
                 const type::row_vec& phi_rad,
                 const type::uint     n_sam,
                 const type::val      d,
                 const type::val      fun_frq,
                 const type::val      frq,
                 const type::val      lamb) {
    out = arma::exp((cnst::i * (type::val) 2 * cnst::mat::pi * (type::val) d * (frq / (fun_frq * lamb))) *
                    (arma::regspace<type::vec>(0, 1, n_sam - 1) * arma::sin(phi_rad)));
}

void frwd_w_fcn (type::cx_vec&    out,
                 const type::val  phi_rad,
                 const type::uint n_sam,
                 const type::val  d,
                 const type::val  fun_frq,
                 const type::val  frq,
                 const type::val  lamb) {
    out = arma::exp((cnst::i * (type::val) 2 * cnst::mat::pi * (type::val) d * (frq / (fun_frq * lamb))) *
                    (arma::regspace<type::vec>(0, 1, n_sam - 1) * std::sin(phi_rad)));
}

void revs_w_fcn (type::vec&         phi_rad,
                 const type::cx_vec rot,
                 const type::val    d,
                 const type::val    fun_frq,
                 const type::val    frq,
                 const type::val    lamb) {
    phi_rad = arma::asin(- arma::arg(rot) /
                         ((type::val) 2 * cnst::mat::pi * (type::val) d * (frq / (fun_frq * lamb))));
}

void music_min_eig (type::vec&          out,
                    const type::cx_mat& in,
                    const type::vec&    phi_rad,
                    const type::val     d,
                    const type::val     c,
                    const type::val     fun_frq,
                    const type::val     frq,
                    const type::val     min_eig) {
    type::uint     n_sam = in.n_cols;
    type::cx_mat   cov;
    sam_cov_zero_mean(cov, in);
    type::cx_vec   eig_val;
    type::cx_mat   eig_vec;
    arma::eig_gen(eig_val, eig_vec, cov);
    type::vec      re_eig_val = arma::real(eig_val);
    type::uint     n_src      = n_sam - arma::sum(re_eig_val < min_eig);
    type::cx_vec   music_vec(phi_rad.n_elem);
    if (n_src > 0) {
        type::val    lamb = c / fun_frq;
        eig_vec           = eig_vec.cols(arma::sort_index(re_eig_val, "descend"));
        type::cx_mat e    = eig_vec.cols(n_src, n_sam - 1);
        e = e * e.t();
        type::cx_vec     w;
        type::cx_row_vec w_t;
        for(type::uint i = 0; i < phi_rad.n_elem; i++) {
            frwd_w_fcn(w, phi_rad(i), n_sam, d, fun_frq, frq, lamb);
            w_t = w.t();
            music_vec(i) = arma::as_scalar(w_t * w)
                         / arma::as_scalar(w_t * e * w);
        }
    }
    else {
        music_vec.zeros();
    }
    out = arma::abs(music_vec);
}

void root_music_min_eig (type::vec&          phi_rad,
                         type::vec&          mag_rot,
                         const type::cx_mat& in,
                         const type::val     d,
                         const type::val     c,
                         const type::val     fun_frq,
                         const type::val     frq,
                         const type::val     min_eig) {
    type::uint     n_sam = in.n_cols;
    type::cx_mat   cov;
    sam_cov_zero_mean(cov, in);
    type::cx_vec   eig_val;
    type::cx_mat   eig_vec;
    arma::eig_gen(eig_val, eig_vec, cov);
    type::vec      re_eig_val = arma::real(eig_val);
    type::uint     n_src      = n_sam - arma::sum(re_eig_val < min_eig);
    if (n_src > 0) {
        type::val    lamb  = c / fun_frq;
        eig_vec            = eig_vec.cols(arma::sort_index(re_eig_val, "descend"));
        type::cx_mat e     = eig_vec.cols(n_src, n_sam - 1);
        type::cx_mat c     = e * e.t();
        type::uint   n_dim = (2 * n_sam) - 1;
        type::sint   k;
        type::cx_vec p(n_dim);
        for (type::uint i = 0; i < n_dim; i++) {
            k = (type::sint) i + 1 - (type::sint) n_sam;
            p(i) = arma::sum(c.diag(k));
        }
        type::cx_vec rot;
        root(rot, p);
        mag_rot                = arma::abs(arma::abs(rot) - 1);
        type::idx_vec idx      = arma::sort_index(mag_rot, "ascend");
        rot                    = rot(idx);
        mag_rot                = mag_rot(idx);
        type::idx_vec skip_idx = arma::regspace<type::idx_vec>(2, 2, 2 * n_src) - 1;
        rot                    = rot(skip_idx);
        mag_rot                = mag_rot(skip_idx);
        revs_w_fcn(phi_rad, rot, d, fun_frq, frq, lamb);
        type::idx_vec fini_idx = find_finite(phi_rad);
        phi_rad                = phi_rad(fini_idx);
        mag_rot                = mag_rot(fini_idx);
    }
    else {
        phi_rad.reset();
        mag_rot.reset();
    }
}
