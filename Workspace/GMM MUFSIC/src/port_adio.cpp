#include <port_adio.h>

template <typename t>
port_adio<t>::~port_adio() {
}

template <typename t>
port_adio<t>::port_adio (const aud_t::uint new_id,
                         const aud_t::uint new_n_ch,
                         const aud_t::uint new_frm_per_buf,
                         const aud_t::uint new_tim_s,
                         const aud_t::uint new_sam_frq) {
    idx_frm     = 0;
    id          = new_id;
    n_ch        = new_n_ch;
    frm_per_buf = new_frm_per_buf;
    tim_s       = new_tim_s;
    sam_frq     = new_sam_frq;
    port_adio<t>::cal_lag();
    port_adio<t>::cal_max_frm();
    port_adio<t>::new_buf_recd();
}

template <typename t>
port_adio<t>::port_adio (const aud_t::uint new_n_ch,
                         const aud_t::uint new_frm_per_buf,
                         const aud_t::uint new_tim_s,
                         const aud_t::uint new_sam_frq) : port_adio(Pa_GetDefaultInputDevice(),
                                                                    new_n_ch,
                                                                    new_frm_per_buf,
                                                                    new_tim_s,
                                                                    new_sam_frq) {
}

template <typename t>
void port_adio<t>::cal_lag (void) {
    lag = Pa_GetDeviceInfo(id)->defaultLowInputLatency;
}

template <typename t>
void port_adio<t>::cal_max_frm (void) {
    max_frm = tim_s * sam_frq;
}

template <typename t>
void port_adio<t>::new_buf_recd (void) {
    buf_recd = new t* [n_ch];
    recd     = new t* [n_ch];
    for (aud_t::uint idx_ch = 0; idx_ch < n_ch; idx_ch++) {
        buf_recd = new t [max_frm];
        recd     = new t [max_frm];
    }
}

template <typename t>
void port_adio<t>::del_buf_recd (void) {
    for (aud_t::uint idx_ch = 0; idx_ch < n_ch; idx_ch++) {
        delete buf_recd[max_frm];
        delete recd[max_frm];
    }
    delete[] buf_recd;
    delete[] recd;
}

template <typename t>
void port_adio<t>::re_new_buf_recd (void) {
    port_adio<t>::del_buf_recd();
    port_adio<t>::new_buf_recd();
}

template <typename t>
aud_t::uint port_adio<t>::get_id (void) {
    return id;
}

template <typename t>
aud_t::uint port_adio<t>::get_n_ch (void) {
    return n_ch;
}

template <typename t>
aud_t::fpt port_adio<t>::get_lag (void) {
    return lag;
}

template <typename t>
aud_t::uint port_adio<t>::get_frm_per_buf (void) {
    return frm_per_buf;
}

template <typename t>
aud_t::uint port_adio<t>::get_tim_s (void) {
    return tim_s;
}

template <typename t>
aud_t::fpt port_adio<t>::get_sam_frq (void) {
    return sam_frq;
}

template <typename t>
aud_t::uint port_adio<t>::get_max_frm (void) {
    return max_frm;
}
