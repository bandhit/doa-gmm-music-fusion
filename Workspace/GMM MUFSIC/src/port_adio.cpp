#include <port_adio.h>

template <typename t>
port_adio<t>::port_adio(void) {
    idx_frm = 0;
}

template <typename t>
port_adio<t>::~port_adio() {
}

template <typename t>
port_adio<t>::port_adio (const aud_t::uint new_id,
                         const aud_t::uint new_n_ch,
                         const aud_t::uint new_frm_per_buf,
                         const aud_t::uint new_tim_s,
                         const aud_t::uint new_sam_frq) {
    idx_frm = 0;
    port_adio<t>::set_id(new_id);
    port_adio<t>::set_n_ch(new_n_ch);
    port_adio<t>::set_frm_per_buf(new_frm_per_buf);
    port_adio<t>::set_tim_s(new_tim_s);
    port_adio<t>::set_sam_frq(new_sam_frq);
    port_adio<t>::cal_lag();
    port_adio<t>::cal_max_frm();
}

template <typename t>
aud_t::fpt port_adio<t>::cal_lag (void) {
    lag = Pa_GetDeviceInfo(id)->defaultLowInputLatency;
    return lag;
}

template <typename t>
aud_t::uint port_adio<t>::cal_max_frm (void) {
    max_frm = tim_s * sam_frq;
    return max_frm;
}

template <typename t>
void port_adio<t>::set_id (const aud_t::uint new_id) {
    id = new_id;
}

template <typename t>
void port_adio<t>::set_dfut_id (void) {
    id = Pa_GetDefaultInputDevice();
}

template <typename t>
void port_adio<t>::set_n_ch (const aud_t::uint new_n_ch) {
    n_ch = new_n_ch;
}

template <typename t>
void port_adio<t>::set_lag (const aud_t::fpt new_lag) {
    lag = new_lag;
}

template <typename t>
void port_adio<t>::set_frm_per_buf (const aud_t::uint new_frm_per_buf) {
    frm_per_buf = new_frm_per_buf;
}

template <typename t>
void port_adio<t>::set_tim_s (const aud_t::uint new_tim_s) {
    tim_s = new_tim_s;
}

template <typename t>
void port_adio<t>::set_sam_frq (const aud_t::fpt new_sam_frq) {
    sam_frq = new_sam_frq;
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
