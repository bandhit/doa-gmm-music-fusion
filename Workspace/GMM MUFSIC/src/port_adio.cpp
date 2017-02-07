#include <port_adio.h>
#include <type_traits>
#include <ctime>

pa_err open_port (void) {
    freopen("/dev/null", "w", stderr); // PortAudio log bug
    pa_err err = Pa_Initialize();
    freopen("/dev/tty", "w", stderr); // reopen global console
    return err;
}

pa_err clos_port (void) {
    return Pa_Terminate();
}

std::string get_err_str (pa_err err) {
    return std::string(Pa_GetErrorText(err));
}

template <typename t>
port_adio<t>::~port_adio(void) {
    del_buf_recd();
}

template <typename t>
port_adio<t>::port_adio (const type::sint new_id,
                         const type::uint new_n_ch,
                         const type::uint new_frm_per_buf,
                         const type::uint new_tim_s,
                         const type::uint new_sam_frq) {
    idx_frm      = 0;
    frm_per_buf  = new_frm_per_buf;
    sam_frq      = new_sam_frq;
    tim_s        = new_tim_s;
    cal_max_frm();
    in_stem_para = new PaStreamParameters;
    in_stem_para->device                    = new_id;
    in_stem_para->channelCount              = new_n_ch;
    in_stem_para->hostApiSpecificStreamInfo = NULL;
    cal_sam_fmt();
    if (new_id >= 0) {
        cal_lag();
    }
    new_buf_recd();

}

template <typename t>
port_adio<t>::port_adio (const type::uint new_n_ch,
                         const type::uint new_frm_per_buf,
                         const type::uint new_tim_s,
                         const type::uint new_sam_frq) : port_adio(get_id(),
                                                                   new_n_ch,
                                                                   new_frm_per_buf,
                                                                   new_tim_s,
                                                                   new_sam_frq) {
}

template <typename t>
port_adio<t>::port_adio (const std::string ch_str,
                         const type::uint  new_n_ch,
                         const type::uint  new_frm_per_buf,
                         const type::uint  new_tim_s,
                         const type::uint  new_sam_frq) : port_adio(get_id(ch_str),
                                                                    new_n_ch,
                                                                    new_frm_per_buf,
                                                                    new_tim_s,
                                                                    new_sam_frq) {
}

template <typename t>
type::sint port_adio<t>::get_id (void) {
    return Pa_GetDefaultInputDevice();
}

template <typename t>
type::sint port_adio<t>::get_id (const std::string str) {
    std::string tmp_str;
    type::uint  str_len;
    type::sint  new_id = -1;
    type::uint  n_id;
    str_len = str.length();
    n_id    = Pa_GetDeviceCount();
    for (type::uint idx_id = 0; idx_id < n_id; idx_id++) {
        tmp_str = Pa_GetDeviceInfo(idx_id)->name;
        if (tmp_str.length() >= str_len) {
            if (str.compare(0, str_len, tmp_str, 0, str_len) == 0) {
                new_id = idx_id;
                break;
            }
            else {
                continue;
            }
        }
        else {
            continue;
        }
    }
    return new_id;
}

template <typename t>
void port_adio<t>::cal_max_frm (void) {
    max_frm = tim_s * sam_frq;
}

template <typename t>
void port_adio<t>::cal_sam_fmt (void) {
    if      (std::is_same<t, type::val>::value) {
        in_stem_para->sampleFormat = paFloat32;
        post_pos_amp               = 1;
        post_pos_off               = 0;
    }
    else if (std::is_same<t, type::sint32>::value) {
        in_stem_para->sampleFormat = paInt32;
        post_pos_amp               = (type::val) std::pow(2, 31);
        post_pos_off               = 0;
    }
    // paInt24: Linux system unsupported type | lazy to implement it
    else if (std::is_same<t, type::sint16>::value) {
        in_stem_para->sampleFormat = paInt16;
        post_pos_amp               = (type::val) std::pow(2, 15);
        post_pos_off               = 0;
    }
    else if (std::is_same<t, type::sint8>::value) {
        in_stem_para->sampleFormat = (type::val) std::pow(2, 7);
        post_pos_off               = 0;
    }
    else if (std::is_same<t, type::uint8>::value) {
        in_stem_para->sampleFormat = paUInt8;
        in_stem_para->sampleFormat = (type::val) std::pow(2, 7);
        post_pos_off               = 128;
    }
    in_stem_para->sampleFormat |= paNonInterleaved;
}

template <typename t>
void port_adio<t>::cal_lag (void) {
    in_stem_para->suggestedLatency = Pa_GetDeviceInfo(in_stem_para->device)->defaultLowInputLatency;
}

template <typename t>
void port_adio<t>::new_buf_recd (void) {
    buf_recd.resize(max_frm, in_stem_para->channelCount);
    recd.resize(max_frm, in_stem_para->channelCount);
}

template <typename t>
void port_adio<t>::del_buf_recd (void) {
    buf_recd.reset();
    recd.reset();
}

template <typename t>
void port_adio<t>::re_new_buf_recd (void) {
    del_buf_recd();
    new_buf_recd();
}

template <typename t>
void port_adio<t>::move_buf_recd (void) {
    #ifndef PORT_ADIO_USING_STD_MEMCPY
    recd = buf_recd;
    #else
    std::memcpy(recd.memptr(),
                buf_recd.memptr(),
                (std::size_t) max_frm * in_stem_para->channelCount * elem_size);
    #endif
}

template <typename t>
const type::abs_mat<t>& port_adio<t>::get_recd (void) {
    return recd;
}

template <typename t>
void port_adio<t>::save (const std::string path) {
    std::time_t raw_tim = std::time(NULL);
    char        buf_bstr[13];
    std::strftime(buf_bstr, sizeof(buf_bstr), "%y%m%d%H%M%S", std::localtime(&raw_tim));
    buf_recd.save(path + buf_bstr + "-BUF" + ".csv", arma::csv_ascii);
    recd    .save(path + buf_bstr + "-REC" + ".csv", arma::csv_ascii);
}

template <typename t>
void port_adio<t>::open_stem (void) {
    err = Pa_OpenStream(&stem,
                        in_stem_para,
                        NULL,
                        sam_frq,
                        frm_per_buf,
                        paClipOff,
                        recd_call_back<t>,
                        this);
}

template <typename t>
int recd_call_back (const void*                     in_buf,
                    void*                           out_buf,
                    unsigned long                   call_frm_per_buf,
                    const PaStreamCallbackTimeInfo* tim_info,
                    PaStreamCallbackFlags           flag,
                    void*                           usr_dat) {
    PaStreamCallbackResult resu;
    port_adio<t>* obj      = (port_adio<t>*) usr_dat;
    const t**     t_in_buf = (const t**) in_buf;
    type::uint    gap_frm  = obj->max_frm - obj->idx_frm;
    type::uint    n_cpy_frm;
    if (gap_frm < call_frm_per_buf) {
        n_cpy_frm = gap_frm;
        resu      = paComplete;
    }
    else {
        n_cpy_frm = call_frm_per_buf;
        resu      = paContinue;
    }
    if(in_buf != NULL) {
        #ifndef PORT_ADIO_USING_STD_MEMCPY
        obj->buf_recd.rows(obj->idx_frm,
                           obj->idx_frm + n_cpy_frm - 1) = type::abs_mat<t>(&t_in_buf[0][0],
                                                                            n_cpy_frm,
                                                                            obj->in_stem_para->channelCount);
        #else
        std::size_t n_cpy_frm_elem = (std::size_t) n_cpy_frm * obj->elem_size;
        for (type::uint idx_ch = 0; idx_ch < obj->in_stem_para->channelCount; idx_ch++) {
            std::memcpy(obj->buf_recd.colptr(idx_ch) + obj->idx_frm,
                        &t_in_buf[idx_ch][0],
                        n_cpy_frm_elem);
        }
        #endif
    }
    else {
        obj->buf_recd.rows(obj->idx_frm, obj->idx_frm + n_cpy_frm - 1).zeros();
    }
    obj->idx_frm += n_cpy_frm;
    return resu;
}

template <typename t>
void port_adio<t>::strt_stem (void) {
    err = Pa_StartStream(stem);
}

template <typename t>
bool port_adio<t>::is_stem_acti (void) {
    return (Pa_IsStreamActive(stem) == 1);
}

template <typename t>
void port_adio<t>::clos_stem (void) {
    err = Pa_CloseStream(stem);
}

template <typename t>
bool port_adio<t>::is_err (void) {
    return (err != pa_no_err);
}

template <typename t>
pa_err port_adio<t>::get_err (void) {
    return err;
}

template <typename t>
type::uint port_adio<t>::get_max_frm (void) {
    return max_frm;
}

template <typename t>
type::uint port_adio<t>::get_n_ch (void) {
    return in_stem_para->channelCount;
}

template <typename t>
type::val port_adio<t>::get_post_pos_amp (void) {
    return post_pos_amp;
}

template <typename t>
type::val port_adio<t>::get_post_pos_off (void) {
    return post_pos_off;
}

template class port_adio<type::val>;
template class port_adio<type::sint32>;
template class port_adio<type::sint16>;
//template class port_adio<type::sint8>; // Armadillo C++ (<= v7.6) unsupported type
template class port_adio<type::uint8>;
