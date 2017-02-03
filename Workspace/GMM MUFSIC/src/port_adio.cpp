#include <port_adio.h>
#include <type_traits>
#include <ctime>

template <typename t>
port_adio<t>::~port_adio(void) {
    del_buf_recd();
}

template <typename t>
port_adio<t>::port_adio (const type::uint new_id,
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
    cal_lag();
    new_buf_recd();
}

template <typename t>
port_adio<t>::port_adio (const type::uint new_n_ch,
                         const type::uint new_frm_per_buf,
                         const type::uint new_tim_s,
                         const type::uint new_sam_frq) : port_adio(Pa_GetDefaultInputDevice(),
                                                         new_n_ch,
                                                         new_frm_per_buf,
                                                         new_tim_s,
                                                         new_sam_frq) {
}

template <typename t>
void port_adio<t>::cal_max_frm (void) {
    max_frm = tim_s * sam_frq;
}

template <typename t>
void port_adio<t>::cal_sam_fmt (void) {
    if      (std::is_same<t, type::val>::value) {
        in_stem_para->sampleFormat = paFloat32;
    }
    else if (std::is_same<t, type::sint32>::value) {
        in_stem_para->sampleFormat = paInt32;
    }
    else if (std::is_same<t, type::sint16>::value) {
        in_stem_para->sampleFormat = paInt16;
    }
    else if (std::is_same<t, type::sint8>::value) {
        in_stem_para->sampleFormat = paInt8;
    }
    else if (std::is_same<t, type::uint8>::value) {
        in_stem_para->sampleFormat = paUInt8;
    }
    else {
        in_stem_para->sampleFormat = paInt16;
    }
    in_stem_para->sampleFormat = in_stem_para->sampleFormat | paNonInterleaved;
}

template <typename t>
void port_adio<t>::cal_lag (void) {
    in_stem_para->suggestedLatency = Pa_GetDeviceInfo(in_stem_para->device)->defaultLowInputLatency;
}

template <typename t>
void port_adio<t>::new_buf_recd (void) {
    type::uint n_ch = in_stem_para->channelCount;
    buf_recd = type::abs_mat<t>(max_frm, n_ch);
    recd     = type::abs_mat<t>(max_frm, n_ch);
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
void port_adio<t>::save (void) {
    time_t now = time(0);
    std::string str(ctime(&now));
    buf_recd.save(str, arma::csv_ascii);
}

template <typename t>
void port_adio<t>::open_port (void) {
    err = Pa_Initialize();
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
    port_adio<t>* obj = (port_adio<t>*) usr_dat;
    const t**     t_in_buf = (const t**) in_buf;
    type::uint    gap_frm  = obj->max_frm - obj->idx_frm;
    type::uint    n_cpy_frm;
    std::size_t   n_cpy_frm_elem;
    if (gap_frm < call_frm_per_buf) {
        n_cpy_frm = gap_frm;
        resu      = paComplete;
    }
    else {
        n_cpy_frm = call_frm_per_buf;
        resu      = paContinue;
    }
    n_cpy_frm_elem = (std::size_t) n_cpy_frm * obj->elem_size;
    if(in_buf != NULL) {
        //buf_recd.rows(idx_frm, idx_frm + n_cpy_frm)->zeros();
        obj->buf_recd.rows(obj->idx_frm, obj->idx_frm + n_cpy_frm).zeros();
    }
    else {
        for (type::uint idx_ch = 0; idx_ch < obj->in_stem_para->channelCount; idx_ch++) {
            //buf_recd.Col(idx_ch)->row(idx_frm)->memptr()
            std::memcpy(obj->buf_recd.colptr(idx_ch) + obj->idx_frm,
                        &t_in_buf[idx_ch][0],
                        n_cpy_frm_elem);
        }
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
void port_adio<t>::clos_port (void) {
    err = Pa_Terminate();
}

template <typename t>
bool port_adio<t>::is_err (void) {
    return (err != paNoError);
}

template class port_adio<type::val>;
template class port_adio<type::sint32>;
//template class port_adio<type::sint16>;
//template class port_adio<type::sint8>;
template class port_adio<type::uint8>;
