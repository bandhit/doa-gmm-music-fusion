#ifndef _PORT_ADIO_H_
#define _PORT_ADIO_H_

#include <portaudio.h>
#include <type_def.hpp>

bool open_port (void);
bool clos_port (void);

template <typename t>
int recd_call_back  (const void*                     in_buf,
                     void*                           out_buf,
                     unsigned long                   call_frm_per_buf,
                     const PaStreamCallbackTimeInfo* tim_info,
                     PaStreamCallbackFlags           flag,
                     void*                           usr_dat);

template <typename t>
class port_adio {
    public:
        virtual       ~port_adio         (void);
                       port_adio         (const type::sint new_id,
                                          const type::uint new_n_ch,
                                          const type::uint new_frm_per_buf,
                                          const type::uint new_tim_s,
                                          const type::uint new_sam_frq);
                       port_adio         (const type::uint new_n_ch,
                                          const type::uint new_frm_per_buf,
                                          const type::uint new_tim_s,
                                          const type::uint new_sam_frq);
        void           new_buf_recd      (void);
        void           re_new_buf_recd   (void);
        // copy record to any buffer
        void           save              (void);
        // save with some path
        void           open_stem         (void);
        void           strt_stem         (void);
        bool           is_stem_acti      (void);
        void           clos_stem         (void);
        bool           is_err            (void);
        std::string    get_err           (void);
        // find id by string
    protected:
        void           cal_max_frm       (void);
        void           cal_sam_fmt       (void);
        void           cal_lag           (void);
        void           del_buf_recd      (void);
        friend int     recd_call_back<t> (const void*                     in_buf,
                                          void*                           out_buf,
                                          unsigned long                   call_frm_per_buf,
                                          const PaStreamCallbackTimeInfo* tim_info,
                                          PaStreamCallbackFlags           flag,
                                          void*                           usr_dat);
        PaStream*           stem;
        PaStreamParameters* in_stem_para;
        unsigned long       frm_per_buf;
        double              sam_frq;
        type::uint          tim_s;
        type::uint          max_frm;
        type::uint          idx_frm;
        type::abs_mat<t>    buf_recd;
        type::abs_mat<t>    recd;
        PaError             err;
        const std::size_t   elem_size = sizeof(t);
};

#endif
