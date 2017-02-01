#ifndef _PORT_ADIO_H_
#define _PORT_ADIO_H_

#include <iostream>
#include <fstream>
#include <portaudio.h>
#include <type_def.hpp>

template <typename t>
class port_adio {
    public:
                       port_adio       (const aud_t::uint new_id,
                                        const aud_t::uint new_n_ch,
                                        const aud_t::uint new_frm_per_buf,
                                        const aud_t::uint new_tim_s,
                                        const aud_t::uint new_sam_frq);
                       port_adio       (const aud_t::uint new_n_ch,
                                        const aud_t::uint new_frm_per_buf,
                                        const aud_t::uint new_tim_s,
                                        const aud_t::uint new_sam_frq);
        virtual       ~port_adio       (void);
        void           new_buf_recd    (void);
        void           re_new_buf_recd (void);
        // copy record to any buffer
        //void           open_port       (void);
        //void           open_stem       (void);
        //void           strt_stem       (void);
        //void           clos_stem       (void);
        //bool           is_stem_acti    (void);
        //void           opn_stem        (void);
        //void           clos_port       (void);
        bool           is_err          (void);
        std::string    get_err         (void);
        // find id by string
        aud_t::uint    get_id          (void);
        aud_t::uint    get_n_ch        (void);
        aud_t::fpt     get_lag         (void);
        aud_t::uint    get_frm_per_buf (void);
        aud_t::uint    get_tim_s       (void);
        aud_t::fpt     get_sam_frq     (void);
        aud_t::uint    get_max_frm     (void);
    protected:
        void           cal_lag         (void);
        void           cal_max_frm     (void);
        void           del_buf_recd    (void);
        PaSampleFormat fmt;
        PaDeviceIndex  id;
        int            n_ch;
        PaTime         lag;
        unsigned long  frm_per_buf;
        double         sam_frq;
        aud_t::uint    tim_s;
        aud_t::uint    max_frm;
        aud_t::uint    idx_frm;
        t**            buf_recd;
        t**            recd;
};

#endif
