#ifdef USING_UNITTEST

/*#include <gtest/gtest.h>
#include <port_adio.h>
#include <misc.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <stack>

TEST(TEST_ORG_PORT_ADIO, INIT_1) {
    PaError err_1;
    PaError err_2;
    freopen("/dev/null", "w", stderr); // PortAudio log bug
    err_1 = Pa_Initialize();
    freopen("/dev/tty", "w", stderr); // reopen global console
    if (err_1 == paNoError) {
        err_2 = Pa_Terminate();
        if (err_2 == paNoError) {
            SUCCEED();
        }
        else {
            FAIL() << Pa_GetErrorText(err_2);
        }
    }
    else {
        FAIL() << Pa_GetErrorText(err_1);
    }
}

typedef short int type_sam;

typedef struct {
    unsigned long idx_frm;
    unsigned long max_frm;
    type_sam*     recd_ch1;
    type_sam*     recd_ch2;
}
paUserData;

int call_back (const void*                     in_buf,
               void*                           out_buf,
               unsigned long                   frm_per_buf,
               const PaStreamCallbackTimeInfo* tim_inf,
               PaStreamCallbackFlags           flag,
               void*                           usr_dat) {
    PaStreamCallbackResult resu;
    paUserData*     dat       = (paUserData*) usr_dat;
    const type_sam* r_ptr     = (const type_sam*) in_buf;
    type_sam*       w_ptr_ch1 = &dat->recd_ch1[dat->idx_frm];
    type_sam*       w_ptr_ch2 = &dat->recd_ch2[dat->idx_frm];
    unsigned long   gap_frm   = dat->max_frm - dat->idx_frm;
    unsigned long   n_cpy_frm;
    if(gap_frm < frm_per_buf) {
        n_cpy_frm = gap_frm;
        resu      = paComplete;
    }
    else {
        n_cpy_frm = frm_per_buf;
        resu      = paContinue;
    }
    if(in_buf != NULL) {
        for(unsigned long idx_frm = 0; idx_frm < n_cpy_frm; idx_frm++) {
            *w_ptr_ch1++ = *r_ptr++;
            *w_ptr_ch2++ = *r_ptr++;
        }
    }
    else {
        for(unsigned long idx_frm = 0; idx_frm < n_cpy_frm; idx_frm++) {
            *w_ptr_ch1++ = (type_sam) 0;
            *w_ptr_ch2++ = (type_sam) 0;
        }
    }
    dat->idx_frm += n_cpy_frm;
    return resu;
}

TEST(TEST_ORG_PORT_ADIO, HELLO_PA) {
    bool        is_err = false;
    PaError     err_open;
    PaError     err_onwk;
    PaError     err_clos;
    std::string err_str;
    freopen("/dev/null", "w", stderr); // PortAudio log bug
    err_open = Pa_Initialize();
    freopen("/dev/tty", "w", stderr); // reopen global console
    if (err_open == paNoError) {
        try {
            PaStreamParameters in_stem_para;
            in_stem_para.device = Pa_GetDefaultInputDevice();
            in_stem_para.channelCount = 2;
            in_stem_para.sampleFormat = paInt16;
            in_stem_para.suggestedLatency = Pa_GetDeviceInfo(in_stem_para.device)->defaultLowInputLatency;
            in_stem_para.hostApiSpecificStreamInfo = NULL;
            unsigned long frm_per_buf = 512;
            unsigned long tim_s       = 3;
            unsigned long sam_frq     = 44100;
            unsigned long max_frm     = tim_s * sam_frq;
            paUserData* usr_dat = new paUserData;
            usr_dat->idx_frm = 0;
            usr_dat->max_frm = max_frm;
            usr_dat->recd_ch1 = (type_sam*) malloc(max_frm * sizeof(type_sam));
            usr_dat->recd_ch2 = (type_sam*) malloc(max_frm * sizeof(type_sam));
            for(unsigned long idx_frm = 0; idx_frm < max_frm; idx_frm++) {
                usr_dat->recd_ch1[idx_frm] = 0;
                usr_dat->recd_ch2[idx_frm] = 0;
            }
            PaStream* stem;
            err_onwk = Pa_OpenStream(&stem,
                                     &in_stem_para,
                                     NULL,
                                     sam_frq,
                                     frm_per_buf,
                                     paClipOff,
                                     call_back,
                                     usr_dat);
            if (err_onwk != paNoError) {
                throw std::string(Pa_GetErrorText(err_onwk));
            }
            err_onwk = Pa_StartStream(stem);
            if (err_onwk != paNoError) {
                throw std::string(Pa_GetErrorText(err_onwk));
            }
            while(Pa_IsStreamActive(stem) == 1);
            err_onwk = Pa_CloseStream(stem);
            if (err_onwk != paNoError) {
                throw std::string(Pa_GetErrorText(err_onwk));
            }
            std::ofstream my_file;
            my_file.open("out/example_port_audio_1.csv");
            for(unsigned long idx_frm = 0; idx_frm < max_frm; idx_frm++) {
                my_file << usr_dat->recd_ch1[idx_frm] << "," << usr_dat->recd_ch2[idx_frm] << std::endl;
            }
            my_file.close();
        }
        catch (std::exception& e) {
            is_err  = true;
            err_str = e.what();
        }
        catch (std::string& e) {
            is_err  = true;
            err_str = e;
        }
        err_clos = Pa_Terminate();
        if ((is_err == false) && (err_clos != paNoError)) {
            is_err  = true;
            err_str = Pa_GetErrorText(err_clos);
        }
    }
    else {
        is_err  = true;
        err_str = Pa_GetErrorText(err_open);
    }
    if (is_err) {
        FAIL() << err_str;
    }
    else {
        SUCCEED();
    }
}

TEST(TEST_PORT_ADIO, INIT_1) {
    bool        is_err   = false;
    pa_err      err_open;
    pa_err      err_onwk = pa_no_err;
    pa_err      err_clos;
    std::string err_str;
    err_open = open_port();
    if (err_open == pa_no_err) {
        try {
            if (err_onwk != pa_no_err) {
                throw get_err_str(err_onwk);
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
        err_clos = clos_port();
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
        FAIL() << err_str;
    }
    else {
        SUCCEED();
    }
}

TEST(TEST_PORT_ADIO, INIT_2) {
    bool        is_err = false;
    pa_err      err_open;
    pa_err      err_clos;
    std::string err_str;
    err_open = open_port();
    if (err_open == pa_no_err) {
        try {
            port_adio<type::sint32>* obj = new port_adio<type::sint32>(2, 16, 3, 44100);
            obj->open_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            obj->strt_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            while(obj->is_stem_acti());
            obj->clos_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            obj->save("out/");
            delete obj;
        }
        catch (std::exception& e) {
            is_err  = true;
            err_str = e.what();
        }
        catch (std::string& e) {
            is_err  = true;
            err_str = e;
        }
        err_clos = clos_port();
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
        FAIL() << err_str;
    }
    else {
        SUCCEED();
    }
}

TEST(TEST_PORT_ADIO, INIT_3) {
    bool        is_err = false;
    pa_err      err_open;
    pa_err      err_clos;
    std::string err_str;
    err_open = open_port();
    if (err_open == pa_no_err) {
        try {
            port_adio<type::sint32>* obj = new port_adio<type::sint32>("UA-101: USB Audio", 8, 16, 3, 44100);
            obj->open_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            obj->strt_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            while(obj->is_stem_acti());
            obj->clos_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            obj->save("out/");
            delete obj;
        }
        catch (std::exception& e) {
            is_err  = true;
            err_str = e.what();
        }
        catch (std::string& e) {
            is_err  = true;
            err_str = e;
        }
        err_clos = clos_port();
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
        FAIL() << err_str;
    }
    else {
        SUCCEED();
    }
}

TEST(TEST_PORT_ADIO, INIT_4) {
    cap_tim     cap_tim_obj;
    bool        is_err = false;
    pa_err      err_open;
    pa_err      err_clos;
    std::string err_str;
    err_open = open_port();
    if (err_open == pa_no_err) {
        try {
            port_adio<type::sint32>* obj = new port_adio<type::sint32>("UA-101: USB Audio", 8, 16, 3, 44100);
            obj->open_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            obj->strt_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            while(obj->is_stem_acti());
            obj->clos_stem();
            if (obj->is_err()) {
                throw get_err_str(obj->get_err());
            }
            cap_tim_obj.tic();
            obj->move_buf_recd();
            type::mat read = arma::conv_to<type::mat>::from(obj->get_recd());
            read -= obj->get_post_pos_off();
            read /= obj->get_post_pos_amp();
            cap_tim_obj.toc();
            obj->save("out/");
            read.save("out/example_port_audio_2.csv", arma::csv_ascii);
            delete obj;
        }
        catch (std::exception& e) {
            is_err  = true;
            err_str = e.what();
        }
        catch (std::string& e) {
            is_err  = true;
            err_str = e;
        }
        err_clos = clos_port();
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
        FAIL() << err_str;
    }
    else {
        std::cout << GLOG() << "Time " << cap_tim_obj.get_tim() << " usec." << std::endl;
        SUCCEED();
    }
}*/

#endif
