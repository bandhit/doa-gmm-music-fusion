#include <gtest/gtest.h>
#include <port_adio.h>
#include <type_def.hpp>
#include <iostream>
#include <fstream>

TEST(TEST_ORG_PORT_ADIO, INIT_1) {
    bool is_not_err = false;
    PaError err_1;
    PaError err_2;
    err_1 = Pa_Initialize();
    if (err_1 == paNoError) {
        err_2 = Pa_Terminate();
        if (err_2 == paNoError) {
            is_not_err = true;
        }
        else {
            //std::cout << "Terminate: " << Pa_GetErrorText(err_2) << std::endl;
        }
    }
    else {
        //std::cout << "Initialize: " << Pa_GetErrorText(err_1) << std::endl;
    }
    EXPECT_TRUE(is_not_err == true);
}

TEST(TEST_ORG_PORT_ADIO, INIT_2) {
    bool is_not_err = true;
    PaError err_1;
    PaError err_2;
    try {
        err_1 = Pa_Initialize();
        if (err_1 != paNoError) {
            throw "Initialization fail";
        }
    }
    catch (std::exception& e) {
        //std::cout << "Exception: " << e.what() << std::endl;
    }
    try {
        err_2 = Pa_Terminate();
        if (err_2 != paNoError) {
            throw "Termination fail";
        }
    }
    catch (std::exception& e) {
        is_not_err = false;
        //std::cout << "Exception: " << e.what() << std::endl;
    }
    EXPECT_TRUE(is_not_err == true);
}

TEST(TEST_ORG_PORT_ADIO, INIT_3) {
    bool is_not_err = true;
    PaError err_1;
    PaError err_2;
    try {
        err_1 = Pa_Initialize();
        if (err_1 != paNoError) {
            throw "Initialization fail";
        }
        PaDeviceIndex n_dev = Pa_GetDeviceCount();
        for (PaDeviceIndex idx_dev = 0; idx_dev < n_dev; idx_dev++) {
            //std::cout << Pa_GetDeviceInfo(idx_dev)->name << std::endl;
        }
        //std::cout << Pa_GetDeviceInfo(Pa_GetDefaultInputDevice())->name << std::endl;
    }
    catch (std::exception& e) {
        //std::cout << "Exception: " << e.what() << std::endl;
    }
    try {
        err_2 = Pa_Terminate();
        if (err_2 != paNoError) {
            throw "Termination fail";
        }
    }
    catch (std::exception& e) {
        is_not_err = false;
        //std::cout << "Exception: " << e.what() << std::endl;
    }
    EXPECT_TRUE(is_not_err == true);
}

typedef short int type_sam;

typedef struct {
    unsigned long idx_frm;
    unsigned long max_frm;
    type_sam*     recd_ch1;
    type_sam*     recd_ch2;
}
paUserData;

int call_back (const void* in_buf,
               void* out_buf,
               unsigned long frm_per_buf,
               const PaStreamCallbackTimeInfo* tim_inf,
               PaStreamCallbackFlags flag,
               void* usr_dat) {
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
            // *w_ptr_chx++ = *r_ptr++; // add this line, depend on number of channel
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
    bool is_not_err = true;
    PaError err;
    try {
        err = Pa_Initialize();
        if (err != paNoError) {
            throw "Initialization fail";
        }
        PaStreamParameters in_stem_para;
        in_stem_para.device = Pa_GetDefaultInputDevice();
        //in_stem_para.device = 4; // Select device by index
        if (in_stem_para.device == paNoDevice) {
            throw "No default input device";
        }
        in_stem_para.channelCount = 2; // number of channel
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
        if((usr_dat->recd_ch1 == NULL) || (usr_dat->recd_ch2 == NULL)) {
            throw "Could not allocate array";
        }
        for(unsigned long idx_frm = 0; idx_frm < max_frm; idx_frm++) {
            usr_dat->recd_ch1[idx_frm] = 0;
            usr_dat->recd_ch2[idx_frm] = 0;
        }
        PaStream* stem;
        err = Pa_OpenStream(&stem,
                            &in_stem_para,
                            NULL,
                            sam_frq,
                            frm_per_buf,
                            paClipOff,
                            call_back,
                            usr_dat);
        if (err != paNoError) {
            throw "Open stream error";
        }
        err = Pa_StartStream(stem);
        if (err != paNoError) {
            throw "Start stream error";
        }
        while(Pa_IsStreamActive(stem) == 1);
        err = Pa_CloseStream(stem);
        if (err != paNoError) {
            throw "Close stream error";
        }
        /*
        std::ofstream my_file;
        my_file.open("ex_1.csv");
        for(unsigned long idx_frm = 0; idx_frm < max_frm; idx_frm++) {
            my_file << usr_dat->recd_ch1[idx_frm] << "," << usr_dat->recd_ch2[idx_frm] << std::endl;
        }
        my_file.close();
        */
    }
    catch (std::exception& e) {
        is_not_err = false;
        //std::cout << "Exception: " << e.what() << std::endl;
    }
    try {
        err = Pa_Terminate();
        if (err != paNoError) {
            throw "Termination Fail";
        }
    }
    catch (std::exception& e) {
        is_not_err = false;
        //std::cout << "Exception: " << e.what() << std::endl;
    }
    EXPECT_TRUE(is_not_err == true);
}

TEST(TEST_PORT_ADIO, HELLO_PA) {
    bool is_not_err = true;
    try {
        if (!open_port())
            throw "Initialization fail";
        port_adio<type::sint32>* obj = new port_adio<type::sint32>(0, 2, 512, 6, 44100);
        obj->open_stem();
        if (obj->is_err())
            throw obj->get_err();
        obj->strt_stem();
        if (obj->is_err())
            throw obj->get_err();
        while(obj->is_stem_acti());
        obj->clos_stem();
        if (obj->is_err())
            throw obj->get_err();
        obj->save();
    }
    catch (std::exception& e) {
        is_not_err = false;
        //std::cout << "Exception: " << e.what() << std::endl;
    }
    try {
        if (!clos_port())
            throw "Termination Fail";
    }
    catch (std::exception& e) {
        is_not_err = false;
        //std::cout << "Exception: " << e.what() << std::endl;
    }
    EXPECT_TRUE(is_not_err == true);
}
