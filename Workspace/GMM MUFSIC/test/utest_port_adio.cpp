#include <gtest/gtest.h>
#include <iostream>
#include <port_adio.h>

TEST(TEST_PORT_ADIO, INIT) {
    PaError err;
    err = Pa_Initialize();
    EXPECT_TRUE(err != paNoError);
}
