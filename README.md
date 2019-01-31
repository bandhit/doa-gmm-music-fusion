# Sound Angle Estimation by Fusion of Gaussian Mixture Model and Multiple Signal classification

This is a c++ library for one-dimensional direction-of-arrival (sound direction) estimation method.
We employs a Gaussian mixture model with a maximum likelihood estimation algorithm and MUltiple SIgnal Classification (MUSIC) algorithm.
The idea can be found in the following paper; https://www.researchgate.net/publication/318555734_Multiple_Frequency_and_Source_Angle_Estimation_by_Gaussian_Mixture_Model_with_Modified_Microphone_Array_Data_Model.
Note that microphones are arranged in a one straight line.

Preferred IDE: Eclipse IDE for C/C++ Developers 4.5.0

Mathematical matrix operations library: fftw-3.3.6-pl1, lapack-3.7.0, openblas_haswellp-r0.2.19, armadillo-7.600.2

Audio library: portaudio (pa_stable_v190600_20161030)

Testing and mocking framework: gtest

Preprocessor notes:\
-DARMA_DONT_USE_WRAPPER\
-DARMA_DONT_PRINT_ERRORS\
-DUSING_UNITTEST % for run unittest, otherwise please don't use it

Library flag notes:\
-llapack\
-lopenblas\
-lgfortran\
-lpthread\
-lfftw3f_threads\
-lfftw3_threads\
-lfftw3f\
-lfftw3\
-lportaudio

Enjoy!
