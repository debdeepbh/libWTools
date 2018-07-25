#include <vector>
#include <iostream>
#include <cmath>

#include "WTools.hpp"

// convolve two arrays
auto WTools::convolve(const double* A, const double* B) -> double* {

    ////////////////////////////////////////
    // we now take FFT of A

    // we create storage for the transform of A
    std::vector<complex<double>> Af; Af.resize(this->length);

    // create a copy of A as dft_r2c_1d destroys input
    double* A_copy = reinterpret_cast<double*>(calloc(this->length, sizeof(double)));
    std::copy(A, A+this->length, A_copy);

    // FFTW cache's FFTW plans of the same size so this should be very quick
    // note; this plan performs FFTs in-place
    const fftw_plan plan_A{fftw_plan_dft_r2c_1d(this->length,
                                                A_copy,
                                                reinterpret_cast<fftw_complex*>(&Af[0]),
                                                FFTW_ESTIMATE)};

    // take FFT of A
    fftw_execute(plan_A);

    ////////////////////////////////////////
    // we now take FFT of B

    // we create storage for the transform of A
    std::vector<complex<double>> Bf; Bf.resize(this->length);

    // create a copy of A as dft_r2c_1d destroys input
    double* B_copy = reinterpret_cast<double*>(calloc(this->length, sizeof(double)));
    std::copy(B, B+this->length, B_copy);

    // FFTW cache's FFTW plans of the same size so this should be very quick
    // note; this plan performs FFTs in-place
    const fftw_plan plan_B{fftw_plan_dft_r2c_1d(this->length,
                                                B_copy,
                                                reinterpret_cast<fftw_complex*>(&Bf[0]),
                                                FFTW_ESTIMATE)};

    // take FFT of A
    fftw_execute(plan_B);

    ////////////////////////////////////////
    // we multiply the spectrums in place
    for (int i = 0; i < this->length; i++) {

        // multiply two arrays together
        // FFTW transforms are un-normalized
        // so we must divide by length
        Af[i] = Af[i]*Bf[i]/double(this->length);
    }

    ////////////////////////////////////////
    // we have convolved, we now do an inverse FFT

    // we create storage for the transform of A
    double* output = reinterpret_cast<double*>(calloc(this->length, sizeof(double)));

    // FFTW cache's FFTW plans of the same size so this should be very quick
    // note; this plan performs FFTs in-place
    const fftw_plan plan_inv{fftw_plan_dft_c2r_1d(this->length,
                                                  reinterpret_cast<fftw_complex*>(&Af[0]),
                                                  output,
                                                  FFTW_ESTIMATE)};

    // take inverse FFT
    fftw_execute(plan_inv);

    // destroy the plans
    fftw_destroy_plan(plan_A);
    fftw_destroy_plan(plan_B);
    fftw_destroy_plan(plan_inv);

    // return a pointer to the output
    return output;
}

// Basic fft
void WTools::fft(int N, fftw_complex* A, fftw_complex* B)
{
		fftw_plan plan_a = fftw_plan_dft_1d(N, A, B, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan_a); /* repeat as needed */
		fftw_destroy_plan(plan_a);
    //fftw_free(B);
    //return B;
}

