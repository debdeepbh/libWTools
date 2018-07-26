#pragma once

#include <fftw3.h>
#include <complex>

using namespace std;

    namespace WTools{
        ///
        /// \brief Convolve time-domain signals `A` and `B`.
        ///

        ///
        /// \brief A wrapper to perform convolution on raw double points
        ///
        /// This function does not check that the memory pointed to by `A` or `B`
        /// is of the right length. Passing `A` or `B` with the wrong length is
        /// undefined behaviour (UB).
        ///
        /// USER MUST FREE RETURNED MEMORY WHEN DONE.
        ///
        //auto convolve(const double* A, const double* B) -> double*;

        ///
        /// \brief Convolve time-domain signals `A` and `B`.
        ///

        ///
        /// \brief A wrapper to perform convolution on raw double points
        ///
        auto deconvolve(const double* input, const double* response) -> double*;

	/// prints a complex vector
        auto cxprint(int length, complex<double>* input) -> void;

	/// basic fft
        auto fft(int length, complex<double>* input, complex<double>* output) -> void;

	/// basic ifft
        auto ifft(int length, complex<double>* input, complex<double>* output) -> void;

	/// circular convolution
        auto convolve(int length, complex<double>* input1, complex<double>* input2, complex<double>* output) -> void;

	// Downsample a vector of even length
	// Length of B must be N/2, where, length of A is N
	auto down(int N, complex<double>* A, complex<double>* B) -> void;

	// folds a vector in half
	// Length of B must be N/2, where, length of A is N
	auto fold(int N, complex<double>* A, complex<double>* B) -> void;

	auto wrec(int z_length, complex<double>* z, int sdim, complex<double>* util, complex<double>* vtil, complex<double>* w) -> void;


    }; 
