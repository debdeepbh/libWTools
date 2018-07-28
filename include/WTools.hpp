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

	// takes a real vector and turns it into a complex vector
	void makeComplex(int N, double* R, complex<double>* C);

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

	// Upsample a vector by zeros in even positions
	// Length of B must be 2N, where, length of A is N
	auto up(int N, complex<double>* A, complex<double>* B) -> void;

	// folds a vector in half
	// Length of B must be N/2, where, length of A is N
	auto fold(int N, complex<double>* A, complex<double>* B) -> void;

	// fast wavelet transform
	// sdim = z_length/2^p for p-th stage wavelet transform
	auto fwt(int z_length, complex<double>* z, int sdim, complex<double>* util, complex<double>* vtil, complex<double>* w) -> void;

	// fast inverse wavelet transform
	// sdim = z_length/2^p for p-th stage wavelet transform
	auto ifwt(int z_length, complex<double>* z, int sdim, complex<double>* u, complex<double>* v, complex<double>* w) -> void;

	// given one parent wavelet, get the other
	 auto getother(int N, complex<double>* u, complex<double>* v) -> void;
	 // get the parent wavelets from string type
	void filt(int N, string type, complex<double>* output_u, complex<double>* output_v, complex<double>* output_util, complex<double>* output_vtil);
	//
	// a test vector of length 512
	void testvec_gen(complex<double>* testvec);


	// write the real part of a vector to data/
	void writeReal(int N, complex<double>* C, string filename);
    }; 
