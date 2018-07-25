#include <vector>
#include <complex>
#include <iostream>
#include <cmath>

#include "WTools.hpp"


// Print a complex vector
void WTools::cxprint(int N, complex<double>* A)
{
	for(int i=0; i<N; i++)
	{
		std::cout << A[i] << '\t';
	}
		std::cout << std::endl;
}

// Basic fft
void WTools::fft(int N, complex<double>* A, complex<double>* B)
{
	fftw_plan plan_a = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(A), reinterpret_cast<fftw_complex*>(B), FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(plan_a);
	fftw_destroy_plan(plan_a);

}

// Basic ifft
void WTools::ifft(int N, complex<double>* A, complex<double>* B)
{
	fftw_plan plan_a = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(A), reinterpret_cast<fftw_complex*>(B), FFTW_BACKWARD, FFTW_ESTIMATE);

	fftw_execute(plan_a);
	fftw_destroy_plan(plan_a);
	for(int i=0;i<N;i++)
	{
		B[i] = B[i]/(complex<double> (N,0));
	}
}

// convolve
//complex<double>* WTools::convolve(int N, complex<double>* A, complex<double>* B)
void WTools::convolve(int N, complex<double>* A, complex<double>* B, complex<double>* C)
{
	complex<double> fA[N];
	complex<double> fB[N];
	complex<double> fC[N];

	// compute the fft of A and B and store in fA and fB
	WTools::fft(N, A, fA);
	WTools::fft(N, B, fB);

	for(int i=0;i<N;i++)
	{
		fC[i] = fA[i]*fB[i];
	}

	// ifft of fC stored in C
	WTools::ifft(N, fC, C);	
}


// Downsample a vector of even length
// Length of B must be N/2, where, length of A is N
void WTools::down(int N, complex<double>* A, complex<double>* B) {
	// check even
	if( N%2 == 1) {
		// error message
		std::cout << "First vector is not of even length" << std::endl;
	} else {
		for(int k=0; k<N/2; k++) {
			B[k] = A[k*2];
		}
	}
}

// Folds a vector A of length N in half and adds it to produce B of length N/2
void WTools::fold(int N, complex<double>* A, complex<double>* B) {
	// check even
	if( N%2 == 1) {
		// error message
		std::cout << "First vector is not of even lenght" << std::endl;
	} else {
		for(int k=0; k<N/2; k++) {
			B[k] = A[k] + A[N/2+k];
		}
	}
}

// recursive implementation of wavelet transform
// % wavelet transform of z wrt parent wavelets u and v with smallest possible dimension sdim
// e.g. for p-th stage wavelets, sdim = N/2^(p-1)
// 



		


