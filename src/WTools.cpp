#include <vector>
#include <complex>
#include <iostream>
#include <cmath>

#include "WTools.hpp"



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
	WTools::fft(N, A, fA);
	WTools::fft(N, B, fB);

	for(int i=0;i<N;i++)
	{
		fC[i] = fA[i]*fB[i];
	}
	WTools::ifft(N, fC, C);	
}

// Print a complex vector
void WTools::cxprint(int N, complex<double>* A)
{
	for(int i=0; i<N; i++)
	{
		std::cout << A[i] << '\t';
	}
		std::cout << std::endl;
}


