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

// Upsample a vector 
// Length of B must be 2N, where, length of A is N
void WTools::up(int N, complex<double>* A, complex<double>* B) {
	for(int k=0; k<N; k++) {
		B[2*k] = A[k];
	}
	for(int k=0; k<N; k++) {
		B[2*k+1] = (complex<double> (0,0));
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
// e.g. for p-th stage wavelets, sdim = length_z/2^(p-1)
// z and  w must have same length
void WTools::fwt(int z_length, complex<double>* z, int sdim, complex<double>* util, complex<double>* vtil, complex<double>* w)
{
	complex<double> first[z_length/2];
	complex<double> second[z_length/2];

	complex<double> convolved_first[z_length];
	// get the first part
	WTools::convolve(z_length, z, vtil, convolved_first);
	WTools::down(z_length, convolved_first, first);

	if( z_length <= 2*sdim )
	{
		complex<double> convolved[z_length];

		// get the second part
		WTools::convolve(z_length, z, util, convolved);
		WTools::down(z_length, convolved, second);
	}
	else
	{
		// fold util and vtil and store in util_folded
		// and vtil_folded
		complex<double> util_folded[z_length/2];
		complex<double> vtil_folded[z_length/2];
		WTools::fold(z_length, util, util_folded);
		WTools::fold(z_length, vtil, vtil_folded);

		// first convolve and then downsample
		complex<double> convolved[z_length];
		complex<double> downed[z_length/2];
		WTools::convolve(z_length, z, util, convolved);
		WTools::down(z_length, convolved, downed);

		// get the second part
		// recursion here in the second part
		WTools::fwt(z_length/2, downed, sdim, util_folded, vtil_folded, second);
	}

	// concatenate the first and second part
	for(int i=0; i<z_length/2; i++)
	{
		w[i] = first[i];
	}
	for(int i=z_length/2; i<z_length; i++)
	{
		w[i] = second[i - z_length/2];
	}
}



// recursive implementation of inverse wavelet transform
void WTools::ifwt(int z_length, complex<double>* z, int sdim, complex<double>* u, complex<double>* v, complex<double>* w)
{
	///// implement upsampling
	complex<double> first[z_length];
	complex<double> second[z_length];

	// break z into two parts
	complex<double> z_first[z_length/2];
	complex<double> z_second[z_length/2];
	for (int i=0; i<z_length/2; i++)
		z_first[i] = z[i];
	for (int i=0; i<z_length/2; i++)
		z_second[i] = z[i+ z_length/2];

	// get the first part
	complex<double> upped_first[z_length];
	WTools::up(z_length/2, z_first, upped_first);
	WTools::convolve(z_length, upped_first, v, first);

	if( z_length <= 2*sdim )
	{
		complex<double> upped_second[z_length];
		WTools::up(z_length/2, z_second, upped_second);
		WTools::convolve(z_length, upped_second, u, second);
	}
	else
	{
		// fold u and v and store in u_folded
		// and v_folded
		complex<double> u_folded[z_length/2];
		complex<double> v_folded[z_length/2];
		WTools::fold(z_length, u, u_folded);
		WTools::fold(z_length, v, v_folded);

		// apply recursion
		complex<double> recur_store[z_length/2];
		WTools::ifwt(z_length/2, z_second, sdim, u_folded, v_folded, recur_store);

		// upsample the recursion
		complex<double> upped_recur[z_length];
		WTools::up(z_length/2, recur_store, upped_recur);
		WTools::convolve(z_length, upped_recur, u, second);
	}

	// add the first and second part
	for(int i=0; i<z_length; i++)
	{
		w[i] = first[i] + second[i];
	}
}
