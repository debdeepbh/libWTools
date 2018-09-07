#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "WTools.hpp"

// Copies the first N elements of one vector into other 
void WTools::copyVector(int N, complex<double>* source, complex<double>* destination)
{
	for(int i=0; i<N; i++)
		destination[i] = source[i];
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
		std::cout << "First vector is not of even length" << std::endl;
	} else {
		for(int k=0; k<N/2; k++) {
			B[k] = A[k] + A[N/2+k];
		}
	}
}

// recursive implementation of wavelet transform
// % wavelet transform of z wrt parent wavelets u and v with smallest possible dimension sdim
// e.g. for p-th stage wavelets, sdim = length_z/2^(p)
// z and  w must have same length
void WTools::fwt(int z_length, complex<double>* z, int sdim, complex<double>* util, complex<double>* vtil, complex<double>* w)
{
	complex<double> first[z_length/2];
	complex<double> second[z_length/2];



	complex<double> convolved_first[z_length];
	// get the first part
	WTools::convolve(z_length, z, vtil, convolved_first);
	//WTools::cxprint(z_length/2, z);
	WTools::down(z_length, convolved_first, first);


	//debug
	if(z_length == 512)
	{
		WTools::writeReal(z_length/2, first, "first");
		WTools::writeReal(z_length, z, "z_here");
		WTools::writeReal(z_length, vtil, "vtil_here");
		WTools::writeReal(z_length, convolved_first, "vtil_here");
	}


	if( z_length <= 2*sdim )
	{
		std::cout << "Stopping length of vector is " << z_length << std::endl;

		complex<double> convolved[z_length];

		// get the second part
		WTools::convolve(z_length, z, util, convolved);
		WTools::down(z_length, convolved, second);
	}
	else
	{
		std::cout << "recursion, z_length "<< z_length << std::endl;


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

// given a parent wavelet u get the other, v
// works for all the Daubechies wavelets from wiki
// u and v are of the same length N
void WTools::getother(int N, complex<double>* u, complex<double>* v)
{
	for(int k=0; k<N; k++)
	{
		//  ((x % N + N ) % N produces positive number for x%N
		int indexU = ((1-k)% N + N) % N;
		v[k] = pow(-1, k-1)*conj(u[indexU]);
	}
}


// takes a real vector and turns it into a complex vector
void WTools::makeComplex(int N, double* R, complex<double>* C)
{
	for(int i=0; i<N; i++)
	{
		C[i] = complex<double> (R[i]);
	}
}


// get the wavelet filter of correct length
// given type and N, get u, v, util and vtil
// util and vtil are used to compute fwt
// u and v are used to compute ifwt
void WTools::filt(int N, string filterType, complex<double>* output_u, complex<double>* output_v, complex<double>* output_util, complex<double>* output_vtil)
{
	int supportLength;
	// wow, the function compare returns 0 if the strings are equal
	// very counter-intuitive
	if(!filterType.compare("meyer"))
	{

	std::cout << "Using Meyer's wavelets." << std::endl;
	// size = 62
	double u[62] = { -1.009999956941423e-12, 8.519459636796214e-09, -1.111944952595278e-08, -1.0798819539621958e-08, 6.066975741351135e-08, -1.0866516536735883e-07, 8.200680650386481e-08, 1.1783004497663934e-07, -5.506340565252278e-07, 1.1307947017916706e-06, -1.489549216497156e-06, 7.367572885903746e-07, 3.20544191334478e-06, -1.6312699734552807e-05, 6.554305930575149e-05, -0.0006011502343516092, -0.002704672124643725, 0.002202534100911002, 0.006045814097323304, -0.006387718318497156, -0.011061496392513451, 0.015270015130934803, 0.017423434103729693, -0.03213079399021176, -0.024348745906078023, 0.0637390243228016, 0.030655091960824263, -0.13284520043622938, -0.035087555656258346, 0.44459300275757724, 0.7445855923188063, 0.44459300275757724, -0.035087555656258346, -0.13284520043622938, 0.030655091960824263, 0.0637390243228016, -0.024348745906078023, -0.03213079399021176, 0.017423434103729693, 0.015270015130934803, -0.011061496392513451, -0.006387718318497156, 0.006045814097323304, 0.002202534100911002, -0.002704672124643725, -0.0006011502343516092, 6.554305930575149e-05, -1.6312699734552807e-05, 3.20544191334478e-06, 7.367572885903746e-07, -1.489549216497156e-06, 1.1307947017916706e-06, -5.506340565252278e-07, 1.1783004497663934e-07, 8.200680650386481e-08, -1.0866516536735883e-07, 6.066975741351135e-08, -1.0798819539621958e-08, -1.111944952595278e-08, 8.519459636796214e-09, -1.009999956941423e-12, 0.0 };
	// treat it as a complex vector
	WTools::makeComplex(62, u, output_u);

	// fill up the rest with zeros
	for(int i=62; i<N; i++)
		output_u[i] = complex<double> (0);

	// get v
	WTools::getother(N, output_u, output_v);

	} else if (!filterType.compare("d8")){
		supportLength = 8;
		double u[supportLength] ={0.32580343, 1.01094572, 0.89220014, -0.03957503, -0.26450717, 0.0436163, 0.0465036, -0.01498699};
		WTools::makeComplex(supportLength, u, output_u);
	} else if (!filterType.compare("d10")){
		supportLength = 10;
		double u[supportLength] = { 0.22641898, 0.85394354, 1.02432694, 0.19576696, -0.34265671, -0.04560113, 0.10970265, -0.0088268, -0.01779187, 0.00471742793};
		WTools::makeComplex(supportLength, u, output_u);

	} else if (!filterType.compare("d12")){
		supportLength = 12;
		double u[supportLength] = {0.15774243, 0.69950381, 1.06226376, 0.44583132, -0.3199866, -0.18351806, 0.13788809, 0.03892321, -0.04466375, 0.000783251152, 0.00675606236, -0.00152353381};
		WTools::makeComplex(supportLength, u, output_u);
	} else if (!filterType.compare("d14")){
		supportLength = 14;
		double u[supportLength] = {0.11009943, 0.56079128, 1.03114849, 0.66437248, -0.20351382, -0.31683501, 0.1008467, 0.11400345, -0.05378245, -0.02343994, 0.01774979, 0.000607514995, -0.00254790472, 0.000500226853};
		WTools::makeComplex(supportLength, u, output_u);
	} else if (!filterType.compare("d16")){
		supportLength = 16;
		double u[supportLength] = {0.07695562, 0.44246725, 0.95548615, 0.82781653, -0.02238574, -0.40165863, 0.000668194092, 0.18207636, -0.0245639, -0.06235021, 0.01977216, 0.01236884, -0.00688771926, -0.000554004549, 0.000955229711, -0.000166137261};
		WTools::makeComplex(supportLength, u, output_u);

	} else if (!filterType.compare("d18")){
		supportLength = 18;
		double u[supportLength] = {0.05385035, 0.3448343, 0.85534906, 0.92954571, 0.18836955, -0.41475176, -0.13695355, 0.21006834, 0.043452675, -0.09564726, 0.000354892813, 0.03162417, -0.00667962023, -0.00605496058, 0.00261296728, 0.000325814671, -0.000356329759, 5.5645514e-05};
		WTools::makeComplex(supportLength, u, output_u);

	} else if (!filterType.compare("d20")){
		supportLength = 20;
		double u[supportLength] = { 0.03771716, 0.26612218, 0.74557507, 0.97362811, 0.39763774, -0.3533362, -0.27710988, 0.18012745, 0.13160299, -0.10096657, -0.04165925, 0.04696981, 0.00510043697, -0.015179, 0.00197332536, 0.00281768659, -0.00096994784, -0.000164709006, 0.000132354367, -1.875841e-05};
		WTools::makeComplex(supportLength, u, output_u);
	} else
	{
		std::cout<< "Wrong filter type" << std::endl;
	}

	// correction using a scaling factor of sqrt(2) for d8 - d20
	// that were copied form wikipedia
	// also find the other vector v using getother()
	if(!filterType.compare("d8") || !filterType.compare("d10") || !filterType.compare("d12") || !filterType.compare("d14") || !filterType.compare("d16") || !filterType.compare("d18") || !filterType.compare("d20"))
	{
		// copy u into the beginning of output_u
		std::cout << "Support of wavelet basis: " << supportLength << std::endl;

		// pad by zero in the end
		for(int i=supportLength; i<N; i++)
		{
			output_u[i] = complex<double> (0);
		}

		// correction scaling: multiplying by sqrt(2)
		for(int i=0; i<N; i++)
		{
			output_u[i] = output_u[i]/(pow(2,0.5));
		}


		//
		// get output_v
		WTools::getother(N, output_u, output_v);
	}

	// get the tildes util and vtil
	// ////////////////////////////
	complex<double> temp[N];
	// compute fft
	WTools::fft(N, output_u, temp);
	// get the conjugate
	for(int i=0;i<N;i++)
		temp[i] = conj(temp[i]);
	// get the ifft
	WTools::ifft(N, temp, output_util);
	//
	// compute fft
	WTools::fft(N, output_v, temp);
	// get the conjugate
	for(int i=0;i<N;i++)
		temp[i] = conj(temp[i]);
	// get the ifft
	WTools::ifft(N, temp, output_vtil);
}

// a test vector of length 512
void WTools::testvec_gen(complex<double>* testvec)
{
	for(int j=0;j<512; j++)
	{
		testvec[j] = complex<double> ((j-256)*exp(-pow((j-256),2)/512));
	}
}

// write the complex part of a complex array C to ./data/filename
void WTools::writeComplex(int N, complex<double>* C, string filename)
{
	ofstream file_stream("data/" + filename);
	for(int i; i<N; i++)
	{
		// force to show the positive sign of the imaginary part
		file_stream << C[i].real() << showpos << imag(C[i]) << "i" << noshowpos << std::endl;
	}
	file_stream.close();
}

// write the real part of a complex array C to ./data/filename
void WTools::writeReal(int N, complex<double>* C, string filename)
{
	ofstream file_stream("data/" + filename);
	for(int i; i<N; i++)
	{
		file_stream << C[i].real() << std::endl;
	}
	file_stream.close();
}

	// reads data/filename with single column of real values to a complex array C
void WTools::readReal(int N, complex<double>* C, string filename)
{
	ifstream file_stream("data/" + filename);
	double tempReal;
	int i=0;
	while(!file_stream.eof())
	{
		file_stream >> tempReal;
		C[i] = complex<double> ( tempReal );
		i++;
	}
	file_stream.close();
}


// A scaled Wiener deconvolution in the Fourier domain
// We supply the fft of the observed signal (fSignal), fft of the impulse response (fImpulse), the standard deviation of the noise (noiseSD) and a scaling constant (0<scaling<1) and it outputs the fft of deconvolved signal (fOutput) and the Fourier shrinkage multipler (multipler) which is needed for the Wavelet based deconvolution method
// Assumption: the signal, impulse response and the output (deconvolved signal) are of the same length

// TODO: implement vector valued noise signal for noise power
//
void WTools::fWienDec(int N, complex<double>* fSignal, complex<double>* fImpulse, double noiseSD, double scaling, complex<double>* fOutput, complex<double>* multiplier)
{

	complex<double> fOriginal;
	double powerOriginalSq;
	double powerImpulseSq;

	// compute the square of the l^2 norm of fImpulse
	double normSqImpulse = 0;
	for(int j=0;j<N; j++)
	{
		normSqImpulse += pow(abs(fImpulse[j]),2);
	}
	
	
	//
	// store a typical element of naive deconvolution
	complex<double> fNaive;

	for(int i=0; i<N; i++)
	{
		// perform a naive deconvolution
		fNaive = fSignal[i]/fImpulse[i];
		
		// a model for the original signal power to be supplied to the algorithm
		// Typically, fOriginal is same as the observed signal fSignal
		// We supply fOriginal = fSignal / ||fImpulse||_2 
		// i.e. we normalize powOriginalSq by the mean of powImpulseSq
		powerOriginalSq = pow(abs(fSignal[i]),2) / normSqImpulse * N;
	//	powerOriginalSq = pow(abs(fSignal[i]),2);
	//
		powerImpulseSq = pow(abs(fImpulse[i]),2);


		// Compute the Fourier multiplier
		multiplier[i] = powerImpulseSq / ( powerImpulseSq + scaling * N * pow(noiseSD,2) / powerOriginalSq );

		// Perform the scaled Wiener Deconvolution 
		fOutput[i] = fNaive * multiplier[i];

	}
}

// inner product, the first one gets conjugates
// i.e. <a,b> = bar(a).b
complex<double> WTools::innerProduct(int N, complex<double>* A, complex<double>* B)
{
	complex<double> output=0;
	for(int i=0; i<N; i++)
	{
		output +=  conj(A[i]) * B[i];
	}
	return output;
}


// circular shift of a vector at a given index
// shift by zero keeps the vector intact
void WTools::circShift(int N, complex<double>* A, int index, complex<double>* B)
{
	// check validity
	if(index > N)
		std::cout << "index bigger than N" << std::endl;
	else
	{
		for(int i=index; i<N; i++)
			B[i] = A[i-index];
		for(int i=0; i<index; i++)
			B[i] = A[N-index+i];
	}
}
	

// given a p-th wavelet transform wt, rule and a threshold vector of size (p+1), apply thesholds to produce output
void WTools::applyThreshold(int N, complex<double>* wt, string thresholdRule, int p, complex<double>* thresholdVector, complex<double>* output, double* ratioThresholded)
{
	int starting = 0;
	int len = N/2;
	int thresholded, ending, Nj;

	// for the first p finer levels
	for(int k=0; k<p+1; k++)
	{
		thresholded = 0;

		if(k == p)
			ending = N; //for the coarsest level
		else
			ending = starting + len; //all finer levels

		// first, from 0 to N/2-1
		for(int l=starting; l<ending; l++)
		{

			if(!thresholdRule.compare("hard"))
			{
				if(abs(wt[l]) < thresholdVector[k].real())
				{
					output[l] = 0;
					thresholded = thresholded + 1;
				}
				else
					output[l] = wt[l];
			}
			else if(!thresholdRule.compare("soft"))
			{
				if(abs(wt[l]) >= thresholdVector[k].real())
				{
					// works if the wavelet transform is real
					// for complex, multiply by the phase term instead of +1 or -1
					output[l] = (complex<double>((wt[l].real() > 0 ? 1 : -1))) * (abs(wt[l].real()) - thresholdVector[k]);
				}
				else
				{
					output[l] = 0;
					thresholded = thresholded + 1;
				}
			}
			else
			{
				std::cout << "Wrong thresholding rule provided" << std::endl;
			}
		}

		if (k == p)
			Nj =  N - starting + 1;
		else
			Nj = len;

		ratioThresholded[k] = (Nj - thresholded)/Nj;

		// for the next loop
		// lenght related to the next level wavelet coeffs
		len = len/2;
		starting = ending +1;
	}

}

		
// get the basis matrix for p-th stage wavelet tansform
// level = 1, 2, ..., (p+1)
// (p+1)th level is the coarsest level
// u, v are the filters
// Structure of the basisMatrix ///////////
// basisMatrix is N*(p+1) length, but we are storing it as a linear array
// A sample basisMatrix: complex<double> sampleMat[(p+1)*N];
// The last row corresponds to the coarsest wavelet level
// To get the j-th row (j=0,...,p) of this matrix, print
// from (j*N) to (j*N + N) i.e. (for i=0; i<N; i++){ jRow[i] = sampleMat[j*N+i];}
// Use getRow(N, sampleMat, k, output) to get the k-th row
void WTools::getBasisMatrix(int N, complex<double>* u, complex<double>* v, int p, complex<double>* basisMatrix)
{
	complex<double> U[p][N];
       	complex<double> V[p][N];	

	// copy u, v into the first row of U and V
	for(int i=0; i<N; i++)
	{
		U[0][i] = u[i];
		V[0][i] = v[i];
	}

	// for U[k][], ,k=2,..,p

	for(int k=2; k<=(p); k++)
	{

		complex<double> u_folded[N];
		complex<double> v_folded[N];
		// copy u to u_folded here
		WTools::copyVector(N, u, u_folded);
		WTools::copyVector(N, v, v_folded);

		// for the upcoming for loop
		complex<double> uTempVec[N];
		complex<double> vTempVec[N];
		int lengthNow = N;

		// fold it k-1 times
		for(int j=1; j<=(k-1); j++)
		{
			// first lengthNow elements of u_folded are folded into tempVec
			// only lengthNow/2 first elements of tempVec are relevant now
			WTools::fold(lengthNow, u_folded, uTempVec);
			WTools::fold(lengthNow, v_folded, vTempVec);
			// copy the first lengthNow/2 elements to u_folded
			WTools::copyVector(lengthNow/2, uTempVec, u_folded);
			WTools::copyVector(lengthNow/2, vTempVec, v_folded);
			// update lengthNow
			lengthNow = lengthNow/2;
		}


		// Now we need to upsample k-1 times
		complex<double> u_upped[N];
		complex<double> v_upped[N];
		// copy final u_folded to u_upped here
		WTools::copyVector(lengthNow, u_folded, u_upped);
		WTools::copyVector(lengthNow, v_folded, v_upped);

		// k-1 times
		for(int j=1; j<=(k-1); j++)
		{
			// first lengthNow elements of u_folded are upsampled into tempVec
			// only lengthNow first elements of tempVec are relevant now
			WTools::up(lengthNow, u_upped, uTempVec);
			WTools::up(lengthNow, v_upped, vTempVec);
			// copy the first lengthNow*2 relevant elements to u_folded
			WTools::copyVector(lengthNow*2, uTempVec, u_upped);
			WTools::copyVector(lengthNow*2, vTempVec, v_upped);
			// update lengthNow, current length of the vector after upsampling
			lengthNow = lengthNow*2;
		}

		//// good idea to check the final length to see if N
		//std::cout << lengthNow << std::endl;

		// copy folded and upped vector to matrix U, V
		for(int j=0; j<N; j++)
		{
			U[k-1][j] = u_upped[j];
			V[k-1][j] = v_upped[j];
			// last index of U is (p-1) i.e p total
		}
	}

	// Now we compute Psi, the first p rows of the basisMatrix
	complex<double> f_old[N];
	complex<double> g_old[N];

	WTools::copyVector(N, v, f_old);
	WTools::copyVector(N, u, g_old);

	// first row of Psi is f_old
	for(int j=0; j<N; j++)
		basisMatrix[j] = f_old[j];

	// to store jth row of matrix U
	complex<double> ujthRow[N];
	complex<double> vjthRow[N];

	// for the loop
	complex<double> f_new[N];
	complex<double> g_new[N];

	// other rows of Psi
	for(int j=1; j<p; j++)
	{
		for(int l=0; l<N; l++)
		{
			ujthRow[l] = U[j][l];
			vjthRow[l] = V[j][l];
		}

		WTools::convolve(N, g_old, vjthRow, f_new);
		WTools::convolve(N, g_old, ujthRow, g_new);

		// store f_new at j-th row of basisMatrix
		for(int l=0; l<N; l++)
			basisMatrix[j*N+l] = f_new[l];

		// loop update
		WTools::copyVector(N, f_new, f_old);
		WTools::copyVector(N, g_new, g_old);
	}

	// last row of the basisMatrix
	for(int l=0; l<N; l++)
		basisMatrix[p*N+l] = g_old[l];
}

// get the k-th row of a matrix stored in a linear way
// when the row-length of the matrix is N
// and store it in output
// as in getBasisMatrix
void WTools::getRow(int N, complex<double>* matrix, int k, complex<double>* output)
{
	for(int j=0; j<N; j++)
		output[j] = matrix[k*N+j];
}















//// Define sampleM as: int * sampleM[m]; //i.e. sampleM[i] is the pointer to the i-th array
//// then: sampleM[i] = new  int
//void WTools::testMat(int **M, int m, int n)
//{
//	for(int i=0; i<m; i++)
//		for(int j=0; j<m; j++)
//			M[i][j] = i+j;
//}




//// feed the Fourier transform of signal fSignal,
//// fft of impulse response fImpulse
//// for p-th stage wavelet based deconvolution
//// noise standard deviation
//// scaling alpha_j, level dependent, j = 1, 2, ..., p+1
//// wavelet threshold parameter rho
//// thresholdMethod - hard, soft
//// store deconvolved signal in output
//complex<double> WTools::wienForwd(int N, complex<double>* signal, comlex<double>* impulse, string filterType, int p, double noiseSD, double* scaling, double rho, string thresholdMethod, complex<double>* output)
//{
//	// compute the fft of given
//	complex<double> fSignal[N];
//	complex<double> fImpulse[N];
//	fft(N, signal, fSignal);
//	fft(N, impulse, fImpulse);
//
//
//	// at j-th level
//	for(int j=0; j<p+1; j++)
//	{
//		complex<double> fDec[N];
//		complex<double> multiplier[N];
//
//		WTools::fWienDec(N, fSignal, complex<double>* fImpulse, noiseSD, scaling[j], fDec, multiplier);
//
//






	
