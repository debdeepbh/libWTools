#include <iostream>
#include <string>
#include <cmath>
#include <WTools.hpp>

using namespace std;

int main()
{

	int N= 4;
	complex<double> A[N] = {{2,1},{1,8},{3,4},{7,0}};
	complex<double> B[N] = {{2,0},{1,0},{3,0},{4,0}};
	complex<double> C[N];
	complex<double> temp[N];
	complex<double> D[N/2];
	complex<double> E[5] = {{2,0},{1,0},{3,0},{4,0},{5,0}};
	complex<double> F[2*N];

	double R[5] = {1, 5, 7, 4, 3};
	complex<double> Rc[5];

	string typestr ("meyer");

	complex<double> mey_u[512];
	complex<double> mey_v[512];
	complex<double> mey_util[512];
	complex<double> mey_vtil[512];
	complex<double> testvec[512];
	complex<double> w[512];	// to store the wavelet transform
	complex<double> z[512];	// to store the inverse wavelet transform
	complex<double> err[512];	// to store the inverse wavelet transform
	complex<double> emptyVec[512];

	// compre return 0 if strings are same
	//if (!typestr.compare("meyer"))
	//	cout << "this is true" << endl;
	//
	//WTools::fft(N,A,temp);
	//WTools::cxprint(N, A);
	//
	//WTools::convolve(N,A,B,C);
	//WTools::cxprint(N, C);
	//
	//WTools::down(N,A,D);
	//WTools::cxprint(N/2, D);
	//
	//WTools::fold(N,A,D);
	//WTools::cxprint(N/2, D);
	//
	//WTools::getother(N, A, C);
	//WTools::cxprint(N, C);
	//
	//WTools::makeComplex(5, R, Rc);
	//WTools::cxprint(5, Rc);

	//WTools::filt(10, "d8", mey_u, mey_v, mey_util, mey_vtil);
	//WTools::cxprint(10, mey_v);
	//
	// get a test vector
	WTools::testvec_gen(testvec);
	WTools::writeReal(512, testvec, "testvec");

	// get the wavelet filters
	WTools::filt(512, "d10", mey_u, mey_v, mey_util, mey_vtil);

	//// debug
	//WTools::writeReal(512, mey_util, "util_test");
	//WTools::writeReal(512, mey_vtil, "vtil_test");

	// perform a forward wavelet transform
	// sdim for 3rd stage wavelet transform
	int sdim = 512/pow(2,3);
	WTools::fwt(512, testvec, sdim, mey_util, mey_vtil, w);
	WTools::writeReal(512, w, "wt");

	// perform an inverse wavelet transform
	WTools::ifwt(512, w, sdim, mey_u, mey_v, z);
	WTools::writeReal(512, z, "check");
	
	// print the error;
	for (int k=0; k<512; k++)
		err[k] =  abs(z[k] - testvec[k]);
	//WTools::writeReal(512, z, "rev");
	//WTools::writeReal(512, err, "err");
	
	//// reading test
	//WTools::readReal(512, emptyVec, "rev");
	//WTools::cxprint(512, emptyVec);
	
	//// circular shift test
	//WTools::circShift(4, B, 1, C);
	//WTools::cxprint(4, B);
	//WTools::cxprint(4, C);
	
	//
	// test applyThreshold on a vector
	//WTools::applyThreshold(int N, complex<double>* wt, string thresholdRule, int p, complex<double>* thresholdVector, complex<double>* output, double* ratioThresholded)
	
	// p = 3
	int ssdim= 512/pow(2,3);
	WTools::fwt(512, testvec, ssdim, mey_util, mey_vtil, w);

	double thresholdVector[4] = {1, 1, 0.1, 10};
	double ratioThresholded[4];
	complex<double> output[512];
	WTools::applyThreshold(512, w, "soft", 3, thresholdVector, output, ratioThresholded);

	WTools::writeReal(512, output, "thresholded");

	complex<double>* Mat[4][5];
	// can fold into a vector or arbitrary length
	//WTools::fold(4, A, Mat[0][]);
	//WTools::testMat(Mat, 4, 5);
	//

	// getBasisMatrix works!
	//int vLength = 512, p=3;
	//complex<double> basisMat[vLength*(p+1)];
	//WTools::getBasisMatrix(vLength, mey_u, mey_v, p, basisMat);
	//for(int i=0; i<p+1; i++)
	//{
	//	WTools::getRow(vLength, basisMat, i, emptyVec);
	//	WTools::writeReal(vLength, emptyVec, "basis" + to_string(i));
	//}
	

//void WTools::fWienDec(int N, complex<double>* fSignal, complex<double>* fImpulse, double noiseSD, double scaling, complex<double>* fOutput, complex<double>* multiplier)
	//int n = 1024;
	//complex<double> signal[n];
	//complex<double> imp[n];
	//complex<double> deconv[n];
	//complex<double> mult[n];
	//complex<double> fsignal[n];
	//complex<double> fimp[n];
	//complex<double> fdeconv[n];
	//complex<double> fdeconvIm[n];
	//WTools::readReal(n, signal, "signal");
	//WTools::readReal(n, imp, "imp");
	//WTools::fft(n, signal, fsignal); 
	//WTools::fft(n, imp, fimp); 
	//WTools::fWienDec(n, fsignal, fimp, 11.2, 1, fdeconv, mult);
	//WTools::writeComplex(n, fdeconv, "fdeconv");
	//for(int i=0; i<n; i++)
	//	fdeconvIm[i] = imag(fdeconv[i]) ;
	//WTools::writeReal(n, fdeconvIm, "fdeconvIm");

	//WTools::ifft(n, fdeconv, deconv);
	//WTools::writeReal(n, deconv, "deconv");
	// scaled wiener deconvolution works!!


		
	return 0;
}
