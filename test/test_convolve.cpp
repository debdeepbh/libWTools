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

	complex<double> thresholdVector[4] = {{1,0},{1,0},{0.1,0},{10,0}}; 
	double ratioThresholded[4];
	complex<double> output[512];
	WTools::applyThreshold(512, w, "hard", 3, thresholdVector, output, ratioThresholded);

	WTools::writeReal(512, output, "thresholded");

	complex<double>* Mat[4][5];
	// can fold into a vector or arbitrary length
	//WTools::fold(4, A, Mat[0][]);
	//WTools::testMat(Mat, 4, 5);
	//


		
	return 0;
}
