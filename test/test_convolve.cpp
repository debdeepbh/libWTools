#include <iostream>
#include <vector>
#include <string>
#include <WTools.hpp>

using namespace std;

int main()
{

	int N= 4;
	complex<double> A[N] = {{2,0},{1,0},{3,0},{7,0}};
	complex<double> B[N] = {{2,0},{1,0},{3,0},{4,0}};
	complex<double> C[N];
	complex<double> temp[N];
	complex<double> D[N/2];
	complex<double> E[N] = {{2,0},{1,0},{3,0},{4,0},{5,0}};
	complex<double> F[2*N];

	double R[5] = {1, 5, 7, 4, 3};
	complex<double> Rc[5];

	string typestr ("meyer");

	complex<double> mey_u[70];
	complex<double> mey_v[70];
	complex<double> mey_util[70];
	complex<double> mey_vtil[70];


	// compre return 0 if strings are same
	//if (!typestr.compare("meyer"))
	//	cout << "this is true" << endl;
	//
	//WTools::fft(N,A,temp);
	//WTools::cxprint(N, temp);
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

	WTools::filt(10, "d8", mey_u, mey_v, mey_util, mey_vtil);
	WTools::cxprint(10, mey_v);

	return 0;
}
