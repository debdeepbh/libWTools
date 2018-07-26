#include <iostream>
#include <vector>
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
	WTools::fold(N,A,D);
	WTools::cxprint(N/2, D);
	//

	return 0;
}
