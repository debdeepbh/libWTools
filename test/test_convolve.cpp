#include <iostream>
#include <vector>
#include <WTools.hpp>

using namespace std;

int main()
{

	int N= 4;
	complex<double> A[N] = {{2,0},{1,0},{3,0},{4,0}};
	complex<double> B[N] = {{2,0},{1,0},{3,0},{4,0}};
	complex<double> C[N];
	complex<double> temp[N];
	//
	//WTools::fft(N,A,temp);
	//
	WTools::convolve(N,A,B,C);
	WTools::cxprint(N, C);

	return 0;
}
