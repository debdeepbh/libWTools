#include <iostream>
#include <WTools.hpp>

using namespace std;

int main()
{
	WTools wtools{4};
	fftw_complex A[4] = {{0,0},{1,0},{0,0},{0,0}};
	fftw_complex* C = wtools.fft(A);
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<2; j++)
		{
			cout << C[i][j] << '\n' << endl;
		}
	}
	return 0;
}
