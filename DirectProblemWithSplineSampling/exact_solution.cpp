#include "stdafx.h"
#include "taskData.h"
#include "exact_solution.h"

using namespace std;

void GetExactSolution(array<array<complex<double>, SPLITTING + 1>, SPLITTING + 1> &xi)
{
	for (size_t i = 0; i <= SPLITTING; ++i)
	{
		for (size_t j = 0; j <= SPLITTING; ++j)
		{
			xi[i][j] = exp(-64.0 * (i * h - 0.6) * (i * h - 0.6) - 64.0 * (j * h - 0.6) * (j * h - 0.6));
		}
	}
}

void WriteSolutionFile(const array<array<complex<double>, SPLITTING + 1>, SPLITTING + 1> & xi)
{
	ofstream file_xi("exact_xi.txt");
	file_xi << fixed << setprecision(6);
	for (size_t i = 0; i <= SPLITTING; ++i)
	{
		for (size_t j = 0; j <= SPLITTING; ++j)
		{
			file_xi << xi[i][j] << " ";
		}
	}
	file_xi.close();
}
