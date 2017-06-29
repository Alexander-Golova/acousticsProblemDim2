#include "stdafx.h"
#include "exact_solution.h"

using namespace std;

void ProjectionXi(vector<vector<complex<double>>> & xi)
{
	const size_t N = (size_t)xi.size();
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			if (real(xi[i][j]) <= 0)
			{
				xi[i][j] = { 0.0, 0.0 };
			}
		}
	}
}

void PrintXi(std::vector<std::vector<std::complex<double>>> & xi, size_t iteration)
{
	const size_t N = static_cast<size_t>(xi.size());
	ofstream f_xi("approximate_xi_" + to_string(iteration + 1) + ".txt");
	f_xi << fixed << setprecision(6);
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			f_xi << real(xi[i][j]) << " ";
		}
	}
	f_xi.close();
}

void Renumbering(const vector<vector<complex<double>>> & xi, vector<complex<double>> & numbered_xi)
{
	size_t ii;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			ii = i * (NUMBER_PARTITION_POINT + 1) + j;
			numbered_xi[ii] = xi[i][j];
		}
	}
}

void InverseRenumbering(const vector<complex<double>> & numbered_xi, vector<vector<complex<double>>> & xi)
{
	size_t ii;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			ii = i * (NUMBER_PARTITION_POINT + 1) + j;
			xi[i][j] = numbered_xi[ii];
		}
	}
}
