#include "stdafx.h"
#include "exact_solution.h"

using namespace std;

void ProjectionXi(std::vector<std::vector<std::complex<float>>> & xi)
{
	const size_t N = (size_t)xi.size();
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			if (real(xi[i][j]) <= 0)
			{
				xi[i][j] = { 0.0f, 0.0f };
			}
		}
	}
}

void Prsize_tXi(std::vector<std::vector<std::complex<float>>> & xi, size_t iteration)
{
	const size_t N = (size_t)xi.size();
	ofstream f_xi("approximate_xi_" + to_string(iteration + 1) + ".txt");
	f_xi << fixed << setprecision(6);
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_xi << real(xi[i][j]) << " ";
		}
	}
	f_xi.close();
}

void Renumbering(const vector<vector<complex<float>>> & xi, vector<complex<float>> & numbered_xi)
{
	size_t ii;
	for (size_t i = 0; i <= NUMBER_PARTITION_POsize_tS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POsize_tS; ++j)
		{
			ii = i * (NUMBER_PARTITION_POsize_tS + 1) + j;
			numbered_xi[ii] = xi[i][j];
		}
	}
}

void InverseRenumbering(const vector<complex<float>> & numbered_xi, vector<vector<complex<float>>> & xi)
{
	size_t ii;
	for (size_t i = 0; i <= NUMBER_PARTITION_POsize_tS; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POsize_tS; ++j)
		{
			ii = i * (NUMBER_PARTITION_POsize_tS + 1) + j;
			xi[i][j] = numbered_xi[ii];
		}
	}
}
