#include "stdafx.h"
#include "exact_solution.h"

using namespace std;

void ProjectionXi(vector<vector<complex<float>>> & xi)
{
	const size_t N = static_cast<size_t>(xi.size());
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			xi[i][j] = real(xi[i][j]);
			if (real(xi[i][j]) <= 0.0f)
			{
				xi[i][j] = { 0.0f, 0.0f };
			}
		}
	}
}

void PrintXi(vector<vector<complex<float>>> & xi, size_t iteration)
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

void RenumberingXi(const vector<vector<complex<float>>> & xi, vector<complex<float>> & numbered_xi)
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

void RenumberingU(const vector<vector<complex<float>>> & u, vector<complex<float>> & numbered_u)
{
	size_t ii;
	for (size_t i = 0; i < NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j < NUMBER_PARTITION_POINT; ++j)
		{
			ii = i * (NUMBER_PARTITION_POINT + 1) + j;
			numbered_u[ii] = u[i][j];
		}
	}
}

void InverseRenumberingXi(const vector<complex<float>> & numbered_xi, vector<vector<complex<float>>> & xi)
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

void InverseRenumberingU(const vector<complex<float>> & numbered_u, vector<vector<complex<float>>> & u)
{
	size_t ii;
	for (size_t i = 0; i < NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j < NUMBER_PARTITION_POINT; ++j)
		{
			ii = i * (NUMBER_PARTITION_POINT + 1) + j;
			u[i][j] = numbered_u[ii];
		}
	}
}
