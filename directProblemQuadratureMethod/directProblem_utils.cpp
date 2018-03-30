#include "stdafx.h"
#include "taskData.h"
#include "directProblem_utils.h"
#include "Sources.h"

using namespace std;

void GetSubstantiveMatrix(const vector<vector<vector<vector<complex<float>>>>> & a,
	const vector<vector<float>> & xi,
	vector<vector<complex<float>>> & substantiveMatrix) noexcept
{
	size_t ii, jj;
	for (ii = 0; ii < N_SQUARED; ++ii)
	{
		for (jj = 0; jj < N_SQUARED; ++jj)
		{
			substantiveMatrix[ii][jj] = { 0.0f, 0.0f };
		}
	}

	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			ii = i * N + j;
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					jj = p * N + q;
					substantiveMatrix[ii][jj] += a[i][j][p][q] * xi[p][q];
				}
			}
			substantiveMatrix[ii][ii] += 1.0f;
		}
	}
}

void GetRightPartEquation(const Source & source, const size_t count, vector<complex<float>> & rightPartEquation) noexcept
{
	size_t ii;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			ii = i * N + j;
			rightPartEquation[ii] = source.Function(source.node[count], step * i, step * j);
		}
	}
}

void InverseRenumbering(const vector<complex<float>> & numbered_u, vector<vector<complex<float>>> & u) noexcept
{
	size_t ii;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			ii = i * N + j;
			u[i][j] = numbered_u[ii];
		}
	}
}

void GetOverlineU(const Source & source, size_t count, const vector<vector<vector<complex<float>>>> & overline_a,
	const vector<vector<float>> & xi, const vector<vector<complex<float>>> & u,
	vector<complex<float>> & overline_u) noexcept
{
	for (size_t i = 0; i < N; ++i)
	{
		overline_u[i] = source.Function(source.node[count], receiver, i * step);
	}
	for (size_t j = 0; j < N; ++j)
	{
		for (size_t p = 0; p < N; ++p)
		{
			for (size_t q = 0; q < N; ++q)
			{
				overline_u[j] -= overline_a[j][p][q] * xi[p][q] * u[p][q];
			}
		}
	}
}
