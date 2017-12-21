#include "stdafx.h"
#include "taskData.h"
#include "directProblem_utils.h"
#include "Sources.h"

using namespace std;

void GetSubstantiveMatrix(const vector<vector<vector<vector<complex<float>>>>> & a,
	const vector<vector<complex<float>>> & b, const vector<vector<float>> & xi,
	vector<vector<complex<float>>> & substantiveMatrix)
{
	size_t ii, jj;
	complex<float> sumOfTheCoefficients;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			ii = i * (NUMBER_PARTITION_POINT + 1) + j;
			sumOfTheCoefficients = complex<float>();
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					jj = p * (NUMBER_PARTITION_POINT + 1) + q;
					if ((i != p) || (q != j))
					{
						substantiveMatrix[ii][jj] += a[i][j][p][q] * xi[p][q];
						sumOfTheCoefficients += a[i][j][p][q];
					}
				}
			}
			substantiveMatrix[ii][ii] += { 1.0f, 0.0f };
			substantiveMatrix[ii][ii] -= sumOfTheCoefficients * xi[i][j];
			substantiveMatrix[ii][ii] += b[i][j] * xi[i][j];
		}
	}
}

void GetRightPartEquation(const Source & source, size_t count, vector<complex<float>> & rightPartEquation)
{
	size_t ii;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			ii = i * (NUMBER_PARTITION_POINT + 1) + j;
			rightPartEquation[ii] = source.Function(source.node[count], i * h, j * h);
		}
	}
}

void InverseRenumbering(const vector<complex<float>> & numbered_u, vector<vector<complex<float>>> & u)
{
	size_t ii;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			ii = i * (NUMBER_PARTITION_POINT + 1) + j;
			u[i][j] = numbered_u[ii];
		}
	}
}

void GetOverlineU(const Source & source, size_t count,
	const vector<vector<vector<complex<float>>>> & overline_a,
	const vector<vector<float>> & xi, const vector<vector<complex<float>>> & u,
	vector<complex<float>> & overline_u)
{
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		overline_u[i] = source.Function(source.node[count], receiver, i * h);
	}
	for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
	{
		for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
		{
			for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
			{
				overline_u[j] -= overline_a[j][p][q] * xi[p][q] * u[p][q];
			}
		}
	}
}
