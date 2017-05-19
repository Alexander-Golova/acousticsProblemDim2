#include "stdafx.h"
#include "initialValue.h"
#include "../directProblemQuadratureSamplingMethodNew/Sources.h"

using namespace std;

void InitialValueU(const size_t numberSource, vector<vector<vector<complex<double>>>> &u,
	vector<vector<vector<complex<double>>>> Source_R)
{
	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t i = 0; i <= NUMBER_PARTITION_POSIZE; ++i)
		{
			for (size_t j = 0; j <= NUMBER_PARTITION_POSIZE; ++j)
			{
				u[count][i][j] = Source_R[count][i][j];
			}
		}
	}
}

void InitialValueXi(vector<vector<complex<double>>> &xi)
{
	for (size_t i = 0; i <= NUMBER_PARTITION_POSIZE; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POSIZE; ++j)
		{
			xi[i][j] = 0.1;
		}
	}
}

