#include "stdafx.h"
#include "initialValue.h"
#include "../directProblemQuadratureMethod/Sources.h"

using namespace std;

void InitialValueU(const size_t numberSource, vector<vector<vector<complex<float>>>> & u,
	vector<vector<vector<complex<float>>>> & Source_R)
{
	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
		{
			for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
			{
				u[count][i][j] = Source_R[count][i][j];
			}
		}
	}
}

void InitialValueXi(vector<vector<complex<float>>> & xi)
{
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			xi[i][j] = static_cast<complex<float>>(0.0f);
		}
	}
}
