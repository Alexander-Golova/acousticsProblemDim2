#include "stdafx.h"
#include "taskData.h"
#include "exact_solution.h"

using namespace std;

void GetExactSolution(vector<vector<float>> & xi) noexcept
{
	float sigma = 64.0f;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			xi[i][j] = 0.1f * exp(-((i * step - 0.6f) * (i * step - 0.6f) + (j * step - 0.6f) * (j * step - 0.6f)) * sigma);
			//xi[i][j] = 0.2f * exp(-((i * step - 0.8f) * (i * step - 0.8f) + (j * step - 0.8f) * (j * step - 0.8f)) * sigma) +
			//	0.4f * exp(-((i * step - 0.2f) * (i * step - 0.2f) + (j * step - 0.2f) * (j * step - 0.2f)) * sigma);
		}
	}
}

void WriteSolutionFile(vector<vector<float>> & xi) noexcept
{
	ofstream file_xi("exact_xi.txt");
	file_xi << fixed << setprecision(6);
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			file_xi << xi[i][j] << " ";
		}
	}
	file_xi.close();
}
