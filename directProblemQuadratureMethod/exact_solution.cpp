#include "stdafx.h"
#include "taskData.h"
#include "exact_solution.h"

using namespace std;

void GetExactSolution(vector<vector<float>> & xi)
{
	float sigma = 25.0f;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			xi[i][j] = exp(-((i * step - 0.6f) * (i * step - 0.6f) + (j * step - 0.6f) * (j * step - 0.6f)) * sigma);
			//xi[i][j] = 0.8 * exp(-((i * h - 0.8) * (i * h - 0.8) + (j * h - 0.8) * (j * h - 0.8)) * sigma) +
				//0.2 * exp(-((i * h - 0.2) * (i * h - 0.2) + (j * h - 0.2) * (j * h - 0.2)) * sigma);
		}
	}
}

void WriteSolutionFile(vector<vector<float>> & xi)
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
