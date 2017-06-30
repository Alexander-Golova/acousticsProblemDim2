#include "stdafx.h"
#include "taskData.h"
#include "exact_solution.h"

using namespace std;

void GetExactSolution(vector<vector<double>> & xi)
{
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			xi[i][j] = 0.8 * exp(-(i * h - 6.0) * (i * h - 6.0) - (j * h - 6.0) * (j * h - 6.0)) +
				0.2 * exp(-(i * h - 2.0) * (i * h - 2.0) - (j * h - 2.0) * (j * h - 2.0));
		}
	}
}

void WriteSolutionFile(vector<vector<double>> & xi)
{
	ofstream file_xi("exact_xi.txt");
	file_xi << fixed << setprecision(6);
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			file_xi << xi[i][j] << " ";
		}
	}
	file_xi.close();
}
