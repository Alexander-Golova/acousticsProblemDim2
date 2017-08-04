#include "stdafx.h"
#include "taskData.h"
#include "exact_solution.h"

using namespace std;

void GetExactSolution(std::vector<std::vector<double>> & xi)
{
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			xi[i][j] = exp(-64.0 * (i * h - 0.6) * (i * h - 0.6) - 64.0 * (j * h - 0.6) * (j * h - 0.6));
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
