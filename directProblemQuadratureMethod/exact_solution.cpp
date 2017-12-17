#include "stdafx.h"
#include "taskData.h"
#include "exact_solution.h"

using namespace std;

void GetExactSolution(vector<vector<double>> & xi)
{
	double sigma = 25;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			xi[i][j] = exp(-((i * h - 0.6) * (i * h - 0.6) + (j * h - 0.6) * (j * h - 0.6)) * sigma);
			//xi[i][j] = 0.8 * exp(-((i * h - 0.8) * (i * h - 0.8) + (j * h - 0.8) * (j * h - 0.8)) * sigma) +
				//0.2 * exp(-((i * h - 0.2) * (i * h - 0.2) + (j * h - 0.2) * (j * h - 0.2)) * sigma);
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
