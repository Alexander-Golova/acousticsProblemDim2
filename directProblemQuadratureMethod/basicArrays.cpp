#include "stdafx.h"
#include "taskData.h"
#include "basicArrays.h"
#include "basicFunctions.h"

using namespace std;

void GetBasicArrays(vector<vector<vector<vector<float>>>> & a,
	vector<vector<vector<vector<float>>>> & b,
	vector<vector<float>> & c,
	vector<vector<vector<float>>> & overline_a,
	vector<vector<vector<float>>> & overline_b)
{
	// счет индексов метода квадратур
	vector<float> index(NUMBER_PARTITION_POINT + 1);
	for (size_t i = 1; i < NUMBER_PARTITION_POINT; ++i)
	{
		if (i % 2 != 0)
		{
			index[i] = 4.0f / 3;
		}
		else
		{
			index[i] = 2.0f / 3;
		}
	}
	index[0] = 1.0f / 3;
	index[NUMBER_PARTITION_POINT] = 1.0f / 3;

	// нахождение массива a
	float dist;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					if ((i != p) && (q != j))
					{
						dist = step * sqrt(pow(static_cast<float>(i - p), 2) + pow(static_cast<float>(j - q), 2));
						a[i][j][p][q] = index[p] * index[q] * 0.25f * N_0(dist) * pow(OMEGA, 2) * pow(step, 2);
					}
					else
					{
						a[i][j][p][q] = 0.0f;
					}
				}
			}
		}
	}

	// нахождение массива b
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					dist = step * sqrt(pow(static_cast<float>(i - p), 2) + pow(static_cast<float>(j - q), 2));
					b[i][j][p][q] = index[p] * index[q] * (-0.25f) * J_0(dist) * pow(OMEGA, 2) * pow(step, 2);
				}
			}
		}
	}

	// нахождение массива c
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			c[i][j] = 0.0f;
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					if ((i != p) && (q != j))
					{
						dist = step * sqrt(pow(static_cast<float>(i - p), 2) + pow(static_cast<float>(j - q), 2));
						c[i][j] += index[p] * index[q] * 0.25f * N_0(dist) * pow(OMEGA, 2) * pow(step, 2);
					}
					else
					{
						a[i][j][p][q] = 0.0f;
					}
				}
			}
		}
	}

	// нахождение массива overline_a
	for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
	{
		for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
		{
			for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
			{
				dist = sqrt(pow(receiver - static_cast<float>(p) * step, 2) + pow(static_cast<float>(j - q) * step, 2));
				overline_a[j][p][q] = index[p] * index[q] * 0.25f * N_0(dist) * pow(OMEGA, 2) * pow(step, 2);
			}
		}
	}

	// нахождение массива overline_b
	for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
	{
		for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
		{
			for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
			{
				dist = sqrt(pow(receiver - static_cast<float>(p) * step, 2) + pow(static_cast<float>(j - q) * step, 2));
				overline_b[j][p][q] = index[p] * index[q] * (-0.25f) * J_0(dist) * pow(OMEGA, 2) * pow(step, 2);
			}
		}
	}

}