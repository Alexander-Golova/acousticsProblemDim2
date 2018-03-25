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
	// ���� �������� ������ ���������
	vector<float> index(N);
	for (size_t i = 1; i < N - 1; ++i)
	{
		index[i] = 1.0f;
	}
	index[0] = 0.5f;
	index[NUMBER_PARTITION_POINT] = 0.5f; // TODO

	// ���������� ������� a
	float dist;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					dist = step * sqrtf(static_cast<float>((i - p) * (i - p) + (j - q) * (j - q)));
					if (dist > 0.000001f) // DBL_EPSILON
					{
						a[i][j][p][q] = index[p] * index[q] * 0.25f * N_0(dist) * pow(OMEGA, 2) * pow(step, 2); // TODO
					}
					else
					{
						a[i][j][p][q] = 0.0f;
					}
				}
			}
		}
	}

	// ���������� ������� b
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					dist = step * sqrtf(static_cast<float>((i - p) * (i - p) + (j - q) * (j - q)));
					b[i][j][p][q] = index[p] * index[q] * (-0.25f) * J_0(dist) * pow(OMEGA, 2) * pow(step, 2);
				}
			}
		}
	}

	// ���������� ������� c
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					dist = step * sqrtf(static_cast<float>((i - p) * (i - p) + (j - q) * (j - q)));
					if (dist > 0.000001f) // DBL_EPSILON
					{
						c[i][j] += index[p] * index[q] * 0.25f * N_0(dist) * pow(OMEGA, 2) * pow(step, 2);
					}
					else
					{
						c[i][j] = 0.0f;
					}
				}
			}
		}
	}

	// ���������� ������� overline_a
	for (size_t j = 0; j < N; ++j)
	{
		for (size_t p = 0; p < N; ++p)
		{
			for (size_t q = 0; q < N; ++q)
			{
				dist = step * sqrt(static_cast<float>((receiver - p) * (receiver - p) + (j - q) * (j - q)));
				overline_a[j][p][q] = index[p] * index[q] * 0.25f * N_0(dist) * pow(OMEGA, 2) * pow(step, 2);
			}
		}
	}

	// ���������� ������� overline_b
	for (size_t j = 0; j < N; ++j)
	{
		for (size_t p = 0; p < N; ++p)
		{
			for (size_t q = 0; q < N; ++q)
			{
				dist = step * sqrt(static_cast<float>((receiver - p) * (receiver - p) + (j - q) * (j - q)));
				overline_b[j][p][q] = index[p] * index[q] * (-0.25f) * J_0(dist) * pow(OMEGA, 2) * pow(step, 2);
			}
		}
	}
}
