#include "stdafx.h"
#include "taskData.h"
#include "basicArrays.h"
#include "basicFunctions.h"

using namespace std;

void GetBasicArrays(vector<vector<vector<vector<complex<float>>>>> & a,
	vector<vector<vector<complex<float>>>> & overline_a) noexcept
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
	Point x_ij, x_pq;
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					x_ij = { step * i, step * j };
					x_pq = { step * p, step * q };
					a[i][j][p][q] = index[p] * index[q] * G(x_ij, x_pq);
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
				x_ij = { receiver, step * j};
				x_pq = { step * p, step * q };
				overline_a[j][p][q] = index[p] * index[q] * G(x_ij, x_pq);
			}
		}
	}
}