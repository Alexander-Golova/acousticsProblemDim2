#include "stdafx.h"
#include "taskData.h"
#include "basicArrays.h"
#include "basicFunctions.h"

using namespace std;

void GetBasicArrays(vector<vector<vector<vector<complex<double>>>>> & a,
     vector<vector<vector<complex<double>>>> & overline_a, vector<vector<complex<double>>> & b)
{
	// ���� �������� ������ ���������
	vector<double> index(NUMBER_PARTITION_POINT + 1);

	for (size_t i = 1; i < NUMBER_PARTITION_POINT; ++i)
	{
		if (i % 2 != 0)
		{
			index[i] = 4.0 / 3;
		}
		else
		{
			index[i] = 2.0 / 3;
		}
	}
	index[0] = 1.0 / 3;
	index[NUMBER_PARTITION_POINT] = 1.0 / 3;

	// ���������� ������� a
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					if ((i != p) || (q != j))
					{
						a[i][j][p][q] = index[p] * index[q];
						a[i][j][p][q] = a[i][j][p][q] * G(i * h, j * h, p * h, q * h);
						a[i][j][p][q] = a[i][j][p][q] * h * h;
					}
				}
			}
		}
	}

	// ���������� ������� overline_a
	for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
	{
		for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
		{
			for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
			{
				overline_a[j][p][q] = index[p] * index[q];
				overline_a[j][p][q] = overline_a[j][p][q] * G(receiver, j * h, p * h, q * h);
				overline_a[j][p][q] = overline_a[j][p][q] * h * h;
			}
		}
	}

	// ���������� ������� b
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					if (i != p)
					{
						b[i][j] += G(i * h, j * h, p * h, q * h);
					}
				}
				b[i][j] *= h * h;
			}
		}
	}
}