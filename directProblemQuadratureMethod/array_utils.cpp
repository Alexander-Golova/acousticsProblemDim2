#include "stdafx.h"

#include "taskData.h"
#include "array_utils.h"

using namespace std;

void WriteBasicArraysFile(vector<vector<vector<vector<float>>>> & a,
	vector<vector<vector<vector<float>>>> & b, vector<vector<float>> & c,
	vector<vector<vector<float>>> & overline_a, vector<vector<vector<float>>> & overline_b)
{
	ofstream f_a("matrix_a.txt");
	f_a << fixed << setprecision(6);
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					f_a << a[i][j][p][q] << " ";
				}
			}
		}
	}
	f_a.close();

	ofstream f_b("matrix_b.txt");
	f_b << fixed << setprecision(6);
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					f_b << b[i][j][p][q] << " ";
				}
			}
		}
	}
	f_b.close();

	ofstream f_c("matrix_c.txt");
	f_c << fixed << setprecision(6);
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			f_c << c[i][j] << " ";
		}
	}
	f_c.close();

	ofstream f_overline_a("matrix_overline_a.txt");
	f_overline_a << fixed << setprecision(6);
	for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
	{
		for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
		{
			for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
			{
				f_overline_a << overline_a[j][p][q] << " ";
			}
		}
	}
	f_overline_a.close();

	ofstream f_overline_b("matrix_overline_b.txt");
	f_overline_b << fixed << setprecision(6);
	for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
	{
		for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
		{
			for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
			{
				f_overline_b << overline_b[j][p][q] << " ";
			}
		}
	}
	f_overline_b.close();
}
