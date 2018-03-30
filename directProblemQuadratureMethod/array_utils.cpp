#include "stdafx.h"

#include "taskData.h"
#include "array_utils.h"

using namespace std;

void WriteBasicArraysFile(vector<vector<vector<vector<complex<float>>>>> & a,
	vector<vector<vector<complex<float>>>> & overline_a) noexcept
{
	ofstream f_a("matrix_a.txt"); 
	f_a << fixed << setprecision(6);
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					f_a << a[i][j][p][q] << " ";
				}
			}
		}
	}
	f_a.close();

	ofstream f_overline_a("matrix_overline_a.txt");
	f_overline_a << fixed << setprecision(6);
	for (size_t j = 0; j < N; ++j)
	{
		for (size_t p = 0; p < N; ++p)
		{
			for (size_t q = 0; q < N; ++q)
			{
				f_overline_a << overline_a[j][p][q] << " ";
			}
		}
	}
	f_overline_a.close();
}
