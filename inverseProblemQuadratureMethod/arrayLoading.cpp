#include "stdafx.h"
#include "arrayLoading.h"

using namespace std;

void ArrayLoadingA(vector<vector<vector<vector<complex<float>>>>> & a)
{
	ifstream f_a("matrix_a.txt");
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
				{
					f_a >> a[i][j][p][q];
				}
			}
		}
	}
	f_a.close();
}

void ArrayLoadingOverlineA(vector<vector<vector<complex<float>>>> & overline_a)
{
	ifstream f_overline_a("matrix_overline_a.txt");
	for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
	{
		for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
		{
			for (size_t q = 0; q < NUMBER_PARTITION_POINT; ++q)
			{
				f_overline_a >> overline_a[j][p][q];
			}
		}
	}
	f_overline_a.close();
}

void ArrayLoadingB(vector<vector<complex<float>>> & b)
{
	ifstream f_b("matrix_b.txt");
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			f_b >> b[i][j];
		}
	}
	f_b.close();
}

void ArrayLoadingSource(const size_t numberSource, vector<vector<vector<complex<float>>>> & Source_R,
	vector<vector<complex<float>>> & Source_X)
{
	ifstream fileSource("Source.txt");
	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
		{
			for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
			{
				fileSource >> Source_R[count][i][j];
			}
		}

		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			fileSource >> Source_X[count][j];
		}
	}
	fileSource.close();
}

void ArrayLoadingOverlineU(const size_t numberSource, vector<vector<complex<float>>> & overline_u)
{
	ifstream file_overline_u("matrix_overline_u.txt");
	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			file_overline_u >> overline_u[count][j];
		}
	}
	file_overline_u.close();
}

void LoadData(const size_t numberSource, vector<vector<vector<vector<complex<float>>>>> & a,
	vector<vector<vector<complex<float>>>> & overline_a, vector<vector<complex<float>>> & b,
	vector<vector<vector<complex<float>>>> & Source_R, vector<vector<complex<float>>> & Source_X,
	vector<vector<complex<float>>> & overline_u)
{
	ArrayLoadingA(a);
	ArrayLoadingOverlineA(overline_a);
	ArrayLoadingB(b);
	ArrayLoadingSource(numberSource, Source_R, Source_X);
	ArrayLoadingOverlineU(numberSource, overline_u);
}
