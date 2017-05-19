#include "stdafx.h"
#include "arrayLoading.h"

using namespace std;

void ArrayLoadingA(vector<vector<vector<vector<complex<double>>>>> & a)
{
	ifstream f_a("matrix_a.txt");
	for (size_t i = 0; i <= NUMBER_PARTITION_POSIZE; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POSIZE; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POSIZE; ++p)
			{
				for (size_t q = 0; q < NUMBER_PARTITION_POSIZE; ++q)
				{
					f_a >> a[i][j][p][q];
				}
			}
		}
	}
	f_a.close();
}

void ArrayLoadingOverlineA(vector<vector<vector<complex<double>>>> & overline_a)
{
	ifstream f_overline_a("matrix_overline_a.txt");
	for (size_t j = 0; j <= NUMBER_PARTITION_POSIZE; ++j)
	{
		for (size_t p = 0; p < NUMBER_PARTITION_POSIZE; ++p)
		{
			for (size_t q = 0; q < NUMBER_PARTITION_POSIZE; ++q)
			{
				f_overline_a >> overline_a[j][p][q];
			}
		}
	}
	f_overline_a.close();
}

void ArrayLoadingB(vector<vector<complex<double>>> & b)
{
	ifstream f_b("matrix_b.txt");
	for (size_t i = 0; i <= NUMBER_PARTITION_POSIZE; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POSIZE; ++j)
		{
			f_b >> b[i][j];
		}
	}
	f_b.close();
}

void ArrayLoadingSource(const size_t numberSource, vector<vector<vector<complex<double>>>> & Source_R,
	vector<vector<complex<double>>>Source_X)
{
	ifstream fileSource("Source.txt");
	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t i = 0; i <= NUMBER_PARTITION_POSIZE; ++i)
		{
			for (size_t j = 0; j <= NUMBER_PARTITION_POSIZE; ++j)
			{
				fileSource >> Source_R[count][i][j];
			}
		}

		for (size_t j = 0; j <= NUMBER_PARTITION_POSIZE; ++j)
		{
			fileSource >> Source_X[count][j];
		}
	}
	fileSource.close();
}

void ArrayLoadingOverlineU(const size_t numberSource, vector<vector<complex<double>>> & overline_u)
{
	ifstream file_overline_u("matrix_overline_u.txt");
	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POSIZE; ++j)
		{
			file_overline_u >> overline_u[count][j];
		}
	}
	file_overline_u.close();
}
