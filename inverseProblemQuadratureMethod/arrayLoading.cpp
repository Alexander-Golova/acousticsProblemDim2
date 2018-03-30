#include "stdafx.h"
#include "arrayLoading.h"

using namespace std;

void ArrayLoadingA(vector<vector<vector<vector<complex<float>>>>> & a) noexcept
{
	ifstream f_a("matrix_a.txt");
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					f_a >> a[i][j][p][q];
				}
			}
		}
	}
	f_a.close();
}

void ArrayLoadingOverlineA(vector<vector<vector<complex<float>>>> & overline_a) noexcept
{
	ifstream f_overline_a("matrix_overline_a.txt");
	for (size_t j = 0; j < N; ++j)
	{
		for (size_t p = 0; p < N; ++p)
		{
			for (size_t q = 0; q < N; ++q)
			{
				f_overline_a >> overline_a[j][p][q];
			}
		}
	}
	f_overline_a.close();
}

void ArrayLoadingSource(const size_t numberSource, vector<vector<vector<complex<float>>>> & Source_R,
	vector<vector<complex<float>>> & Source_X) noexcept
{
	ifstream fileSource("Source.txt");
	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				fileSource >> Source_R[count][i][j];
			}
		}

		for (size_t j = 0; j < N; ++j)
		{
			fileSource >> Source_X[count][j];
		}
	}
	fileSource.close();
}

void ArrayLoadingOverlineU(const size_t numberSource, vector<vector<complex<float>>> & overline_u) noexcept
{
	ifstream file_overline_u("matrix_overline_u.txt");
	for (size_t count = 0; count < numberSource; ++count)
	{
		for (size_t j = 0; j < N; ++j)
		{
			file_overline_u >> overline_u[count][j];
		}
	}
	file_overline_u.close();
}

void LoadData(const size_t numberSource,
	vector<vector<vector<vector<complex<float>>>>> & a,
	vector<vector<vector<complex<float>>>> & overline_a,
	vector<vector<vector<complex<float>>>> & Source_R,
	vector<vector<complex<float>>> & Source_X,
	vector<vector<complex<float>>> & overline_u) noexcept
{
	ArrayLoadingA(a);
	ArrayLoadingOverlineA(overline_a);

	ArrayLoadingSource(numberSource, Source_R, Source_X);
	ArrayLoadingOverlineU(numberSource, overline_u);
}
