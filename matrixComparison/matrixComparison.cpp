#include "stdafx.h"

using namespace std;

// возвращает евклидову норму матрицы в квадрате
double GetEuclideanNorm(vector<vector<double>> & matrix)
{
	const size_t dim1 = (size_t)matrix.size();
	const size_t dim2 = (size_t)matrix[0].size();

	double euclideanNorm = 0.0;

	for (size_t row = 0; row < dim1; ++row)
	{
		for (size_t col = 0; col < dim2; ++col)
		{
			euclideanNorm += norm(matrix[row][col]);
		}
	}
	return euclideanNorm;
}

// Вычитание квадратных матриц lhs = lhs - rhs
void SubtractionOfSquareMatrices(vector<vector<double>> & lhs, const vector<vector<double>> & rhs)
{
	const size_t dim = (size_t)lhs.size();

	for (size_t row = 0; row < dim; ++row)
	{
		for (size_t col = 0; col < dim; ++col)
		{
			lhs[row][col] -= rhs[row][col];
		}
	}
}

int main()
{
	size_t N = 50;
	vector<vector<double>> a(N + 1, vector<double>(N + 1, 0.0));
	vector<vector<double>> b(N + 1, vector<double>(N + 1, 0.0));

	ifstream f_a("a.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_a >> a[i][j];
		}
	}
	f_a.close();
	ifstream f_b("b.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_b >> a[i][j];
		}
	}
	f_b.close();

	SubtractionOfSquareMatrices(a, b);

	cout << GetEuclideanNorm(a) << endl;

    return 0;
}

