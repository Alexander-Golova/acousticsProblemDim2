#include "stdafx.h"

using namespace std;

// возвращает евклидову норму матрицы в квадрате
float GetEuclideanNorm(vector<vector<float>> & matrix)
{
	const size_t dim1 = (size_t)matrix.size();
	const size_t dim2 = (size_t)matrix[0].size();

	float euclideanNorm = 0.0f;

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
void SubtractionOfSquareMatrices(vector<vector<float>> & lhs, const vector<vector<float>> & rhs)
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
	vector<vector<float>> a(N + 1, vector<float>(N + 1, 0.0f));
	vector<vector<float>> b(N + 1, vector<float>(N + 1, 0.0f));

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

