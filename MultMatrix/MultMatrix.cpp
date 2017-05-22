#include "stdafx.h"

using namespace std;

const size_t dim = 500;

void PrintMatrix(const vector<vector<complex<double>>> & matrix)
{
	const size_t dim1 = (size_t)matrix.size();
	const size_t dim2 = (size_t)matrix[0].size();

	for (size_t row = 0; row < dim1; ++row)
	{
		for (size_t col = 0; col < dim2; ++col)
		{
			cout << matrix[row][col] << " ";
		}
		cout << endl;
	}
}

void GetNullMatrix(vector<vector<complex<double>>> & matrix)
{
	const size_t dim1 = (size_t)matrix.size();
	const size_t dim2 = (size_t)matrix[0].size();

	for (size_t row = 0; row < dim1; ++row)
	{
		for (size_t col = 0; col < dim2; ++col)
		{
			matrix[row][col] = { 0.0, 0.0 };
		}
	}
}

// Обычное умножение матриц
void MultMatrix(const vector<vector<complex<double>>> & lhs, const vector<vector<complex<double>>> & rhs,
	vector<vector<complex<double>>> & result)
{
	const size_t dim1 = (size_t)lhs.size();
	const size_t dim2 = (size_t)lhs[0].size();
	const size_t dim3 = (size_t)rhs[0].size();

	GetNullMatrix(result);

	for (size_t row = 0; row < dim1; ++row)
	{
		for (size_t col = 0; col < dim3; ++col)
		{
			for (size_t inner = 0; inner < dim2; ++inner)
			{
				result[row][col] += lhs[row][inner] * rhs[inner][col];
			}
		}
	}
}

// умножение матриц, предварительно транспонировав вторую матрицу
void MultTransMatrix(const vector<vector<complex<double>>> & lhs, const vector<vector<complex<double>>> & rhs,
	vector<vector<complex<double>>> & result)
{
	const size_t dim1 = (size_t)lhs.size();
	const size_t dim2 = (size_t)lhs[0].size();
	const size_t dim3 = (size_t)rhs[0].size();

	GetNullMatrix(result);

	vector<complex<double>> thatColumn(dim3);
	vector<complex<double>> thisRow(dim2);
	complex<double> summand;

	for (size_t col = 0; col < dim3; ++col)
	{
		for (size_t inner = 0; inner < dim2; ++inner)
		{
			thatColumn[inner] = rhs[inner][col];
		}
		for (size_t row = 0; row < dim1; ++row)
		{
			thisRow = lhs[row];
			summand = { 0.0, 0.0 };
			for (size_t inner = 0; inner < dim2; ++inner)
			{
				summand += thisRow[inner] * thatColumn[inner];
			}
			result[row][col] = summand;
		}
	}
}


// блочное умножение матриц 
void MultMatrixBlock(const vector<vector<complex<double>>> & lhs, const vector<vector<complex<double>>> & rhs,
	vector<vector<complex<double>>> & result)
{
	const size_t dim1 = (size_t)lhs.size();
	const size_t dim2 = (size_t)lhs[0].size();
	const size_t dim3 = (size_t)rhs[0].size();

	GetNullMatrix(result);

	const register size_t SIZE_BLOCK = 16;
	complex<double> r;

	for (size_t jj = 0; jj < dim3; jj += SIZE_BLOCK)
	{
		for (size_t kk = 0; kk < dim2; kk += SIZE_BLOCK)
		{
			for (size_t row = 0; row < dim1; ++row)
			{
				vector<complex<double>> &ptr = result[row];

				for (size_t col = jj; col < min(jj + SIZE_BLOCK, dim3); ++col)
				{
					r = { 0.0, 0.0 };
					for (size_t inner = kk; inner < min(kk + SIZE_BLOCK, dim2); ++inner)
					{
						r += lhs[row][inner] * rhs[inner][col];
					}
					ptr[col] += r;
				}
			}
		}
	}
}



// возвращает евклидову норму матрицы в квадрате
double GetEuclideanNorm(vector<vector<complex<double>>> & matrix)
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
void SubtractionOfSquareMatrices(vector<vector<complex<double>>> & lhs,
	const vector<vector<complex<double>>> & rhs)
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
	vector<vector<complex<double>>> matrix1(dim, vector<complex<double>>(dim, complex<double>()));
	vector<vector<complex<double>>> matrix2(dim, vector<complex<double>>(dim, complex<double>()));
	vector<vector<complex<double>>> matrix3(dim, vector<complex<double>>(dim, complex<double>()));
	vector<vector<complex<double>>> matrix4(dim, vector<complex<double>>(dim, complex<double>()));

	//matrix1 = { { {1.0f, 0.0f}, {1.0f, 0.0f}, {0.0f, 0.0f} },
	//			{ {0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 0.0f} } };

	//matrix2 = { { {1.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {1.0f, 0.0f} },
	//			{ {1.0f, 0.0f}, {0.0f, 0.0f}, {1.0f, 0.0f}, {2.0f, 0.0f} },
	//			{ {1.0f, 0.0f}, {0.0f, 0.0f}, {2.0f, 0.0f}, {3.0f, 0.0f} } };

	for (size_t row = 0; row < dim; ++row)
	{
		for (size_t col = 0; col < dim; ++col)
		{
			matrix1[row][col] = (1.0 * (rand() / 1000), 1.0 * (rand() / 1000));
			matrix2[row][col] = (1.0 * (rand() / 1000), 1.0 * (rand() / 1000));
		}
	}
	
	// начало счета времени
	double d;
	clock_t timeStart, timeFinish;
	timeStart = clock();

	MultMatrix(matrix1, matrix2, matrix3);

	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "MultMatrix " << d << endl;
	timeStart = clock();

	MultTransMatrix(matrix1, matrix2, matrix3);
	
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "MultTransMatrix " << d << endl;
	timeStart = clock();

	MultMatrixBlock(matrix1, matrix2, matrix4);

	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "MultMatrixBlock " << d << endl;

	//PrintMatrix(matrix3);
	//cout << "" << endl;
	//PrintMatrix(matrix4);
	SubtractionOfSquareMatrices(matrix3, matrix4);
	
	//cout << "" << endl;
	//PrintMatrix(matrix3);
	
	cout << GetEuclideanNorm(matrix3) << endl;

	return 0;
}

