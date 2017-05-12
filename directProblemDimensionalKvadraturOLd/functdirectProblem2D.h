#pragma once

#include <math.h>
#include <complex>

using namespace std;

complex<float> I(0.0f, 1.0f);

float omega = 1.0f;
float c_0 = 1.0f;
float PI = 3.1415926f;

// количество источников
size_t numberSource = 5;
// количество квадратиков по каждому измерению
size_t numberPartitionPosize_ts_N = 50;
// размер квадрата в котором находится неоднородность
float domainInHomogeneity_R = 10.0f;

inline size_t min(size_t a, size_t b)
{
	if (a < b) return a;
	return b;
}

// функция Ханкеля
complex<float> Hankel(float x)
{
	return ((float)_j0(x), (float)_y0(x));
}

// функция Грина
complex<float> G(float x_1, float x_2, float y_1, float y_2)
{
	float dist = sqrt(pow(x_1 - y_1, 2) + pow(x_2 - y_2, 2));
	return -0.25f * I * omega * omega * Hankel(omega * dist / c_0);
}

// задаём 1 источник
complex<float> f_01(float x_1, float x_2)
{
	//физическое местоположение 1 источника
	float q_1 = -1.0f; float q_2 = 0.0f;
	float dist = sqrt(pow(x_1 - q_1, 2) + pow(x_2 - q_2, 2));
	return -0.25f * I * Hankel(omega * dist / c_0);
}
// задаём 2 источник
complex<float> f_02(float x_1, float x_2)
{
	float q_1 = -1.0f; float q_2 = 2.5f;                      //физическое местоположение 1 источника
	float dist = sqrt(pow(x_1 - q_1, 2) + pow(x_2 - q_2, 2));
	return -0.25f * I * Hankel(omega * dist / c_0);
}
// задаём 3 источник
complex<float> f_03(float x_1, float x_2)
{
	float q_1 = -1.0f; float q_2 = 5.0f;                      //физическое местоположение 1 источника
	float dist = sqrt(pow(x_1 - q_1, 2) + pow(x_2 - q_2, 2));
	return -0.25f * I * Hankel(omega * dist / c_0);
}
// задаём 4 источник
complex<float> f_04(float x_1, float x_2)
{
	float q_1 = -1.0f; float q_2 = 7.5f;                      //физическое местоположение 1 источника
	float dist = sqrt(pow(x_1 - q_1, 2) + pow(x_2 - q_2, 2));
	return -0.25f * I * Hankel(omega * dist / c_0);
}
// задаём 5 источник
complex<float> f_05(float x_1, float x_2)
{
	float q_1 = -1.0f; float q_2 = 10.0f;                      //физическое местоположение 1 источника
	float dist = sqrt(pow(x_1 - q_1, 2) + pow(x_2 - q_2, 2));
	return -0.25f * I * Hankel(omega * dist / c_0);
}


// решение системы методом Гаусса
void SolveSlauGaussa(complex<float> **matrix_A, int dimensionMatrix, complex<float> *matrix_b, complex<float> *exactSolution)
{
	complex<double> **basicMatrix, *g_b;
	basicMatrix = new complex<double> *[dimensionMatrix];
	for (int i = 0;i<dimensionMatrix;i++)
	{
		basicMatrix[i] = new complex<double>[dimensionMatrix];
	}
	g_b = new complex<double>[dimensionMatrix];
	for (int i = 0; i < dimensionMatrix; i++)
	{
		for (int j = 0; j< dimensionMatrix; j++)
		{
			basicMatrix[i][j] = static_cast<complex<double>>(matrix_A[i][j]);
		}
	}
	for (int i = 0; i<dimensionMatrix; i++)
	{
		g_b[i] = static_cast<complex<double>>(matrix_b[i]);
	}
	for (int k = 0; k<dimensionMatrix; k++)
	{
		double maxString = abs(basicMatrix[k][k]);
		int maxNomerInString = k;
		for (int i = k + 1;i<dimensionMatrix;i++)
		{
			if (abs(basicMatrix[i][k])>maxString)
			{
				maxString = abs(basicMatrix[i][k]);
				maxNomerInString = i;
			}
		}
		for (int j = 0; j<dimensionMatrix;j++)
		{
			complex<double> c = basicMatrix[k][j];
			basicMatrix[k][j] = basicMatrix[maxNomerInString][j];
			basicMatrix[maxNomerInString][j] = c;
		}
		complex<double> c = g_b[k];
		g_b[k] = g_b[maxNomerInString];
		g_b[maxNomerInString] = c;
		for (int i = k + 1;i<dimensionMatrix;i++)
		{
			complex<double> M = basicMatrix[i][k] / basicMatrix[k][k];
			for (int j = k; j<dimensionMatrix;j++)
			{
				basicMatrix[i][j] -= M*basicMatrix[k][j];
			}
			g_b[i] -= M*g_b[k];
		}
	}
	for (int i = dimensionMatrix - 1; i >= 0; i--)
	{
		complex<double> s = 0;
		for (int j = i + 1; j<dimensionMatrix; j++)
		{
			s += basicMatrix[i][j] * static_cast<complex<double>>(exactSolution[j]);
		}
		exactSolution[i] = static_cast<complex<float>>((g_b[i] - s) / basicMatrix[i][i]);
	}
	for (int i = 0;i<dimensionMatrix;i++)
	{
		delete[] basicMatrix[i];
	}
	delete[] basicMatrix;
	delete[] g_b;
}



/*
void SolveSlauGaussa(complex<float> **matrix_A, size_t dimensionMatrix, complex<float> *matrix_b, complex<float> *exactSolution)
{
	complex<float> **basicMatrix, *g_b;
	basicMatrix = new complex<float> *[dimensionMatrix];
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		basicMatrix[i] = new complex<float>[dimensionMatrix];
	}
	g_b = new complex<float>[dimensionMatrix];


	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		for (size_t j = 0; j < dimensionMatrix; ++j)
		{
			basicMatrix[i][j] = matrix_A[i][j];
		}
	}
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		g_b[i] = matrix_b[i];
	}
	
	float maxString;
	size_t maxNomerInString;
	complex<float> c;
	complex<float> M;
	complex<float> s;

	for (size_t k = 0; k < dimensionMatrix; ++k)
	{
		maxString = abs(basicMatrix[k][k]);
		maxNomerInString = k;

		for (size_t i = k + 1; i < dimensionMatrix; ++i)
		{
			if (abs(basicMatrix[i][k]) > maxString)
			{
				maxString = abs(basicMatrix[i][k]);
				maxNomerInString = i;
			}
		}
		for (size_t j = 0; j < dimensionMatrix; ++j)
		{
			c = basicMatrix[k][j];
			basicMatrix[k][j] = basicMatrix[maxNomerInString][j];
			basicMatrix[maxNomerInString][j] = c;
		}
		c = g_b[k];
		g_b[k] = g_b[maxNomerInString];
		g_b[maxNomerInString] = c;
		for (size_t i = k + 1; i < dimensionMatrix; ++i)
		{
			M = basicMatrix[i][k] / basicMatrix[k][k];
			for (size_t j = k; j < dimensionMatrix; ++j)
			{
				basicMatrix[i][j] -= M * basicMatrix[k][j];
			}
			g_b[i] -= M * g_b[k];
		}
	}

	for (size_t i = dimensionMatrix - 1; i > 0; --i)
	{
		s = { 0.0f, 0.0f };
		for (size_t j = i + 1; j < dimensionMatrix; ++j)
		{
			s += basicMatrix[i][j] * exactSolution[j];
		}
		exactSolution[i] = (g_b[i] - s) / basicMatrix[i][i];
	}
	s = { 0.0f, 0.0f };
	for (size_t j = 1; j < dimensionMatrix; ++j)
	{
		s += basicMatrix[0][j] * exactSolution[j];
	}
	exactSolution[0] = (g_b[0] - s) / basicMatrix[0][0];

	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		delete[] basicMatrix[i];
	}
	delete[] basicMatrix;
	delete[] g_b;
}
*/
// сложение двух квадратных матриц - результат записывается в первую матрицу

void AdditionOfSquareMatrices(size_t dimensionMatrix, complex<float> **matrix_A, complex<float> **matrix_B)
{
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		for (size_t j = 0; j < dimensionMatrix; ++j)
		{
			matrix_A[i][j] = matrix_A[i][j] + matrix_B[i][j];
		}
	}

}


// блочное умножение матриц 
// первый размер матриц -  размер результирующей матрицы
// второй размер 
void MultiplicationMatrixBlock(size_t dimensionMatrix_1, size_t dimensionMatrix_2, size_t dimensionMatrix_3, complex<float> **matrix_A, complex<float> **matrix_B, complex<float> **matrix_C)
{
	for (size_t i = 0; i < dimensionMatrix_1; ++i)
	{
		for (size_t j = 0; j < dimensionMatrix_3; ++j)
		{
			matrix_C[i][j] = 0.0f;
		}
	}
	register size_t sizeBlockParallelMultiplicationMatrix = 32; // размер блока при умножении матриц
	for (size_t jj = 0; jj < dimensionMatrix_3; jj = jj + sizeBlockParallelMultiplicationMatrix)
	{
		for (size_t kk = 0; kk < dimensionMatrix_2; kk = kk + sizeBlockParallelMultiplicationMatrix)
		{
			for (size_t i = 0; i < dimensionMatrix_1; ++i)
			{
				complex<float> *t = matrix_C[i];
				for (size_t j = jj; j < min(jj + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_3); ++j)
				{
					complex<float> r = 0.0f;
					for (size_t k = kk; k < min(kk + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_2); ++k)
					{
						r = r + matrix_A[i][k] * matrix_B[k][j];
					}
					t[j] += r;
				}
			}
		}
	}
}

// блочное умножение транспонированной матрицы на обычную матрицу
void MultiplicationTransposedMatrix(size_t dimensionMatrix_1, size_t dimensionMatrix_2, size_t dimensionMatrix_3, complex<float> **matrix_A, complex<float> **matrix_B, complex<float> **matrix_C)
{
	for (size_t i = 0; i < dimensionMatrix_1; ++i)
	{
		for (size_t j = 0; j < dimensionMatrix_3; ++j)
		{
			matrix_C[i][j] = 0.0f;
		}
	}
	register size_t sizeBlockParallelMultiplicationMatrix = 32; // размер блока при умножении матриц
	for (size_t jj = 0; jj < dimensionMatrix_3; jj = jj + sizeBlockParallelMultiplicationMatrix)
	{
		for (size_t kk = 0; kk < dimensionMatrix_2; kk = kk + sizeBlockParallelMultiplicationMatrix)
		{
			for (size_t i = 0; i < dimensionMatrix_1; ++i)
			{
				complex<float> *t = matrix_C[i];
				for (size_t j = jj; j < min(jj + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_3); ++j)
				{
					complex<float> r = 0.0f;
					for (size_t k = kk; k < min(kk + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_2); ++k)
					{
						r = r + conj(matrix_A[k][i])*matrix_B[k][j];
					}
					t[j] += r;
				}
			}
		}
	}
}


// блочное умножение обычной матрицы на транспонированную матрицу
void MultiplicationMatrixTransposed(size_t dimensionMatrix_1, size_t dimensionMatrix_2, size_t dimensionMatrix_3, complex<float> **matrix_A, complex<float> **matrix_B, complex<float> **matrix_C)
{
	for (size_t i = 0; i < dimensionMatrix_1; ++i)
	{
		for (size_t j = 0; j < dimensionMatrix_3; ++j)
		{
			matrix_C[i][j] = 0.0f;
		}
	}
	register size_t sizeBlockParallelMultiplicationMatrix = 32; // размер блока при умножении матриц
	for (size_t jj = 0; jj < dimensionMatrix_3; jj = jj + sizeBlockParallelMultiplicationMatrix)
	{
		for (size_t kk = 0; kk < dimensionMatrix_2; kk = kk + sizeBlockParallelMultiplicationMatrix)
		{
			for (size_t i = 0; i < dimensionMatrix_1; ++i)
			{
				complex<float> *t = matrix_C[i];
				for (size_t j = jj; j < min(jj + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_3); ++j)
				{
					complex<float> r = 0.0f;
					for (size_t k = kk; k < min(kk + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_2); ++k)
					{
						r = r + matrix_A[i][k] * conj(matrix_B[j][k]);
					}
					t[j] += r;
				}
			}
		}
	}
}


// умножение транспонированной матрицы на вектор
void MultiplicationTransposedMatrixVector(size_t dimensionMatrix_1, size_t dimensionMatrix_2, complex<float> **matrix_A, complex<float> *matrix_b, complex<float> *matrix_c)
{
	for (size_t i = 0; i < dimensionMatrix_1; ++i)
	{
		matrix_c[i] = 0.0f;
	}
	for (size_t i = 0; i < dimensionMatrix_1; ++i)
	{
		for (size_t j = 0; j < dimensionMatrix_2; ++j)
		{
			matrix_c[i] += conj(matrix_A[j][i])*matrix_b[j];
		}
	}
}


// умножение матрицы на вектор
void MultiplicationMatrixVector(size_t dimensionMatrix_1, size_t dimensionMatrix_2, complex<float> **matrix_A, complex<float> *matrix_b, complex<float> *matrix_c)
{
	for (size_t i = 0; i < dimensionMatrix_1; ++i)
	{
		matrix_c[i] = 0.0f;
	}
	for (size_t i = 0; i < dimensionMatrix_1; ++i)
	{
		for (size_t j = 0; j < dimensionMatrix_2; ++j)
		{
			matrix_c[i] += matrix_A[i][j] * matrix_b[j];
		}
	}
}


// умножение вектора на вектор
void MultiplicationVectorVector(size_t dimensionMatrix, complex<float> *matrix_x, complex<float> *matrix_y, complex<float> *matrix_z)
{
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		matrix_z[i] = 0.0f;
	}
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		matrix_z[i] += matrix_x[i] * matrix_y[i];
	}
}


// сложение двух векторов - результат записывается в первый вектор
void AdditionOfSquareVector(size_t dimensionMatrix, complex<float> *matrix_x, complex<float> *matrix_y)
{
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		matrix_x[i] = matrix_x[i] + matrix_y[i];
	}
}


// вычитание двух векторов - результат записывается в первый вектор
void SubtractionOfSquareVector(size_t dimensionMatrix, complex<float> *matrix_x, complex<float> *matrix_y)
{
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		matrix_x[i] = matrix_x[i] - matrix_y[i];
	}
}


// вычитание двух квадратных матриц - результат записывается в первую матрицу
void SubtractionOfSquareMatrices(size_t dimensionMatrix, complex<float> **matrix_A, complex<float> **matrix_B)
{
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		for (size_t j = 0; j < dimensionMatrix; ++j)
		{
			matrix_A[i][j] = matrix_A[i][j] - matrix_B[i][j];
		}
	}

}


//нахождение обратной матрицы методом Жордана-Гаусса
void ComputeInverseMatrixGaussJordan(size_t dimensionMatrix, complex<float> **matrix_A, complex<float> **inverseMatrix_A)
{
	complex<float> **basicMatrix;
	basicMatrix = new complex<float> *[dimensionMatrix];
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		basicMatrix[i] = new complex<float>[dimensionMatrix];
	}
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		for (size_t j = 0; j< dimensionMatrix; ++j)
		{
			basicMatrix[i][j] = matrix_A[i][j];
		}
	}
	complex<float> temp;
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		for (size_t j = 0; j < dimensionMatrix; ++j)
		{
			inverseMatrix_A[i][j] = 0.0f;
		}
		inverseMatrix_A[i][i] = 1.0f;
	}
	for (size_t k = 0; k < dimensionMatrix; ++k)
	{
		temp = basicMatrix[k][k];

		for (size_t j = 0; j < dimensionMatrix; ++j)
		{
			basicMatrix[k][j] /= temp;
			inverseMatrix_A[k][j] /= temp;
		}
		for (size_t i = k + 1; i < dimensionMatrix; ++i)
		{
			temp = basicMatrix[i][k];
			for (size_t j = 0; j < dimensionMatrix; ++j)
			{
				basicMatrix[i][j] -= basicMatrix[k][j] * temp;
				inverseMatrix_A[i][j] -= inverseMatrix_A[k][j] * temp;
			}
		}
	}
	for (size_t k = dimensionMatrix - 1; k > 0; --k)
	{
		for (size_t i = k - 1; i >= 0; --i)
		{
			temp = basicMatrix[i][k];
			for (size_t j = 0; j < dimensionMatrix; ++j)
			{
				basicMatrix[i][j] -= basicMatrix[k][j] * temp;
				inverseMatrix_A[i][j] -= inverseMatrix_A[k][j] * temp;
			}
		}
	}
	for (size_t i = 0; i < dimensionMatrix; ++i)
	{
		delete[] basicMatrix[i];
	}
	delete[] basicMatrix;
}
