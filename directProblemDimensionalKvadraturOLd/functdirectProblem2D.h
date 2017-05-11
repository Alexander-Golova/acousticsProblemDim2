#pragma once

#include <math.h>
#include <complex>

using namespace std;

complex<double> I(0.0, 1.0);

double omega = 1.0;
double c_0 = 1.0;
double PI = 3.1415926;

int numberSource = 5;												  // количество источников
int numberPartitionPoints_N = 50;		                          // количество квадратиков по каждому измерению
double domainInHomogeneity_R = 10.0;	                              // размер квадрата в котором находится неоднородность

inline int min(int a, int b)
{
	if (a < b) return a;
	return b;
}

// функция Ханкеля
complex<double> Hankel(double x)
{
	return _j0(x) + I*_y0(x);
}

// функция Грина
complex<double> G(double x_1, double x_2, double y_1, double y_2)
{
	double dist = sqrt(pow(x_1 - y_1, 2) + pow(x_2 - y_2, 2));
	return -0.25*I*omega*omega*Hankel(omega*dist / c_0);
}


// задаём 1 источник
complex<double> f_01(double x_1, double x_2)
{
	double q_1 = -1.0; double q_2 = 0.0;                      //физическое местоположение 1 источника
	double dist = sqrt(pow(x_1 - q_1, 2) + pow(x_2 - q_2, 2));
	return -0.25*I*Hankel(omega*dist / c_0);
}
// задаём 2 источник
complex<double> f_02(double x_1, double x_2)
{
	double q_1 = -1.0; double q_2 = 2.5;                      //физическое местоположение 1 источника
	double dist = sqrt(pow(x_1 - q_1, 2) + pow(x_2 - q_2, 2));
	return -0.25*I*Hankel(omega*dist / c_0);
}
// задаём 3 источник
complex<double> f_03(double x_1, double x_2)
{
	double q_1 = -1.0; double q_2 = 5.0;                      //физическое местоположение 1 источника
	double dist = sqrt(pow(x_1 - q_1, 2) + pow(x_2 - q_2, 2));
	return -0.25*I*Hankel(omega*dist / c_0);
}
// задаём 4 источник
complex<double> f_04(double x_1, double x_2)
{
	double q_1 = -1.0; double q_2 = 7.5;                      //физическое местоположение 1 источника
	double dist = sqrt(pow(x_1 - q_1, 2) + pow(x_2 - q_2, 2));
	return -0.25*I*Hankel(omega*dist / c_0);
}
// задаём 5 источник
complex<double> f_05(double x_1, double x_2)
{
	double q_1 = -1.0; double q_2 = 10.0;                      //физическое местоположение 1 источника
	double dist = sqrt(pow(x_1 - q_1, 2) + pow(x_2 - q_2, 2));
	return -0.25*I*Hankel(omega*dist / c_0);
}


// решение системы методом Гаусса
void SolveSlauGaussa(complex<double> **matrix_A, int dimensionMatrix, complex<double> *matrix_b, complex<double> *exactSolution)
{
	complex<double> **basicMatrix, *g_b;
	basicMatrix = new complex<double> *[dimensionMatrix];
	for (int i = 0;i<dimensionMatrix;i++)
	{
		basicMatrix[i] = new complex<double>[dimensionMatrix];
	}
	g_b = new complex<double>[dimensionMatrix];
	for (int i = 0; i<dimensionMatrix; i++)
	{
		for (int j = 0; j< dimensionMatrix; j++)
		{
			basicMatrix[i][j] = matrix_A[i][j];
		}
	}
	for (int i = 0; i<dimensionMatrix; i++)
	{
		g_b[i] = matrix_b[i];
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
			s += basicMatrix[i][j] * exactSolution[j];
		}
		exactSolution[i] = (g_b[i] - s) / basicMatrix[i][i];
	}
	for (int i = 0;i<dimensionMatrix;i++)
	{
		delete[] basicMatrix[i];
	}
	delete[] basicMatrix;
	delete[] g_b;
}

// сложение двух квадратных матриц - результат записывается в первую матрицу

void AdditionOfSquareMatrices(int dimensionMatrix, complex<double> **matrix_A, complex<double> **matrix_B)
{
	for (int i = 0; i < dimensionMatrix; i++)
	{
		for (int j = 0; j < dimensionMatrix; j++)
		{
			matrix_A[i][j] = matrix_A[i][j] + matrix_B[i][j];
		}
	}

}



// блочное умножение матриц 
// первый размер матриц -  размер результирующей матрицы
// второй размер 
void MultiplicationMatrixBlock(int dimensionMatrix_1, int dimensionMatrix_2, int dimensionMatrix_3, complex<double> **matrix_A, complex<double> **matrix_B, complex<double> **matrix_C)
{
	for (int i = 0; i < dimensionMatrix_1; i++)
	{
		for (int j = 0; j < dimensionMatrix_3; j++)
		{
			matrix_C[i][j] = 0.0;
		}
	}
	register int sizeBlockParallelMultiplicationMatrix = 32; // размер блока при умножении матриц
	for (int jj = 0; jj < dimensionMatrix_3; jj = jj + sizeBlockParallelMultiplicationMatrix)
	{
		for (int kk = 0; kk < dimensionMatrix_2; kk = kk + sizeBlockParallelMultiplicationMatrix)
		{
			for (int i = 0; i < dimensionMatrix_1; i++)
			{
				complex<double> *t = matrix_C[i];
				for (int j = jj; j < min(jj + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_3); j++)
				{
					complex<double> r = 0;
					for (int k = kk; k < min(kk + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_2); k++)
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
void MultiplicationTransposedMatrix(int dimensionMatrix_1, int dimensionMatrix_2, int dimensionMatrix_3, complex<double> **matrix_A, complex<double> **matrix_B, complex<double> **matrix_C)
{
	for (int i = 0; i < dimensionMatrix_1; i++)
	{
		for (int j = 0; j < dimensionMatrix_3; j++)
		{
			matrix_C[i][j] = 0.0;
		}
	}
	register int sizeBlockParallelMultiplicationMatrix = 32; // размер блока при умножении матриц
	for (int jj = 0; jj < dimensionMatrix_3; jj = jj + sizeBlockParallelMultiplicationMatrix)
	{
		for (int kk = 0; kk < dimensionMatrix_2; kk = kk + sizeBlockParallelMultiplicationMatrix)
		{
			for (int i = 0; i < dimensionMatrix_1; i++)
			{
				complex<double> *t = matrix_C[i];
				for (int j = jj; j < min(jj + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_3); j++)
				{
					complex<double> r = 0;
					for (int k = kk; k < min(kk + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_2); k++)
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
void MultiplicationMatrixTransposed(int dimensionMatrix_1, int dimensionMatrix_2, int dimensionMatrix_3, complex<double> **matrix_A, complex<double> **matrix_B, complex<double> **matrix_C)
{
	for (int i = 0; i < dimensionMatrix_1; i++)
	{
		for (int j = 0; j < dimensionMatrix_3; j++)
		{
			matrix_C[i][j] = 0.0;
		}
	}
	register int sizeBlockParallelMultiplicationMatrix = 32; // размер блока при умножении матриц
	for (int jj = 0; jj < dimensionMatrix_3; jj = jj + sizeBlockParallelMultiplicationMatrix)
	{
		for (int kk = 0; kk < dimensionMatrix_2; kk = kk + sizeBlockParallelMultiplicationMatrix)
		{
			for (int i = 0; i < dimensionMatrix_1; i++)
			{
				complex<double> *t = matrix_C[i];
				for (int j = jj; j < min(jj + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_3); j++)
				{
					complex<double> r = 0;
					for (int k = kk; k < min(kk + sizeBlockParallelMultiplicationMatrix, dimensionMatrix_2); k++)
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
void MultiplicationTransposedMatrixVector(int dimensionMatrix_1, int dimensionMatrix_2, complex<double> **matrix_A, complex<double> *matrix_b, complex<double> *matrix_c)
{
	for (int i = 0; i < dimensionMatrix_1; i++)
	{
		matrix_c[i] = 0.0;
	}
	for (int i = 0; i < dimensionMatrix_1; i++)
	{
		for (int j = 0; j < dimensionMatrix_2; j++)
		{
			matrix_c[i] += conj(matrix_A[j][i])*matrix_b[j];
		}
	}
}

// умножение матрицы на вектор
void MultiplicationMatrixVector(int dimensionMatrix_1, int dimensionMatrix_2, complex<double> **matrix_A, complex<double> *matrix_b, complex<double> *matrix_c)
{
	for (int i = 0; i < dimensionMatrix_1; i++)
	{
		matrix_c[i] = 0.0;
	}
	for (int i = 0; i < dimensionMatrix_1; i++)
	{
		for (int j = 0; j < dimensionMatrix_2; j++)
		{
			matrix_c[i] += matrix_A[i][j] * matrix_b[j];
		}
	}
}



// умножение вектора на вектор
void MultiplicationVectorVector(int dimensionMatrix, complex<double> *matrix_x, complex<double> *matrix_y, complex<double> *matrix_z)
{
	for (int i = 0; i < dimensionMatrix; i++)
	{
		matrix_z[i] = 0.0;
	}
	for (int i = 0; i < dimensionMatrix; i++)
	{
		matrix_z[i] += matrix_x[i] * matrix_y[i];
	}
}

// сложение двух векторов - результат записывается в первый вектор

void AdditionOfSquareVector(int dimensionMatrix, complex<double> *matrix_x, complex<double> *matrix_y)
{
	for (int i = 0; i < dimensionMatrix; i++)
	{
		matrix_x[i] = matrix_x[i] + matrix_y[i];
	}
}

// вычитание двух векторов - результат записывается в первый вектор

void SubtractionOfSquareVector(int dimensionMatrix, complex<double> *matrix_x, complex<double> *matrix_y)
{
	for (int i = 0; i < dimensionMatrix; i++)
	{
		matrix_x[i] = matrix_x[i] - matrix_y[i];
	}
}


// вычитание двух квадратных матриц - результат записывается в первую матрицу

void SubtractionOfSquareMatrices(int dimensionMatrix, complex<double> **matrix_A, complex<double> **matrix_B)
{
	for (int i = 0; i < dimensionMatrix; i++)
	{
		for (int j = 0; j < dimensionMatrix; j++)
		{
			matrix_A[i][j] = matrix_A[i][j] - matrix_B[i][j];
		}
	}

}



//нахождение обратной матрицы методом Жордана-Гаусса
void ComputeInverseMatrixGaussJordan(int dimensionMatrix, complex<double> **matrix_A, complex<double> **inverseMatrix_A)
{
	complex<double> **basicMatrix;
	basicMatrix = new complex<double> *[dimensionMatrix];
	for (int i = 0;i<dimensionMatrix;i++)
	{
		basicMatrix[i] = new complex<double>[dimensionMatrix];
	}
	for (int i = 0; i<dimensionMatrix; i++)
	{
		for (int j = 0; j< dimensionMatrix; j++)
		{
			basicMatrix[i][j] = matrix_A[i][j];
		}
	}
	complex<double> temp;
	for (int i = 0; i < dimensionMatrix; i++)
	{
		for (int j = 0; j < dimensionMatrix; j++)
		{
			inverseMatrix_A[i][j] = 0.0;
		}
		inverseMatrix_A[i][i] = 1.0;
	}
	for (int k = 0; k < dimensionMatrix; k++)
	{
		temp = basicMatrix[k][k];

		for (int j = 0; j < dimensionMatrix; j++)
		{
			basicMatrix[k][j] /= temp;
			inverseMatrix_A[k][j] /= temp;
		}
		for (int i = k + 1; i < dimensionMatrix; i++)
		{
			temp = basicMatrix[i][k];
			for (int j = 0; j < dimensionMatrix; j++)
			{
				basicMatrix[i][j] -= basicMatrix[k][j] * temp;
				inverseMatrix_A[i][j] -= inverseMatrix_A[k][j] * temp;
			}
		}
	}
	for (int k = dimensionMatrix - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = basicMatrix[i][k];
			for (int j = 0; j < dimensionMatrix; j++)
			{
				basicMatrix[i][j] -= basicMatrix[k][j] * temp;
				inverseMatrix_A[i][j] -= inverseMatrix_A[k][j] * temp;
			}
		}
	}
	for (int i = 0;i<dimensionMatrix;i++)
	{
		delete[] basicMatrix[i];
	}
	delete[] basicMatrix;
}


