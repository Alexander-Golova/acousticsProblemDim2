#include "stdafx.h"
#include "../directProblemQuadratureSamplingMethodNew/basicFunctions.h"
#include "../directProblemQuadratureSamplingMethodNew/Sources.h"
#include "../directProblemQuadratureSamplingMethodNew/taskData.h"
#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"

using namespace std;

int main()
{
	size_t numberOfIterations;
	cout << "Enter the number of iterations ";
	cin >> numberOfIterations;

	float alpha;
	cout << "Enter alpha ";
	cin >> alpha;

	float q;
	cout << "Enter q ";
	cin >> q;

	const Source source;

	const size_t N = NUMBER_PARTITION_POINTS;
	const float h = (float)DOMAIN_IN_HOMOGENEITY / N;

	// выделяем память под основные матрицы
	vector<vector<vector<vector<complex<float>>>>> a(N + 1,
		vector<vector<vector<complex<float>>>>(N + 1, vector<vector<complex<float>>>(N + 1,
			vector<complex<float>>(N + 1, complex<float>()))));

	vector<vector<vector<complex<float>>>> overline_a(N + 1, vector<vector<complex<float>>>(N + 1,
		vector<complex<float>>(N + 1, complex<float>())));

	vector<vector<complex<float>>> b(N + 1, vector<complex<float>>(N + 1, complex<float>()));

	// начало счета времени
	clock_t timeStart, timeFinish, timeBegin;
	timeBegin = clock();
	timeStart = clock();
	// Начало загрузки основных матриц
	ifstream f_a("matrix_a.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
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

	ifstream f_overline_a("matrix_overline_a.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		for (size_t p = 0; p < N; ++p)
		{
			for (size_t q = 0; q < N; ++q)
			{
				f_a >> overline_a[j][p][q];
			}
		}
	}
	f_overline_a.close();

	ifstream f_b("matrix_b.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_b >> b[i][j];
		}
	}
	f_b.close();

	//печатаем время работы
	timeFinish = clock();
	float d;
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Download time major arrays " << d << endl;
	timeStart = clock();

	// выделяем память для значений источников
	vector<vector<vector<complex<float>>>> Source_R(source.numberSource, vector<vector<complex<float>>>(N + 1, vector<complex<float>>(N + 1, complex<float>())));
	vector<vector<complex<float>>> Source_X(source.numberSource, vector<complex<float>>(N + 1, complex<float>()));

	// загружаем значения источников
	ifstream fileSource("Source.txt");
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource >> Source_R[count][i][j];
			}
		}

		for (size_t j = 0; j <= N; ++j)
		{
			fileSource >> Source_X[count][j];
		}
	}
	fileSource.close();

	//печатаем время работы
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Download Time arrays sources " << d << endl;
	timeStart = clock();

	// Выделяем память для поля в приемниках
	vector<vector<complex<float>>> overline_u(source.numberSource, vector<complex<float>>(N + 1, complex<float>()));
	
	// Загрузка акустического поля в приёмнике
	ifstream file_overline_u("matrix_overline_u_1.txt");
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			file_overline_u >> overline_u[count][j];
		}
	}
	file_overline_u.close();

	//печатаем время работы
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Loading time of the acoustic field in the receivers " << d << endl;
	timeStart = clock();

	// вспомогательные переменные
	const size_t N_squared = (N + 1) * (N + 1);

	//выделение памяти под массивы производных  F_1, F_2, ...
	vector<vector<vector<complex<float>>>> F_odd(source.numberSource, vector<vector<complex<float>>>(N_squared, vector<complex<float>>(N_squared, complex<float>())));
	vector<vector<vector<complex<float>>>> F_even(source.numberSource, vector<vector<complex<float>>>(N + 1, vector<complex<float>>(N_squared, complex<float>())));

	//выделение памяти под массивы A и B
	vector<vector<vector<complex<float>>>> A(source.numberSource + 1, vector<vector<complex<float>>>(N + 1, vector<complex<float>>(N_squared, complex<float>())));
	vector<vector<complex<float>>> B(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> inverseMatrixB(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> auxiliaryMatrix(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> secondAuxiliaryMatrix(N_squared, vector<complex<float>>(N_squared, complex<float>()));

	// память для хранения значений основного оператора
	vector<vector<complex<float>>> F_part_odd(source.numberSource, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_part_even(source.numberSource, vector<complex<float>>(N + 1, complex<float>()));

	// память для b_0, b_1,...
	vector<vector<complex<float>>> b_right(source.numberSource + 1, vector<complex<float>>(N_squared, complex<float>()));

	// память для u^(1), u^(2), u^(3)
	vector<vector<vector<complex<float>>>> u(source.numberSource + 1, vector<vector<complex<float>>>(N + 1, vector<complex<float>>(N + 1, complex<float>())));

	// память для xi
	vector<vector<complex<float>>> xi(N + 1, vector<complex<float>>(N + 1, complex<float>()));

	// память для перенумерованных переменных и вспомогательного вектора
	vector<vector<complex<float>>> numbered_u(source.numberSource, vector<complex<float>>(N_squared, complex<float>()));
	vector<complex<float>> numbered_xi(N_squared, complex<float>());
	vector<complex<float>> supportingVector_square(N_squared, complex<float>());
	vector<complex<float>> secondSupportingVector_square(N_squared, complex<float>());
	vector<complex<float>> supportingVector(N_squared, complex<float>());

	//печатаем время работы
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Time allocation " << d << endl;
	timeStart = clock();

	// Начало вычислительной части
	// начальное приближение u
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				u[count][i][j] = Source_R[count][i][j];
			}
		}
	}
	
	// начальное приближение xi
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			xi[i][j] = (0.0f, 0.0f);
		}
	}

	size_t coordinate_x;
	size_t coordinate_y;
	size_t ii, jj;
	complex<float> sumOfTheCoefficients;

	// начало основных итераций
	for (size_t iteration = 0; iteration < numberOfIterations; ++iteration)
	{
		cout << endl;
		cout << "Iteration number " << (iteration + 1) << endl;
		cout << "alpha= " << alpha << endl;
		timeStart = clock();
		//////////////////////////////////////////////////////////////////////////
		//строим левую часть СЛАУ основного метода Ньютона
		//////////////////////////////////////////////////////////////////////////
		//
		// определяем матрицы F_01, F_03, F_05, ..., F_09 Якобиана
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				sumOfTheCoefficients = (0.0f, 0.0f);
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t q = 0; q < N; ++q)
					{
						jj = p * (N + 1) + q;
						if ((i != p) || (q != j))
						{
							sumOfTheCoefficients += a[i][j][p][q];
							for (size_t count = 1; count <= source.numberSource; ++count)
							{
								F_odd[count][ii][jj] = a[i][j][p][q] * u[count - 1][p][q];
							}
							F_odd[0][ii][jj] = a[i][j][p][q] * xi[p][q];
						}
					}
				}
				for (size_t count = 1; count <= source.numberSource; ++count)
				{
					F_odd[count][ii][ii] = u[count - 1][i][j] * (b[i][j] - sumOfTheCoefficients);
				}
				F_odd[0][ii][ii] = (1.0f, 0.0f) + xi[i][j] * (b[i][j] - sumOfTheCoefficients);
			}
		}

		// определяем матрицы F_02, F_04, F_06, ..., F_10 Якобиана
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					jj = p * (N + 1) + q;
					for (size_t count = 1; count <= source.numberSource; ++count)
					{
						F_even[count][j][jj] = overline_a[j][p][q] * u[count - 1][p][q];
					}
					F_even[0][j][jj] = overline_a[j][p][q] * xi[p][q];
				}
			}
		}

		// находим матрицы А и В
		// находим матрицы А
		MultTransposedMatrix(F_odd[1], F_odd[1], A[0]);
		for (size_t count = 2; count <= source.numberSource; ++count)
		{
			MultTransposedMatrix(F_odd[count], F_odd[count], auxiliaryMatrix);
			AddSquareMatrices(A[0], auxiliaryMatrix);
		}
		for (size_t count = 1; count <= source.numberSource; ++count)
		{
			MultTransposedMatrix(F_even[count], F_even[count], auxiliaryMatrix);
			AddSquareMatrices(A[0], auxiliaryMatrix);
		}

		for (size_t count = 1; count <= source.numberSource; ++count)
		{
			MultTransposedMatrix(F_odd[count], F_odd[0], A[count]);
			MultTransposedMatrix(F_even[count], F_even[0], auxiliaryMatrix);
			AddSquareMatrices(A[count], auxiliaryMatrix);
		}

		// находим матрицу В
		MultTransposedMatrix(F_odd[0], F_odd[0], B);
		MultTransposedMatrix(F_even[0], F_even[0], auxiliaryMatrix);
		AddSquareMatrices(B, auxiliaryMatrix);
		
		//добавляем alpha к диагонали
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			A[0][ii][ii] += alpha;
			B[ii][ii] += alpha;
		}

		timeFinish = clock();
		d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
		cout << "Calculation time left side " << d << endl;
		timeStart = clock();
		//////////////////////////////////////////////////////////////////////////
		//строим правую часть СЛАУ основного метода Ньютона
		//////////////////////////////////////////////////////////////////////////
		//
		// находим значения основного оператора F
		//
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				for (size_t count = 0; count < source.numberSource; ++count)
				{
					F_part_odd[count][ii] = u[count][i][j];
				}
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t q = 0; q < N; ++q)
					{
						if ((i != p) || (q != j))
						{
							for (size_t count = 0; count < source.numberSource; ++count)
							{
								F_part_odd[count][ii] += a[i][j][p][q] * (xi[p][q] * u[count][p][q] - xi[i][j] * u[count][i][j]);
							}
						}
					}
				}
				for (size_t count = 0; count < source.numberSource; ++count)
				{
					F_part_odd[count][ii] += b[i][j] * xi[i][j] * u[count][i][j];
					F_part_odd[count][ii] -= Source_R[count][i][j];
				}
			}
		}
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t count = 0; count < source.numberSource; ++count)
			{
				F_part_even[count][j] = overline_u[count][j] - Source_X[count][j];
			}
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					for (size_t count = 0; count < source.numberSource; ++count)
					{
						F_part_even[count][j] += overline_a[j][p][q] * xi[p][q] * u[count][p][q];
					}
				}
			}
		}

		// находим F^{\prime}(x)*x и вычитаем F(x) 
		// перенумерация xi
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				numbered_xi[ii] = xi[i][j];
			}
		}

		for (size_t count = 1; count <= source.numberSource; ++count)
		{
			for (size_t i = 0; i <= N; ++i)
			{
				for (size_t j = 0; j <= N; ++j)
				{
					ii = i * (N + 1) + j;
					numbered_u[count - 1][ii] = u[count - 1][i][j];
				}
			}
			MultMatrixVector(F_odd[count], numbered_xi, supportingVector_square);
			for (size_t ii = 0; ii < N_squared; ++ii)
			{
				F_part_odd[count - 1][ii] = supportingVector_square[ii] - F_part_odd[count - 1][ii];
			}
			MultMatrixVector(F_odd[0], numbered_u[count - 1], supportingVector_square);
			for (size_t ii = 0; ii < N_squared; ++ii)
			{
				F_part_odd[count - 1][ii] += supportingVector_square[ii];
			}
			MultMatrixVector(F_even[count], numbered_xi, supportingVector);
			for (size_t ii = 0; ii <= N; ++ii)
			{
				F_part_even[count - 1][ii] = supportingVector_square[ii] - F_part_even[count - 1][ii];
			}
			MultMatrixVector(F_even[0], numbered_u[count - 1], supportingVector);
			for (size_t ii = 0; ii <= N; ++ii)
			{
				F_part_even[count - 1][ii] += supportingVector_square[ii];
			}
		}

		// находим окончательно правую часть b0, b1,...
		MultTransposedMatrixVector(F_odd[1], F_part_odd[0], b_right[0]);
		MultTransposedMatrixVector(F_even[1], F_part_even[0], supportingVector_square);
		AddVectors(b_right[0], supportingVector_square);
		for (size_t count = 2; count <= source.numberSource; ++count)
		{
			MultTransposedMatrixVector(F_odd[count], F_part_odd[count - 1], supportingVector_square);
			AddVectors(b_right[0], supportingVector_square);
			MultTransposedMatrixVector(F_even[count], F_part_even[count - 1], supportingVector_square);
			AddVectors(b_right[0], supportingVector_square);
		}
		
		for (size_t count = 1; count <= source.numberSource; ++count)
		{
			MultTransposedMatrixVector(F_odd[0], F_part_odd[count - 1], b_right[count]);
			MultTransposedMatrixVector(F_even[0], F_part_even[count - 1], supportingVector_square);
			AddVectors(b_right[count], supportingVector_square);
		}

		timeFinish = clock();
		d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
		cout << "Calculation time right side " << d << endl;
		timeStart = clock();

		// составляем основное уравнение для метода Ньютона

		// Находим xi
		// находим обратную матрицу для B
		InvertMatrix(B, inverseMatrixB);

		//для левой части уравнения с xi все складываем в A_00
		for (size_t count = 1; count <= source.numberSource; ++count)
		{
			MultMatrix(A[count], inverseMatrixB, auxiliaryMatrix);
			MultMatrixTransposed(auxiliaryMatrix, A[count], secondAuxiliaryMatrix);
			SubSquareMatrices(A[0], secondAuxiliaryMatrix);
		}

		//для правой части уравнения с xi все складываем в b0
		for (size_t count = 1; count <= source.numberSource; ++count)
		{
			MultMatrixVector(inverseMatrixB, b_right[count], supportingVector_square);
			MultMatrixVector(A[count], supportingVector_square, secondSupportingVector_square);
			SubVectors(b_right[0], secondSupportingVector_square);
		}
		
		// находим xi
		SolveSlauGaussa(A[0], b_right[0], numbered_xi);

		//находим u^{i}
		for (size_t count = 0; count < source.numberSource; ++count)
		{
			MultTransposedMatrixVector(A[count + 1], numbered_xi, supportingVector_square);
			SubVectors(b_right[count + 1], supportingVector_square);
			MultMatrixVector(inverseMatrixB, b_right[count + 1], numbered_u[count]);
		}

		// изменяем alpha для следующей итерации
		alpha = alpha * q;

		// Обратная перенумерация u^{i}, xi
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			coordinate_x = ii / (N + 1);
			coordinate_y = ii % (N + 1);
			for (size_t count = 0; count < source.numberSource; ++count)
			{
				u[count][coordinate_x][coordinate_y] = numbered_u[count][ii];
			}
			xi[coordinate_x][coordinate_y] = numbered_xi[ii];
		}

		// проекция xi >= 0
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				if (real(xi[i][j]) <= 0)
				{
					xi[i][j] = (0.0f, 0.0f);
				}
			}
		}

		// печать результатов итераций
		string s((char*)&iteration);
		ofstream f_xi("approximate_xi_" + s + ".txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				f_xi << fixed << setprecision(6) << real(xi[i][j]) << " ";
			}
		}
		f_xi.close();
	}

	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Calculation time solutions " << d << endl;

	return 0;
}
