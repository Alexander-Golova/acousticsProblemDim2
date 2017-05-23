#include "stdafx.h"
#include "../directProblemDimensionalKvadraturOLd/functdirectProblem2D.h"

using namespace std;

int main()
{
	size_t numberOfIterations; //количество итераций

	cout << "Enter the number of iterations ";
	cin >> numberOfIterations;
	++numberOfIterations;

	double alpha;
	cout << "Enter alpha ";
	cin >> alpha;

	double multiplier;
	cout << "Enter q ";
	cin >> multiplier;

	size_t N = NUMBER_PARTITION_POINT;

	//
	// выделяем память под основные матрицы
	//
	complex<double> ****a;
	a = new complex<double> ***[N + 1];
	for (size_t i = 0; i <= N; ++i)
	{
		a[i] = new complex<double>**[N + 1];
		for (size_t j = 0; j <= N; ++j)
		{
			a[i][j] = new complex<double>*[N + 1];
			{
				for (size_t p = 0; p <= N; ++p)
				{
					a[i][j][p] = new complex<double>[N + 1];
				}
			}
		}
	}
	complex<double> **b;
	b = new complex<double> *[N + 1];
	for (size_t i = 0; i <= N; ++i)
	{
		b[i] = new complex<double>[N + 1];
	}
	complex<double> ***overline_a;
	overline_a = new complex<double> **[N + 1];
	for (size_t i = 0; i <= N; ++i)
	{
		overline_a[i] = new complex<double>*[N + 1];
		for (size_t j = 0; j <= N; ++j)
		{
			overline_a[i][j] = new complex<double>[N + 1];
		}
	}

	// Начало загрузки основных матриц
	//
	clock_t timeStart, timeFinish; // для замера времени
	timeStart = clock();
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
	double d;
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Download time major arrays " << d << endl;
	timeStart = clock();
	//
	// память для источников
	//
	complex<double> **f_Source_01, **f_Source_03, **f_Source_05, **f_Source_07, **f_Source_09;
	f_Source_01 = new complex<double> *[N + 1];
	f_Source_03 = new complex<double> *[N + 1];
	f_Source_05 = new complex<double> *[N + 1];
	f_Source_07 = new complex<double> *[N + 1];
	f_Source_09 = new complex<double> *[N + 1];
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		f_Source_01[i] = new complex<double>[N + 1];
		f_Source_03[i] = new complex<double>[N + 1];
		f_Source_05[i] = new complex<double>[N + 1];
		f_Source_07[i] = new complex<double>[N + 1];
		f_Source_09[i] = new complex<double>[N + 1];
	}
	complex<double> *f_Source_02, *f_Source_04, *f_Source_06, *f_Source_08, *f_Source_10;
	f_Source_02 = new complex<double>[N + 1];
	f_Source_04 = new complex<double>[N + 1];
	f_Source_06 = new complex<double>[N + 1];
	f_Source_08 = new complex<double>[N + 1];
	f_Source_10 = new complex<double>[N + 1];
	//
	// загружаем значения источника
	// первый источник
	ifstream file_Source_01("f_Source_01.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			file_Source_01 >> f_Source_01[i][j];
		}
	}
	file_Source_01.close();
	ifstream file_Source_02("f_Source_02.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		file_Source_02 >> f_Source_02[i];
	}
	file_Source_02.close();
	// второй источник
	ifstream file_Source_03("f_Source_03.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			file_Source_03 >> f_Source_03[i][j];
		}
	}
	file_Source_03.close();
	ifstream file_Source_04("f_Source_04.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		file_Source_04 >> f_Source_04[i];
	}
	file_Source_04.close();
	// третий источник
	ifstream file_Source_05("f_Source_05.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			file_Source_05 >> f_Source_05[i][j];
		}
	}
	file_Source_05.close();
	ifstream file_Source_06("f_Source_06.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		file_Source_06 >> f_Source_06[i];
	}
	file_Source_06.close();
	// четвертый источник
	ifstream file_Source_07("f_Source_07.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			file_Source_07 >> f_Source_07[i][j];
		}
	}
	file_Source_07.close();
	ifstream file_Source_08("f_Source_08.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		file_Source_08 >> f_Source_08[i];
	}
	file_Source_08.close();
	// пятый источник
	ifstream file_Source_09("f_Source_09.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			file_Source_09 >> f_Source_09[i][j];
		}
	}
	file_Source_09.close();
	ifstream file_Source_10("f_Source_10.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		file_Source_10 >> f_Source_10[i];
	}
	file_Source_10.close();

	//печатаем время работы
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Download Time arrays sources " << d << endl;
	timeStart = clock();
	//
	// Выделяем память для поля в приемниках
	//
	complex<double> *overline_u_1;
	complex<double> *overline_u_2;
	complex<double> *overline_u_3;
	complex<double> *overline_u_4;
	complex<double> *overline_u_5;
	overline_u_1 = new complex<double>[N + 1];
	overline_u_2 = new complex<double>[N + 1];
	overline_u_3 = new complex<double>[N + 1];
	overline_u_4 = new complex<double>[N + 1];
	overline_u_5 = new complex<double>[N + 1];

	// Загрузка акустического поля в приёмнике
	//
	ifstream file_overline_u_1("matrix_overline_u_1.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		file_overline_u_1 >> overline_u_1[j];
	}
	file_overline_u_1.close();
	ifstream file_overline_u_2("matrix_overline_u_2.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		file_overline_u_2 >> overline_u_2[j];
	}
	file_overline_u_2.close();
	ifstream file_overline_u_3("matrix_overline_u_3.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		file_overline_u_3 >> overline_u_3[j];
	}
	file_overline_u_3.close();
	ifstream file_overline_u_4("matrix_overline_u_4.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		file_overline_u_4 >> overline_u_4[j];
	}
	file_overline_u_4.close();
	ifstream file_overline_u_5("matrix_overline_u_5.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		file_overline_u_5 >> overline_u_5[j];
	}
	file_overline_u_5.close();
	//печатаем время работы
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Loading time of the acoustic field in the receivers " << d << endl;
	timeStart = clock();

	// вспомогательные переменные
	size_t N_squared;
	N_squared = (N + 1) * (N + 1);
	//выделение памяти под массивы производных  F_1, F_2, ...
	complex<double> **F_01, **F_02, **F_03, **F_04, **F_05, **F_06, **F_07, **F_08, **F_09, **F_10, **F_o, **F_oo;
	F_01 = new complex<double> *[N_squared];
	F_03 = new complex<double> *[N_squared];
	F_05 = new complex<double> *[N_squared];
	F_07 = new complex<double> *[N_squared];
	F_09 = new complex<double> *[N_squared];
	F_o = new complex<double> *[N_squared];
	for (size_t i = 0; i < N_squared; ++i)
	{
		F_01[i] = new complex<double>[N_squared];
		F_03[i] = new complex<double>[N_squared];
		F_05[i] = new complex<double>[N_squared];
		F_07[i] = new complex<double>[N_squared];
		F_09[i] = new complex<double>[N_squared];
		F_o[i] = new complex<double>[N_squared];
	}
	F_02 = new complex<double> *[N + 1];
	F_04 = new complex<double> *[N + 1];
	F_06 = new complex<double> *[N + 1];
	F_08 = new complex<double> *[N + 1];
	F_10 = new complex<double> *[N + 1];
	F_oo = new complex<double> *[N + 1];
	for (size_t i = 0; i <= N; ++i)
	{
		F_02[i] = new complex<double>[N_squared];
		F_04[i] = new complex<double>[N_squared];
		F_06[i] = new complex<double>[N_squared];
		F_08[i] = new complex<double>[N_squared];
		F_10[i] = new complex<double>[N_squared];
		F_oo[i] = new complex<double>[N_squared];
	}
	//выделение памяти под массивы A и B
	complex<double> **A_00, **A_01, **A_02, **A_03, **A_04, **A_05, **B, **inverseMatrixB, **auxiliaryMatrix, **secondAuxiliaryMatrix;
	A_00 = new complex<double> *[N_squared];
	A_01 = new complex<double> *[N_squared];
	A_02 = new complex<double> *[N_squared];
	A_03 = new complex<double> *[N_squared];
	A_04 = new complex<double> *[N_squared];
	A_05 = new complex<double> *[N_squared];
	B = new complex<double> *[N_squared];
	auxiliaryMatrix = new complex<double> *[N_squared];
	secondAuxiliaryMatrix = new complex<double> *[N_squared];
	inverseMatrixB = new complex<double> *[N_squared];
	for (size_t i = 0; i < N_squared; ++i)
	{
		A_00[i] = new complex<double>[N_squared];
		A_01[i] = new complex<double>[N_squared];
		A_02[i] = new complex<double>[N_squared];
		A_03[i] = new complex<double>[N_squared];
		A_04[i] = new complex<double>[N_squared];
		A_05[i] = new complex<double>[N_squared];
		B[i] = new complex<double>[N_squared];
		auxiliaryMatrix[i] = new complex<double>[N_squared];
		secondAuxiliaryMatrix[i] = new complex<double>[N_squared];
		inverseMatrixB[i] = new complex<double>[N_squared];
	}
	// память для хранения значений основного оператора
	complex<double> *F_part_01, *F_part_02, *F_part_03, *F_part_04, *F_part_05, *F_part_06, *F_part_07, *F_part_08, *F_part_09, *F_part_10;
	F_part_01 = new complex<double>[N_squared];
	F_part_02 = new complex<double>[N + 1];
	F_part_03 = new complex<double>[N_squared];
	F_part_04 = new complex<double>[N + 1];
	F_part_05 = new complex<double>[N_squared];
	F_part_06 = new complex<double>[N + 1];
	F_part_07 = new complex<double>[N_squared];
	F_part_08 = new complex<double>[N + 1];
	F_part_09 = new complex<double>[N_squared];
	F_part_10 = new complex<double>[N + 1];
	// память для b_0, b_1,... 
	complex<double> *b0, *b1, *b2, *b3, *b4, *b5;
	b0 = new complex<double>[N_squared];
	b1 = new complex<double>[N_squared];
	b2 = new complex<double>[N_squared];
	b3 = new complex<double>[N_squared];
	b4 = new complex<double>[N_squared];
	b5 = new complex<double>[N_squared];
	// память для u^(1), u^(2), u^(3)
	complex<double> **u_1, **u_2, **u_3, **u_4, **u_5;
	u_1 = new complex<double> *[N + 1];
	u_2 = new complex<double> *[N + 1];
	u_3 = new complex<double> *[N + 1];
	u_4 = new complex<double> *[N + 1];
	u_5 = new complex<double> *[N + 1];
	for (size_t i = 0; i <= N; ++i)
	{
		u_1[i] = new complex<double>[N + 1];
		u_2[i] = new complex<double>[N + 1];
		u_3[i] = new complex<double>[N + 1];
		u_4[i] = new complex<double>[N + 1];
		u_5[i] = new complex<double>[N + 1];
	}
	// память для xi
	complex<double> **xi;
	xi = new complex<double> *[N + 1];
	for (size_t i = 0; i <= N; ++i)
	{
		xi[i] = new complex<double>[N + 1];
	}
	// память для перенумерованных переменных и вспомогательного вектора
	complex<double> *numbered_u_1, *numbered_u_2, *numbered_u_3, *numbered_u_4, *numbered_u_5, *numbered_xi, *supportingVector_square, *secondSupportingVector_square, *supportingVector;
	numbered_u_1 = new complex<double>[N_squared];
	numbered_u_2 = new complex<double>[N_squared];
	numbered_u_3 = new complex<double>[N_squared];
	numbered_u_4 = new complex<double>[N_squared];
	numbered_u_5 = new complex<double>[N_squared];
	numbered_xi = new complex<double>[N_squared];
	supportingVector_square = new complex<double>[N_squared];
	secondSupportingVector_square = new complex<double>[N_squared];
	supportingVector = new complex<double>[N + 1];

	// Начало вычислительной части
	//
	// начальное приближение u
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			u_1[i][j] = f_Source_01[i][j];
			u_2[i][j] = f_Source_03[i][j];
			u_3[i][j] = f_Source_05[i][j];
			u_3[i][j] = f_Source_07[i][j];
			u_3[i][j] = f_Source_09[i][j];
		}
	}
	// начальное приближение xi
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			xi[i][j] = 0.1;
		}
	}

	// начало основных итераций
	//
	for (size_t iteration = 0; iteration < numberOfIterations; ++iteration)
	{
		cout << endl;
		cout << "Iteration number " << (iteration + 1) << endl;
		cout << "alpha= " << alpha << endl;
		timeStart = clock();
		//////////////////////////////////////////////////////////////////////////
		//строим левую часть СЛАУ основаного метода Ньютона
		//////////////////////////////////////////////////////////////////////////
		size_t ii, jj;
		//
		// определяем матрицы F_01, F_03, F_05, ..., F_09 Якобиана
		//
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				complex<double> sumOfTheCoefficients = 0.0;
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t q = 0; q < N; ++q)
					{
						jj = p * (N + 1) + q;
						if ((i != p) || (q != j))
						{
							sumOfTheCoefficients += a[i][j][p][q];
							F_01[ii][jj] = a[i][j][p][q] * u_1[p][q];
							F_03[ii][jj] = a[i][j][p][q] * u_2[p][q];
							F_05[ii][jj] = a[i][j][p][q] * u_3[p][q];
							F_07[ii][jj] = a[i][j][p][q] * u_4[p][q];
							F_09[ii][jj] = a[i][j][p][q] * u_5[p][q];
							F_o[ii][jj] = a[i][j][p][q] * xi[p][q];
						}
					}
				}
				F_01[ii][ii] = u_1[i][j] * (b[i][j] - sumOfTheCoefficients);
				F_03[ii][ii] = u_2[i][j] * (b[i][j] - sumOfTheCoefficients);
				F_05[ii][ii] = u_3[i][j] * (b[i][j] - sumOfTheCoefficients);
				F_07[ii][ii] = u_4[i][j] * (b[i][j] - sumOfTheCoefficients);
				F_09[ii][ii] = u_5[i][j] * (b[i][j] - sumOfTheCoefficients);
				F_o[ii][ii] = 1.0 + xi[i][j] * (b[i][j] - sumOfTheCoefficients);
			}
		}
		//
		// определяем матрицы F_02, F_04, F_06, ..., F_10 Якобиана
		//
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					jj = p * (N + 1) + q;
					F_02[j][jj] = overline_a[j][p][q] * u_1[p][q];
					F_04[j][jj] = overline_a[j][p][q] * u_2[p][q];
					F_06[j][jj] = overline_a[j][p][q] * u_3[p][q];
					F_08[j][jj] = overline_a[j][p][q] * u_4[p][q];
					F_10[j][jj] = overline_a[j][p][q] * u_5[p][q];
					F_oo[j][jj] = overline_a[j][p][q] * xi[p][q];
				}
			}
		}
		//
		// находим матрицы А и В
		//
		// находим матрицы А_0
		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_01, F_01, A_00);
		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_02, F_02, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_00, auxiliaryMatrix);

		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_03, F_03, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_00, auxiliaryMatrix);

		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_04, F_04, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_00, auxiliaryMatrix);

		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_05, F_05, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_00, auxiliaryMatrix);

		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_06, F_06, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_00, auxiliaryMatrix);

		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_07, F_07, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_00, auxiliaryMatrix);

		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_08, F_08, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_00, auxiliaryMatrix);

		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_09, F_09, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_00, auxiliaryMatrix);

		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_10, F_10, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_00, auxiliaryMatrix);

		// находим матрицу А_1
		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_01, F_o, A_01);
		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_02, F_oo, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_01, auxiliaryMatrix);

		// находим матрицу А_2
		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_03, F_o, A_02);
		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_04, F_oo, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_05, auxiliaryMatrix);

		// находим матрицу А_3
		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_05, F_o, A_03);
		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_06, F_oo, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_03, auxiliaryMatrix);

		// находим матрицу А_4
		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_07, F_o, A_04);
		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_08, F_oo, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_04, auxiliaryMatrix);

		// находим матрицу А_5
		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_09, F_o, A_05);
		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_10, F_oo, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, A_05, auxiliaryMatrix);

		// находим матрицу В
		MultiplicationTransposedMatrix(N_squared, N_squared, N_squared, F_o, F_o, B);
		MultiplicationTransposedMatrix(N_squared, N + 1, N_squared, F_oo, F_oo, auxiliaryMatrix);
		AdditionOfSquareMatrices(N_squared, B, auxiliaryMatrix);

		//добавляем alpha к диагонали
		for (size_t l = 0; l < N_squared; ++l)
		{
			A_00[l][l] += alpha;
			B[l][l] += alpha;
		}
		timeFinish = clock();
		d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
		cout << "Calculation time left side " << d << endl;
		timeStart = clock();
		//////////////////////////////////////////////////////////////////////////
		//строим правую часть СЛАУ основаного метода Ньютона
		//////////////////////////////////////////////////////////////////////////
		//
		// находим значения основного оператора F
		//
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				F_part_01[ii] = u_1[i][j];
				F_part_03[ii] = u_2[i][j];
				F_part_05[ii] = u_3[i][j];
				F_part_07[ii] = u_4[i][j];
				F_part_09[ii] = u_5[i][j];
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t q = 0; q < N; ++q)
					{
						if ((i != p) || (q != j))
						{
							F_part_01[ii] += a[i][j][p][q] * (xi[p][q] * u_1[p][q] - xi[i][j] * u_1[i][j]);
							F_part_03[ii] += a[i][j][p][q] * (xi[p][q] * u_2[p][q] - xi[i][j] * u_2[i][j]);
							F_part_05[ii] += a[i][j][p][q] * (xi[p][q] * u_3[p][q] - xi[i][j] * u_3[i][j]);
							F_part_07[ii] += a[i][j][p][q] * (xi[p][q] * u_4[p][q] - xi[i][j] * u_4[i][j]);
							F_part_09[ii] += a[i][j][p][q] * (xi[p][q] * u_5[p][q] - xi[i][j] * u_5[i][j]);

						}
					}
				}
				F_part_01[ii] += b[i][j] * xi[i][j] * u_1[i][j];
				F_part_01[ii] -= f_Source_01[i][j];
				F_part_03[ii] += b[i][j] * xi[i][j] * u_2[i][j];
				F_part_03[ii] -= f_Source_03[i][j];
				F_part_05[ii] += b[i][j] * xi[i][j] * u_3[i][j];
				F_part_05[ii] -= f_Source_05[i][j];
				F_part_07[ii] += b[i][j] * xi[i][j] * u_4[i][j];
				F_part_07[ii] -= f_Source_07[i][j];
				F_part_09[ii] += b[i][j] * xi[i][j] * u_5[i][j];
				F_part_09[ii] -= f_Source_09[i][j];
			}
		}
		for (size_t j = 0; j <= N; ++j)
		{
			F_part_02[j] = overline_u_1[j] - f_Source_02[j];
			F_part_04[j] = overline_u_2[j] - f_Source_04[j];
			F_part_06[j] = overline_u_3[j] - f_Source_06[j];
			F_part_08[j] = overline_u_4[j] - f_Source_08[j];
			F_part_10[j] = overline_u_5[j] - f_Source_10[j];
			for (size_t p = 0; p<N; p++)
			{
				for (size_t q = 0; q < N; ++q)
				{
					F_part_02[j] += overline_a[j][p][q] * xi[p][q] * u_1[p][q];
					F_part_04[j] += overline_a[j][p][q] * xi[p][q] * u_2[p][q];
					F_part_06[j] += overline_a[j][p][q] * xi[p][q] * u_3[p][q];
					F_part_08[j] += overline_a[j][p][q] * xi[p][q] * u_4[p][q];
					F_part_10[j] += overline_a[j][p][q] * xi[p][q] * u_5[p][q];
				}
			}
		}
		//
		// находим F^{\prime}(x)*x и вычитаем F(x) 
		//
		// перенумерация xi
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				size_t l = i * (N + 1) + j;
				numbered_xi[l] = xi[i][j];
			}
		}
		//
		//первый источник 	
		//
		// перенумерация u^{1}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				size_t l = i * (N + 1) + j;
				numbered_u_1[l] = u_1[i][j];
			}
		}
		MultiplicationMatrixVector(N_squared, N_squared, F_01, numbered_xi, supportingVector_square);
		for (size_t l = 0; l < N_squared; ++l)
		{
			F_part_01[l] = supportingVector_square[l] - F_part_01[l];
		}
		MultiplicationMatrixVector(N_squared, N_squared, F_o, numbered_u_1, supportingVector_square);
		for (size_t l = 0; l < N_squared; ++l)
		{
			F_part_01[l] += supportingVector_square[l];
		}
		MultiplicationMatrixVector(N + 1, N_squared, F_02, numbered_xi, supportingVector);
		for (size_t l = 0; l <= N; ++l)
		{
			F_part_02[l] = supportingVector_square[l] - F_part_02[l];
		}
		MultiplicationMatrixVector(N + 1, N_squared, F_oo, numbered_u_1, supportingVector);
		for (size_t l = 0; l <= N; ++l)
		{
			F_part_02[l] = F_part_02[l] + supportingVector_square[l];
		}
		//
		//второй источник
		//
		// перенумерация u^{2}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				size_t l = i * (N + 1) + j;
				numbered_u_2[l] = u_2[i][j];
			}
		}
		MultiplicationMatrixVector(N_squared, N_squared, F_03, numbered_xi, supportingVector_square);
		for (size_t l = 0; l < N_squared; ++l)
		{
			F_part_03[l] = supportingVector_square[l] - F_part_03[l];
		}
		MultiplicationMatrixVector(N_squared, N_squared, F_o, numbered_u_2, supportingVector_square);
		for (size_t l = 0; l < N_squared; ++l)
		{
			F_part_03[l] += supportingVector_square[l];
		}
		MultiplicationMatrixVector(N + 1, N_squared, F_04, numbered_xi, supportingVector);
		for (size_t l = 0; l <= N; ++l)
		{
			F_part_04[l] = supportingVector_square[l] - F_part_04[l];
		}
		MultiplicationMatrixVector(N + 1, N_squared, F_oo, numbered_u_2, supportingVector);
		for (size_t l = 0; l <= N; ++l)
		{
			F_part_04[l] = F_part_04[l] + supportingVector_square[l];
		}
		//
		//третий источник
		//
		// перенумерация u^{3}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				size_t l = i * (N + 1) + j;
				numbered_u_3[l] = u_3[i][j];
			}
		}
		MultiplicationMatrixVector(N_squared, N_squared, F_05, numbered_xi, supportingVector_square);
		for (size_t l = 0; l < N_squared; ++l)
		{
			F_part_05[l] = supportingVector_square[l] - F_part_05[l];
		}
		MultiplicationMatrixVector(N_squared, N_squared, F_o, numbered_u_3, supportingVector_square);
		for (size_t l = 0; l <N_squared; ++l)
		{
			F_part_05[l] += supportingVector_square[l];
		}
		MultiplicationMatrixVector(N + 1, N_squared, F_06, numbered_xi, supportingVector);
		for (size_t l = 0; l <= N; ++l)
		{
			F_part_06[l] = supportingVector_square[l] - F_part_06[l];
		}
		MultiplicationMatrixVector(N + 1, N_squared, F_oo, numbered_u_3, supportingVector);
		for (size_t l = 0; l <= N; ++l)
		{
			F_part_06[l] = F_part_06[l] + supportingVector_square[l];
		}
		//
		// четвертый источник
		//
		// перенумерация u^{4}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				size_t l = i * (N + 1) + j;
				numbered_u_4[l] = u_4[i][j];
			}
		}
		MultiplicationMatrixVector(N_squared, N_squared, F_07, numbered_xi, supportingVector_square);
		for (size_t l = 0; l < N_squared; ++l)
		{
			F_part_07[l] = supportingVector_square[l] - F_part_07[l];
		}
		MultiplicationMatrixVector(N_squared, N_squared, F_o, numbered_u_4, supportingVector_square);
		for (size_t l = 0; l < N_squared; ++l)
		{
			F_part_07[l] += supportingVector_square[l];
		}
		MultiplicationMatrixVector(N + 1, N_squared, F_08, numbered_xi, supportingVector);
		for (size_t l = 0; l <= N; ++l)
		{
			F_part_08[l] = supportingVector_square[l] - F_part_08[l];
		}
		MultiplicationMatrixVector(N + 1, N_squared, F_oo, numbered_u_4, supportingVector);
		for (size_t l = 0; l <= N; ++l)
		{
			F_part_08[l] = F_part_08[l] + supportingVector_square[l];
		}
		//
		// пятый источник
		//
		// перенумерация u^{5}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				size_t l = i * (N + 1) + j;
				numbered_u_5[l] = u_5[i][j];
			}
		}
		MultiplicationMatrixVector(N_squared, N_squared, F_09, numbered_xi, supportingVector_square);
		for (size_t l = 0; l < N_squared; ++l)
		{
			F_part_09[l] = supportingVector_square[l] - F_part_09[l];
		}
		MultiplicationMatrixVector(N_squared, N_squared, F_o, numbered_u_5, supportingVector_square);
		for (size_t l = 0; l < N_squared; ++l)
		{
			F_part_09[l] += supportingVector_square[l];
		}
		MultiplicationMatrixVector(N + 1, N_squared, F_10, numbered_xi, supportingVector);
		for (size_t l = 0; l <= N; ++l)
		{
			F_part_10[l] = supportingVector_square[l] - F_part_10[l];
		}
		MultiplicationMatrixVector(N + 1, N_squared, F_oo, numbered_u_5, supportingVector);
		for (size_t l = 0; l <= N; ++l)
		{
			F_part_10[l] = F_part_10[l] + supportingVector_square[l];
		}

		// находим окончательно правую часть b0, b1,...
		// b0
		MultiplicationTransposedMatrixVector(N_squared, N_squared, F_01, F_part_01, b0);
		MultiplicationTransposedMatrixVector(N_squared, N + 1, F_02, F_part_02, supportingVector_square);
		AdditionOfSquareVector(N_squared, b0, supportingVector_square);

		MultiplicationTransposedMatrixVector(N_squared, N_squared, F_03, F_part_03, supportingVector_square);
		AdditionOfSquareVector(N_squared, b0, supportingVector_square);
		MultiplicationTransposedMatrixVector(N_squared, N + 1, F_04, F_part_04, supportingVector_square);
		AdditionOfSquareVector(N_squared, b0, supportingVector_square);

		MultiplicationTransposedMatrixVector(N_squared, N_squared, F_05, F_part_05, supportingVector_square);
		AdditionOfSquareVector(N_squared, b0, supportingVector_square);
		MultiplicationTransposedMatrixVector(N_squared, N + 1, F_06, F_part_06, supportingVector_square);
		AdditionOfSquareVector(N_squared, b0, supportingVector_square);

		MultiplicationTransposedMatrixVector(N_squared, N_squared, F_07, F_part_07, supportingVector_square);
		AdditionOfSquareVector(N_squared, b0, supportingVector_square);
		MultiplicationTransposedMatrixVector(N_squared, N + 1, F_08, F_part_08, supportingVector_square);
		AdditionOfSquareVector(N_squared, b0, supportingVector_square);

		MultiplicationTransposedMatrixVector(N_squared, N_squared, F_09, F_part_09, supportingVector_square);
		AdditionOfSquareVector(N_squared, b0, supportingVector_square);
		MultiplicationTransposedMatrixVector(N_squared, N + 1, F_10, F_part_10, supportingVector_square);
		AdditionOfSquareVector(N_squared, b0, supportingVector_square);

		//b1
		MultiplicationTransposedMatrixVector(N_squared, N_squared, F_o, F_part_01, b1);
		MultiplicationTransposedMatrixVector(N_squared, N + 1, F_oo, F_part_02, supportingVector_square);
		AdditionOfSquareVector(N_squared, b1, supportingVector_square);

		//b2
		MultiplicationTransposedMatrixVector(N_squared, N_squared, F_o, F_part_03, b2);
		MultiplicationTransposedMatrixVector(N_squared, N + 1, F_oo, F_part_04, supportingVector_square);
		AdditionOfSquareVector(N_squared, b2, supportingVector_square);

		//b3
		MultiplicationTransposedMatrixVector(N_squared, N_squared, F_o, F_part_05, b3);
		MultiplicationTransposedMatrixVector(N_squared, N + 1, F_oo, F_part_06, supportingVector_square);
		AdditionOfSquareVector(N_squared, b3, supportingVector_square);

		//b4
		MultiplicationTransposedMatrixVector(N_squared, N_squared, F_o, F_part_07, b4);
		MultiplicationTransposedMatrixVector(N_squared, N + 1, F_oo, F_part_08, supportingVector_square);
		AdditionOfSquareVector(N_squared, b4, supportingVector_square);

		//b5
		MultiplicationTransposedMatrixVector(N_squared, N_squared, F_o, F_part_09, b5);
		MultiplicationTransposedMatrixVector(N_squared, N + 1, F_oo, F_part_10, supportingVector_square);
		AdditionOfSquareVector(N_squared, b5, supportingVector_square);
		timeFinish = clock();
		d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
		cout << "Calculation time right side " << d << endl;
		timeStart = clock();

		// составляем основное уравнение для метода Ньютона
		// Нахождим xi
		// находим обратную матрицу для B
		ComputeInverseMatrixGaussJordan(N_squared, B, inverseMatrixB);
		//для левой части уравнения с xi все складываем в A_00

		MultiplicationMatrixBlock(N_squared, N_squared, N_squared, A_01, inverseMatrixB, auxiliaryMatrix);
		MultiplicationMatrixTransposed(N_squared, N_squared, N_squared, auxiliaryMatrix, A_01, secondAuxiliaryMatrix);
		SubtractionOfSquareMatrices(N_squared, A_00, secondAuxiliaryMatrix);

		MultiplicationMatrixBlock(N_squared, N_squared, N_squared, A_02, inverseMatrixB, auxiliaryMatrix);
		MultiplicationMatrixTransposed(N_squared, N_squared, N_squared, auxiliaryMatrix, A_02, secondAuxiliaryMatrix);
		SubtractionOfSquareMatrices(N_squared, A_00, secondAuxiliaryMatrix);

		MultiplicationMatrixBlock(N_squared, N_squared, N_squared, A_03, inverseMatrixB, auxiliaryMatrix);
		MultiplicationMatrixTransposed(N_squared, N_squared, N_squared, auxiliaryMatrix, A_03, secondAuxiliaryMatrix);
		SubtractionOfSquareMatrices(N_squared, A_00, secondAuxiliaryMatrix);

		MultiplicationMatrixBlock(N_squared, N_squared, N_squared, A_04, inverseMatrixB, auxiliaryMatrix);
		MultiplicationMatrixTransposed(N_squared, N_squared, N_squared, auxiliaryMatrix, A_04, secondAuxiliaryMatrix);
		SubtractionOfSquareMatrices(N_squared, A_00, secondAuxiliaryMatrix);

		MultiplicationMatrixBlock(N_squared, N_squared, N_squared, A_05, inverseMatrixB, auxiliaryMatrix);
		MultiplicationMatrixTransposed(N_squared, N_squared, N_squared, auxiliaryMatrix, A_05, secondAuxiliaryMatrix);
		SubtractionOfSquareMatrices(N_squared, A_00, secondAuxiliaryMatrix);

		//для правой части уравнения с xi все складываем в b0 
		MultiplicationMatrixVector(N_squared, N_squared, inverseMatrixB, b1, supportingVector_square);
		MultiplicationMatrixVector(N_squared, N_squared, A_01, supportingVector_square, secondSupportingVector_square);
		SubtractionOfSquareVector(N_squared, b0, secondSupportingVector_square);

		MultiplicationMatrixVector(N_squared, N_squared, inverseMatrixB, b2, supportingVector_square);
		MultiplicationMatrixVector(N_squared, N_squared, A_02, supportingVector_square, secondSupportingVector_square);
		SubtractionOfSquareVector(N_squared, b0, secondSupportingVector_square);

		MultiplicationMatrixVector(N_squared, N_squared, inverseMatrixB, b3, supportingVector_square);
		MultiplicationMatrixVector(N_squared, N_squared, A_03, supportingVector_square, secondSupportingVector_square);
		SubtractionOfSquareVector(N_squared, b0, secondSupportingVector_square);

		MultiplicationMatrixVector(N_squared, N_squared, inverseMatrixB, b4, supportingVector_square);
		MultiplicationMatrixVector(N_squared, N_squared, A_04, supportingVector_square, secondSupportingVector_square);
		SubtractionOfSquareVector(N_squared, b0, secondSupportingVector_square);

		MultiplicationMatrixVector(N_squared, N_squared, inverseMatrixB, b5, supportingVector_square);
		MultiplicationMatrixVector(N_squared, N_squared, A_05, supportingVector_square, secondSupportingVector_square);
		SubtractionOfSquareVector(N_squared, b0, secondSupportingVector_square);

		//  находим xi
		SolveSlauGaussa(A_00, N_squared, b0, numbered_xi);

		//находим u^{1}
		MultiplicationTransposedMatrixVector(N_squared, N_squared, A_01, numbered_xi, supportingVector_square);
		SubtractionOfSquareVector(N_squared, b1, supportingVector_square);
		MultiplicationMatrixVector(N_squared, N_squared, inverseMatrixB, b1, numbered_u_1);

		//находим u^{2}
		MultiplicationTransposedMatrixVector(N_squared, N_squared, A_02, numbered_xi, supportingVector_square);
		SubtractionOfSquareVector(N_squared, b2, supportingVector_square);
		MultiplicationMatrixVector(N_squared, N_squared, inverseMatrixB, b2, numbered_u_2);

		//находим u^{3}
		MultiplicationTransposedMatrixVector(N_squared, N_squared, A_03, numbered_xi, supportingVector_square);
		SubtractionOfSquareVector(N_squared, b3, supportingVector_square);
		MultiplicationMatrixVector(N_squared, N_squared, inverseMatrixB, b3, numbered_u_3);

		//находим u^{4}
		MultiplicationTransposedMatrixVector(N_squared, N_squared, A_04, numbered_xi, supportingVector_square);
		SubtractionOfSquareVector(N_squared, b4, supportingVector_square);
		MultiplicationMatrixVector(N_squared, N_squared, inverseMatrixB, b4, numbered_u_4);

		//находим u^{5}
		MultiplicationTransposedMatrixVector(N_squared, N_squared, A_05, numbered_xi, supportingVector_square);
		SubtractionOfSquareVector(N_squared, b5, supportingVector_square);
		MultiplicationMatrixVector(N_squared, N_squared, inverseMatrixB, b5, numbered_u_5);

		timeFinish = clock();
		d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
		cout << "Calculation time solutions " << d << endl;

		// изменяем alpha для следующей итерации
		alpha = alpha * multiplier;
		// Обратная перенумерация u^{1}, u^{2}, u^{3}, u^{4}, u^{5}, xi
		for (size_t l = 0; l < N_squared; ++l)
		{
			size_t coordinate_x = l / (N + 1);
			size_t coordinate_y = l % (N + 1);
			u_1[coordinate_x][coordinate_y] = numbered_u_1[l];
			u_2[coordinate_x][coordinate_y] = numbered_u_2[l];
			u_3[coordinate_x][coordinate_y] = numbered_u_3[l];
			u_4[coordinate_x][coordinate_y] = numbered_u_4[l];
			u_5[coordinate_x][coordinate_y] = numbered_u_5[l];
			xi[coordinate_x][coordinate_y] = numbered_xi[l];
		}
		// проекция xi>=0
		for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
		{
			for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
			{
				if (real(xi[i][j]) <= 0)
				{
					xi[i][j] = 0.0;
				}
			}
		}
		if (iteration == 0)
		{
			ofstream f_xi("approximate_xi_1.txt");
			for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
			{
				for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
				{
					f_xi << fixed << setprecision(6) << real(xi[i][j]) << " ";
				}
			}
			f_xi.close();
		}
		if (iteration == 1)
		{
			ofstream f_xi("approximate_xi_2.txt");
			for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
			{
				for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
				{
					f_xi << fixed << setprecision(6) << real(xi[i][j]) << " ";
				}
			}
			f_xi.close();
		}
		if (iteration == 2)
		{
			ofstream f_xi("approximate_xi_3.txt");
			for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
			{
				for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
				{
					f_xi << fixed << setprecision(6) << real(xi[i][j]) << " ";
				}
			}
			f_xi.close();
		}
		if (iteration == 3)
		{
			ofstream f_xi("approximate_xi_4.txt");
			for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
			{
				for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
				{
					f_xi << fixed << setprecision(6) << real(xi[i][j]) << " ";
				}
			}
			f_xi.close();
		}
		if (iteration == 4)
		{
			ofstream f_xi("approximate_xi_5.txt");
			for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
			{
				for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
				{
					f_xi << fixed << setprecision(6) << real(xi[i][j]) << " ";
				}
			}
			f_xi.close();
		}
	} // конец итераций

	//печатаем приближенное решение в файл
	ofstream f_xi("approximate_xi.txt");
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			f_xi << fixed << setprecision(6) << real(xi[i][j]) << " ";
		}
	}
	f_xi.close();



	//
	// освобождение памяти	
	//
	// освобождаем основные массивы
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p <= N; ++p)
			{
				delete[] a[i][j][p];
			}
			delete[] a[i][j];
		}
		delete[] a[i];
	}
	delete[] a;

	for (size_t i = 0; i <= N; ++i)
	{
		delete[] b[i];
	}
	delete[] b;

	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			delete[] overline_a[i][j];
		}
		delete[] overline_a[i];
	}
	delete[] overline_a;

	// очистка памяти для левой части
	for (size_t i = 0; i <= N; ++i)
	{
		delete[] f_Source_01[i];
		delete[] f_Source_03[i];
		delete[] f_Source_05[i];
		delete[] f_Source_07[i];
		delete[] f_Source_09[i];
	}
	delete[] f_Source_01;
	delete[] f_Source_03;
	delete[] f_Source_05;
	delete[] f_Source_07;
	delete[] f_Source_09;

	delete[] f_Source_02;
	delete[] f_Source_04;
	delete[] f_Source_06;
	delete[] f_Source_08;
	delete[] f_Source_10;

	// очистка памяти поля в приемниках
	delete[] overline_u_1;
	delete[] overline_u_2;
	delete[] overline_u_3;
	delete[] overline_u_4;
	delete[] overline_u_5;

	// очистка памяти от производных
	for (size_t i = 0; i<N_squared; ++i)
	{
		delete[] F_01[i];
		delete[] F_03[i];
		delete[] F_05[i];
		delete[] F_07[i];
		delete[] F_09[i];
		delete[] F_o[i];
	}
	delete[] F_01;
	delete[] F_03;
	delete[] F_05;
	delete[] F_07;
	delete[] F_09;
	delete[] F_o;

	for (size_t i = 0; i <= N; ++i)
	{
		delete[] F_02[i];
		delete[] F_04[i];
		delete[] F_06[i];
		delete[] F_08[i];
		delete[] F_10[i];
		delete[] F_oo[i];
	}
	delete[] F_02;
	delete[] F_04;
	delete[] F_06;
	delete[] F_08;
	delete[] F_10;
	delete[] F_oo;

	for (size_t i = 0; i < N_squared; ++i)
	{
		delete[] A_00[i];
		delete[] A_01[i];
		delete[] A_02[i];
		delete[] A_03[i];
		delete[] A_04[i];
		delete[] A_05[i];
		delete[] B[i];
		delete[] auxiliaryMatrix[i];
		delete[] secondAuxiliaryMatrix[i];
		delete[] inverseMatrixB[i];
	}
	delete[] A_00;
	delete[] A_01;
	delete[] A_02;
	delete[] A_03;
	delete[] A_04;
	delete[] A_05;
	delete[] B;
	delete[] auxiliaryMatrix;
	delete[] secondAuxiliaryMatrix;
	delete[] inverseMatrixB;

	delete[] F_part_01;
	delete[] F_part_02;
	delete[] F_part_03;
	delete[] F_part_04;
	delete[] F_part_05;
	delete[] F_part_06;
	delete[] F_part_07;
	delete[] F_part_08;
	delete[] F_part_09;
	delete[] F_part_10;

	delete[] b0;
	delete[] b1;
	delete[] b2;
	delete[] b3;
	delete[] b4;
	delete[] b5;

	for (size_t i = 0; i <= N; ++i)
	{
		delete[] u_1[i];
		delete[] u_2[i];
		delete[] u_3[i];
		delete[] u_4[i];
		delete[] u_5[i];
	}
	delete[] u_1;
	delete[] u_2;
	delete[] u_3;
	delete[] u_4;
	delete[] u_5;

	for (size_t i = 0; i <= N; ++i)
	{
		delete[] xi[i];
	}
	delete[] xi;

	delete[] numbered_u_1;
	delete[] numbered_u_2;
	delete[] numbered_u_3;
	delete[] numbered_u_4;
	delete[] numbered_u_5;
	delete[] numbered_xi;
	delete[] supportingVector_square;
	delete[] secondSupportingVector_square;
	delete[] supportingVector;
}
