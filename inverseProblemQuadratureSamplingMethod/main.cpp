#include "stdafx.h"
#include "../directProblemQuadratureSamplingMethod/basicFunctions.h"
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
	vector<vector<complex<float>>> f_Source_01(N + 1, vector<complex<float>>(N + 1, complex<float>()));
	vector<complex<float>> f_Source_02(N + 1, complex<float>());
	vector<vector<complex<float>>> f_Source_03(N + 1, vector<complex<float>>(N + 1, complex<float>()));
	vector<complex<float>> f_Source_04(N + 1, complex<float>());
	vector<vector<complex<float>>> f_Source_05(N + 1, vector<complex<float>>(N + 1, complex<float>()));
	vector<complex<float>> f_Source_06(N + 1, complex<float>());
	vector<vector<complex<float>>> f_Source_07(N + 1, vector<complex<float>>(N + 1, complex<float>()));
	vector<complex<float>> f_Source_08(N + 1, complex<float>());
	vector<vector<complex<float>>> f_Source_09(N + 1, vector<complex<float>>(N + 1, complex<float>()));
	vector<complex<float>> f_Source_10(N + 1, complex<float>());

	// загружаем значения источника
	// первый источник
	{
		ifstream fileSource01("Source_01.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource01 >> f_Source_01[i][j];
			}
		}
		fileSource01.close();
		ifstream fileSource02("Source_02.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			fileSource02 >> f_Source_02[i];
		}
		fileSource02.close();
	}
	
	// второй источник
	{
		ifstream fileSource03("Source_03.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource03 >> f_Source_03[i][j];
			}
		}
		fileSource03.close();
		ifstream fileSource04("Source_04.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			fileSource04 >> f_Source_04[i];
		}
		fileSource04.close();
	}

	// третий источник
	{
		ifstream fileSource05("Source_05.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource05 >> f_Source_05[i][j];
			}
		}
		fileSource05.close();
		ifstream fileSource06("Source_06.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			fileSource06 >> f_Source_06[i];
		}
		fileSource06.close();
	}

	// четвертый источник
	{
		ifstream fileSource07("Source_07.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource07 >> f_Source_07[i][j];
			}
		}
		fileSource07.close();
		ifstream fileSource08("Source_08.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			fileSource08 >> f_Source_08[i];
		}
		fileSource08.close();
	}

	// пятый источник
	{
		ifstream fileSource09("Source_09.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource09 >> f_Source_09[i][j];
			}
		}
		fileSource09.close();
		ifstream fileSource10("Source_10.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			fileSource10 >> f_Source_10[i];
		}
		fileSource10.close();
	}

	//печатаем время работы
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Download Time arrays sources " << d << endl;
	timeStart = clock();

	// Выделяем память для поля в приемниках
	vector<complex<float>> overline_u_1(N + 1, complex<float>());
	vector<complex<float>> overline_u_2(N + 1, complex<float>());
	vector<complex<float>> overline_u_3(N + 1, complex<float>());
	vector<complex<float>> overline_u_4(N + 1, complex<float>());
	vector<complex<float>> overline_u_5(N + 1, complex<float>());

	// Загрузка акустического поля в приёмнике
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
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Loading time of the acoustic field in the receivers " << d << endl;
	timeStart = clock();

	// вспомогательные переменные
	const size_t N_squared = (N + 1) * (N + 1);

	//выделение памяти под массивы производных  F_1, F_2, ...
	vector<vector<complex<float>>> F_01(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_03(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_05(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_07(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_09(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_o(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_02(N + 1, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_04(N + 1, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_06(N + 1, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_08(N + 1, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_10(N + 1, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_oo(N + 1, vector<complex<float>>(N_squared, complex<float>()));

	//выделение памяти под массивы A и B
	vector<vector<complex<float>>> A_00(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> A_01(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> A_02(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> A_03(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> A_04(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> A_05(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> B(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> inverseMatrixB(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> auxiliaryMatrix(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> secondAuxiliaryMatrix(N_squared, vector<complex<float>>(N_squared, complex<float>()));

	// память для хранения значений основного оператора
	vector<complex<float>> F_part_01(N_squared, complex<float>());
	vector<complex<float>> F_part_02(N + 1, complex<float>());
	vector<complex<float>> F_part_03(N_squared, complex<float>());
	vector<complex<float>> F_part_04(N + 1, complex<float>());
	vector<complex<float>> F_part_05(N_squared, complex<float>());
	vector<complex<float>> F_part_06(N + 1, complex<float>());
	vector<complex<float>> F_part_07(N_squared, complex<float>());
	vector<complex<float>> F_part_08(N + 1, complex<float>());
	vector<complex<float>> F_part_09(N_squared, complex<float>());
	vector<complex<float>> F_part_10(N + 1, complex<float>());

	// память для b_0, b_1,...
	vector<complex<float>> b0(N_squared, complex<float>());
	vector<complex<float>> b1(N_squared, complex<float>());
	vector<complex<float>> b2(N_squared, complex<float>());
	vector<complex<float>> b3(N_squared, complex<float>());
	vector<complex<float>> b4(N_squared, complex<float>());
	vector<complex<float>> b5(N_squared, complex<float>());

	// память для u^(1), u^(2), u^(3)
	vector<vector<complex<float>>> u_1(N + 1, vector<complex<float>>(N + 1, complex<float>()));
	vector<vector<complex<float>>> u_2(N + 1, vector<complex<float>>(N + 1, complex<float>()));
	vector<vector<complex<float>>> u_3(N + 1, vector<complex<float>>(N + 1, complex<float>()));
	vector<vector<complex<float>>> u_4(N + 1, vector<complex<float>>(N + 1, complex<float>()));
	vector<vector<complex<float>>> u_5(N + 1, vector<complex<float>>(N + 1, complex<float>()));

	// память для xi
	vector<vector<complex<float>>> xi(N + 1, vector<complex<float>>(N + 1, complex<float>()));

	// память для перенумерованных переменных и вспомогательного вектора
	vector<complex<float>> numbered_u_1(N_squared, complex<float>());
	vector<complex<float>> numbered_u_2(N_squared, complex<float>());
	vector<complex<float>> numbered_u_3(N_squared, complex<float>());
	vector<complex<float>> numbered_u_4(N_squared, complex<float>());
	vector<complex<float>> numbered_u_5(N_squared, complex<float>());

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
	for(size_t i = 0; i <= N; ++i)
	{
		for(size_t j = 0; j <= N; ++j)
		{
			u_1[i][j] = f_Source_01[i][j];
			u_2[i][j] = f_Source_03[i][j];
			u_3[i][j] = f_Source_05[i][j];
			u_3[i][j] = f_Source_07[i][j];
			u_3[i][j] = f_Source_09[i][j];
		}
	}
	// начальное приближение xi
	for(size_t i = 0; i <= N; ++i)
	{
		for(size_t j = 0; j <= N; ++j)
		{
			xi[i][j] = 0.1f;
		}
	}

	size_t coordinate_x;
	size_t coordinate_y;
	size_t ii, jj;
	complex<float> sumOfTheCoefficients;

	// начало основных итераций
	for(size_t iteration = 0; iteration < numberOfIterations; ++iteration)
	{
		cout << endl;
		cout << "Iteration number " << (iteration + 1) << endl;
		cout << "alpha= " << alpha << endl;
		timeStart = clock();
		//////////////////////////////////////////////////////////////////////////
		//строим левую часть СЛАУ основаного метода Ньютона
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
				F_o[ii][ii] = (1.0f, 0.0f) + xi[i][j] * (b[i][j] - sumOfTheCoefficients);
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
					F_02[j][jj] = overline_a[j][p][q] * u_1[p][q];
					F_04[j][jj] = overline_a[j][p][q] * u_2[p][q];
					F_06[j][jj] = overline_a[j][p][q] * u_3[p][q];
					F_08[j][jj] = overline_a[j][p][q] * u_4[p][q];
					F_10[j][jj] = overline_a[j][p][q] * u_5[p][q];
					F_oo[j][jj] = overline_a[j][p][q] * xi[p][q];
				}
			}
		}

		// находим матрицы А и В
		// находим матрицы А_0
		MultTransposedMatrix(F_01, F_01, A_00);
		MultTransposedMatrix(F_02, F_02, auxiliaryMatrix);
		AddSquareMatrices(A_00, auxiliaryMatrix);

		MultTransposedMatrix(F_03, F_03, auxiliaryMatrix);
		AddSquareMatrices(A_00, auxiliaryMatrix);

		MultTransposedMatrix(F_04, F_04, auxiliaryMatrix);
		AddSquareMatrices(A_00, auxiliaryMatrix);

		MultTransposedMatrix(F_05, F_05, auxiliaryMatrix);
		AddSquareMatrices(A_00, auxiliaryMatrix);

		MultTransposedMatrix(F_06, F_06, auxiliaryMatrix);
		AddSquareMatrices(A_00, auxiliaryMatrix);

		MultTransposedMatrix(F_07, F_07, auxiliaryMatrix);
		AddSquareMatrices(A_00, auxiliaryMatrix);

		MultTransposedMatrix(F_08, F_08, auxiliaryMatrix);
		AddSquareMatrices(A_00, auxiliaryMatrix);

		MultTransposedMatrix(F_09, F_09, auxiliaryMatrix);
		AddSquareMatrices(A_00, auxiliaryMatrix);

		MultTransposedMatrix(F_10, F_10, auxiliaryMatrix);
		AddSquareMatrices(A_00, auxiliaryMatrix);

		// находим матрицу А_1
		MultTransposedMatrix(F_01, F_o, A_01);
		MultTransposedMatrix(F_02, F_oo, auxiliaryMatrix);
		AddSquareMatrices(A_01, auxiliaryMatrix);

		// находим матрицу А_2
		MultTransposedMatrix(F_03, F_o, A_02);
		MultTransposedMatrix(F_04, F_oo, auxiliaryMatrix);
		AddSquareMatrices(A_05, auxiliaryMatrix);

		// находим матрицу А_3
		MultTransposedMatrix(F_05, F_o, A_03);
		MultTransposedMatrix(F_06, F_oo, auxiliaryMatrix);
		AddSquareMatrices(A_03, auxiliaryMatrix);

		// находим матрицу А_4
		MultTransposedMatrix(F_07, F_o, A_04);
		MultTransposedMatrix(F_08, F_oo, auxiliaryMatrix);
		AddSquareMatrices(A_04, auxiliaryMatrix);

		// находим матрицу А_5
		MultTransposedMatrix(F_09, F_o, A_05);
		MultTransposedMatrix(F_10, F_oo, auxiliaryMatrix);
		AddSquareMatrices(A_05, auxiliaryMatrix);

		// находим матрицу В
		MultTransposedMatrix(F_o, F_o, B);
		MultTransposedMatrix(F_oo, F_oo, auxiliaryMatrix);
		AddSquareMatrices(B, auxiliaryMatrix);

		//добавляем alpha к диагонали
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			A_00[ii][ii] += alpha;
			B[ii][ii] += alpha;
		}
		timeFinish = clock();
		d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
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
			for (size_t p = 0; p < N; ++p)
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

		//первый источник 
		// перенумерация u^{1}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				numbered_u_1[ii] = u_1[i][j];
			}
		}
		MultMatrixVector(F_01, numbered_xi, supportingVector_square);
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			F_part_01[ii] = supportingVector_square[ii] - F_part_01[ii];
		}
		MultMatrixVector(F_o, numbered_u_1, supportingVector_square);
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			F_part_01[ii] += supportingVector_square[ii];
		}
		MultMatrixVector(F_02, numbered_xi, supportingVector);
		for (size_t ii = 0; ii <= N; ++ii)
		{
			F_part_02[ii] = supportingVector_square[ii] - F_part_02[ii];
		}
		MultMatrixVector(F_oo, numbered_u_1, supportingVector);
		for (size_t ii = 0; ii <= N; ++ii)
		{
			F_part_02[ii] += supportingVector_square[ii];
		}

		//второй источник
		// перенумерация u^{2}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				numbered_u_2[ii] = u_2[i][j];
			}
		}
		MultMatrixVector(F_03, numbered_xi, supportingVector_square);
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			F_part_03[ii] = supportingVector_square[ii] - F_part_03[ii];
		}
		MultMatrixVector(F_o, numbered_u_2, supportingVector_square);
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			F_part_03[ii] += supportingVector_square[ii];
		}
		MultMatrixVector(F_04, numbered_xi, supportingVector);
		for (size_t ii = 0; ii <= N; ++ii)
		{
			F_part_04[ii] = supportingVector_square[ii] - F_part_04[ii];
		}
		MultMatrixVector(F_oo, numbered_u_2, supportingVector);
		for (size_t ii = 0; ii <= N; ++ii)
		{
			F_part_04[ii] += supportingVector_square[ii];
		}

		//третий источник
		// перенумерация u^{3}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				numbered_u_3[ii] = u_3[i][j];
			}
		}
		MultMatrixVector(F_05, numbered_xi, supportingVector_square);
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			F_part_05[ii] = supportingVector_square[ii] - F_part_05[ii];
		}
		MultMatrixVector(F_o, numbered_u_3, supportingVector_square);
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			F_part_05[ii] += supportingVector_square[ii];
		}
		MultMatrixVector(F_06, numbered_xi, supportingVector);
		for (size_t ii = 0; ii <= N; ++ii)
		{
			F_part_06[ii] = supportingVector_square[ii] - F_part_06[ii];
		}
		MultMatrixVector(F_oo, numbered_u_3, supportingVector);
		for (size_t ii = 0; ii <= N; ++ii)
		{
			F_part_06[ii] += supportingVector_square[ii];
		}

		// четвертый источник
		// перенумерация u^{4}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				numbered_u_4[ii] = u_4[i][j];
			}
		}
		MultMatrixVector(F_07, numbered_xi, supportingVector_square);
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			F_part_07[ii] = supportingVector_square[ii] - F_part_07[ii];
		}
		MultMatrixVector(F_o, numbered_u_4, supportingVector_square);
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			F_part_07[ii] += supportingVector_square[ii];
		}
		MultMatrixVector(F_08, numbered_xi, supportingVector);
		for (size_t ii = 0; ii <= N; ++ii)
		{
			F_part_08[ii] = supportingVector_square[ii] - F_part_08[ii];
		}
		MultMatrixVector(F_oo, numbered_u_4, supportingVector);
		for (size_t ii = 0; ii <= N; ++ii)
		{
			F_part_08[ii] += supportingVector_square[ii];
		}

		// пятый источник
		// перенумерация u^{5}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				numbered_u_5[ii] = u_5[i][j];
			}
		}
		MultMatrixVector(F_09, numbered_xi, supportingVector_square);
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			F_part_09[ii] = supportingVector_square[ii] - F_part_09[ii];
		}
		MultMatrixVector(F_o, numbered_u_5, supportingVector_square);
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			F_part_09[ii] += supportingVector_square[ii];
		}
		MultMatrixVector(F_10, numbered_xi, supportingVector);
		for (size_t ii = 0; ii <= N; ++ii)
		{
			F_part_10[ii] = supportingVector_square[ii] - F_part_10[ii];
		}
		MultMatrixVector(F_oo, numbered_u_5, supportingVector);
		for (size_t ii = 0; ii <= N; ++ii)
		{
			F_part_10[ii] += supportingVector_square[ii];
		}

		// находим окончательно правую часть b0, b1,...
		// b0
		MultTransposedMatrixVector(F_01, F_part_01, b0);
		MultTransposedMatrixVector(F_02, F_part_02, supportingVector_square);
		AddVectors(b0, supportingVector_square);

		MultTransposedMatrixVector(F_03, F_part_03, supportingVector_square);
		AddVectors(b0, supportingVector_square);
		MultTransposedMatrixVector(F_04, F_part_04, supportingVector_square);
		AddVectors(b0, supportingVector_square);

		MultTransposedMatrixVector(F_05, F_part_05, supportingVector_square);
		AddVectors(b0, supportingVector_square);
		MultTransposedMatrixVector(F_06, F_part_06, supportingVector_square);
		AddVectors(b0, supportingVector_square);

		MultTransposedMatrixVector(F_07, F_part_07, supportingVector_square);
		AddVectors(b0, supportingVector_square);
		MultTransposedMatrixVector(F_08, F_part_08, supportingVector_square);
		AddVectors(b0, supportingVector_square);

		MultTransposedMatrixVector(F_09, F_part_09, supportingVector_square);
		AddVectors(b0, supportingVector_square);
		MultTransposedMatrixVector(F_10, F_part_10, supportingVector_square);
		AddVectors(b0, supportingVector_square);

		//b1
		MultTransposedMatrixVector(F_o, F_part_01, b1);
		MultTransposedMatrixVector(F_oo, F_part_02, supportingVector_square);
		AddVectors(b1, supportingVector_square);

		//b2
		MultTransposedMatrixVector(F_o, F_part_03, b2);
		MultTransposedMatrixVector(F_oo, F_part_04, supportingVector_square);
		AddVectors(b2, supportingVector_square);

		//b3
		MultTransposedMatrixVector(F_o, F_part_05, b3);
		MultTransposedMatrixVector(F_oo, F_part_06, supportingVector_square);
		AddVectors(b3, supportingVector_square);

		//b4
		MultTransposedMatrixVector(F_o, F_part_07, b4);
		MultTransposedMatrixVector(F_oo, F_part_08, supportingVector_square);
		AddVectors(b4, supportingVector_square);

		//b5
		MultTransposedMatrixVector(F_o, F_part_09, b5);
		MultTransposedMatrixVector(F_oo, F_part_10, supportingVector_square);
		AddVectors(b5, supportingVector_square);
		timeFinish = clock();
		d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
		cout << "Calculation time right side " << d << endl;
		timeStart = clock();

		// составляем основное уравнение для метода Ньютона

		// Находим xi
		// находим обратную матрицу для B
		InvertMatrix(B, inverseMatrixB);
		//для левой части уравнения с xi все складываем в A_00
		MultMatrix(A_01, inverseMatrixB, auxiliaryMatrix);
		MultMatrixTransposed(auxiliaryMatrix, A_01, secondAuxiliaryMatrix);
		SubSquareMatrices(A_00, secondAuxiliaryMatrix);

		MultMatrix(A_02, inverseMatrixB, auxiliaryMatrix);
		MultMatrixTransposed(auxiliaryMatrix, A_02, secondAuxiliaryMatrix);
		SubSquareMatrices(A_00, secondAuxiliaryMatrix);

		MultMatrix(A_03, inverseMatrixB, auxiliaryMatrix);
		MultMatrixTransposed(auxiliaryMatrix, A_03, secondAuxiliaryMatrix);
		SubSquareMatrices(A_00, secondAuxiliaryMatrix);

		MultMatrix(A_04, inverseMatrixB, auxiliaryMatrix);
		MultMatrixTransposed(auxiliaryMatrix, A_04, secondAuxiliaryMatrix);
		SubSquareMatrices(A_00, secondAuxiliaryMatrix);

		MultMatrix(A_05, inverseMatrixB, auxiliaryMatrix);
		MultMatrixTransposed(auxiliaryMatrix, A_05, secondAuxiliaryMatrix);
		SubSquareMatrices(A_00, secondAuxiliaryMatrix);
		//
		//для правой части уравнения с xi все складываем в b0 
		//
		MultMatrixVector(inverseMatrixB, b1, supportingVector_square);
		MultMatrixVector(A_01, supportingVector_square, secondSupportingVector_square);
		SubVectors(b0, secondSupportingVector_square);

		MultMatrixVector(inverseMatrixB, b2, supportingVector_square);
		MultMatrixVector(A_02, supportingVector_square, secondSupportingVector_square);
		SubVectors(b0, secondSupportingVector_square);

		MultMatrixVector(inverseMatrixB, b3, supportingVector_square);
		MultMatrixVector(A_03, supportingVector_square, secondSupportingVector_square);
		SubVectors(b0, secondSupportingVector_square);

		MultMatrixVector(inverseMatrixB, b4, supportingVector_square);
		MultMatrixVector(A_04, supportingVector_square, secondSupportingVector_square);
		SubVectors(b0, secondSupportingVector_square);

		MultMatrixVector(inverseMatrixB, b5, supportingVector_square);
		MultMatrixVector(A_05, supportingVector_square, secondSupportingVector_square);
		SubVectors(b0, secondSupportingVector_square);

		// находим xi
		SolveSlauGaussa(A_00, b0, numbered_xi);

		//находим u^{1}
		MultTransposedMatrixVector(A_01, numbered_xi, supportingVector_square);
		SubVectors(b1, supportingVector_square);
		MultMatrixVector(inverseMatrixB, b1, numbered_u_1);

		//находим u^{2}
		MultTransposedMatrixVector(A_02, numbered_xi, supportingVector_square);
		SubVectors(b2, supportingVector_square);
		MultMatrixVector(inverseMatrixB, b2, numbered_u_2);

		//находим u^{3}
		MultTransposedMatrixVector(A_03, numbered_xi, supportingVector_square);
		SubVectors(b3, supportingVector_square);
		MultMatrixVector(inverseMatrixB, b3, numbered_u_3);

		//находим u^{4}
		MultTransposedMatrixVector(A_04, numbered_xi, supportingVector_square);
		SubVectors(b4, supportingVector_square);
		MultMatrixVector(inverseMatrixB, b4, numbered_u_4);

		//находим u^{5}
		MultTransposedMatrixVector(A_05, numbered_xi, supportingVector_square);
		SubVectors(b5, supportingVector_square);
		MultMatrixVector(inverseMatrixB, b5, numbered_u_5);

		timeFinish = clock();
		d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
		cout << "Calculation time solutions " << d << endl;

		// изменяем alpha для следующей итерации
		alpha = alpha * q;
		// Обратная перенумерация u^{1}, u^{2}, u^{3}, u^{4}, u^{5}, xi
		for (size_t ii = 0; ii < N_squared; ++ii)
		{
			coordinate_x = ii / (N + 1);
			coordinate_y = ii % (N + 1);
			u_1[coordinate_x][coordinate_y] = numbered_u_1[ii];
			u_2[coordinate_x][coordinate_y] = numbered_u_2[ii];
			u_3[coordinate_x][coordinate_y] = numbered_u_3[ii];
			u_4[coordinate_x][coordinate_y] = numbered_u_4[ii];
			u_5[coordinate_x][coordinate_y] = numbered_u_5[ii];
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
		ofstream f_xi("approximate_xi_" + to_string(iteration + 1) + ".txt");
		f_xi << fixed << setprecision(6);
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				f_xi << real(xi[i][j]) << " ";
			}
		}
		f_xi.close();
	}

	d = (float)(timeFinish - timeBegin) / CLOCKS_PER_SEC;
	cout << "The total time of the program " << d << endl;

	return 0;
}
