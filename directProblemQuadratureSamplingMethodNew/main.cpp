﻿#include "stdafx.h"
#include "basicFunctions.h"
#include "Sources.h"
#include "taskData.h"
#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"

using namespace std;

int main()
{
	const size_t N = NUMBER_PARTITION_POINTS;
	const float h = (float)DOMAIN_IN_HOMOGENEITY / N;

	const Source source;

	// выделение памяти для акустического поля u
	vector<vector<complex<float>>> u(N + 1, vector<complex<float>>(N + 1, complex<float>()));

	// задание точного решения \xi
	vector<vector<float>> xi(N + 1, vector<float>(N + 1, 0.0f));
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			xi[i][j] = exp(-(i * h - 6.0f) * (i * h - 6.0f) - (j * h - 6.0f) * (j * h - 6.0f));
		}
	}
	//печатаем точное решение в файл
	ofstream f_xi("exact_xi.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_xi << fixed << setprecision(6) << xi[i][j] << " ";
		}
	}
	f_xi.close();

	// 
	// Начало вычислений основных матриц
	//
	// начало счета времени
	clock_t timeStart, timeFinish, timeBegin;
	timeBegin = clock();
	timeStart = clock();


	// выделение памяти под 4-х мерный "квадратный" комплексный массив
	vector<vector<vector<vector<complex<float>>>>> a(N + 1,
		vector<vector<vector<complex<float>>>>(N + 1, vector<vector<complex<float>>>(N + 1,
			vector<complex<float>>(N + 1, complex<float>()))));

	// выделение памяти под 3-х мерный "квадратный" комплексный массив
	vector<vector<vector<complex<float>>>> overline_a(N + 1, vector<vector<complex<float>>>(N + 1,
		vector<complex<float>>(N + 1, complex<float>())));

	// выделение памяти под 2-х мерный квадратный комплексный массив
	vector<vector<complex<float>>> b(N + 1, vector<complex<float>>(N + 1, complex<float>()));

	// счет индексов метода квадратур
	vector<float> index(N + 1);
	for (size_t i = 1; i < N; ++i)
	{
		if (i % 2 != 0)
		{
			index[i] = (float)4 / 3;
		}
		else
		{
			index[i] = (float)2 / 3;
		}
	}
	index[0] = (float)1 / 3;
	index[N] = (float)1 / 3;

	// нахождение массива a
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					if ((i != p) || (q != j))
					{
						a[i][j][p][q] = index[p] * index[q];
						a[i][j][p][q] = a[i][j][p][q] * G(i * h, j * h, p * h, q * h);
						a[i][j][p][q] = a[i][j][p][q] * h * h;
					}
				}
			}
		}
	}

	// нахождение массива overline_a
	for (size_t j = 0; j <= N; ++j)
	{
		for (size_t p = 0; p < N; ++p)
		{
			for (size_t q = 0; q < N; ++q)
			{
				overline_a[j][p][q] = index[p] * index[q];
				overline_a[j][p][q] = overline_a[j][p][q] * G(receiver, j * h, p * h, q * h);
				overline_a[j][p][q] = overline_a[j][p][q] * h * h;
			}
		}
	}

	// нахождение массива b
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					if (i != p)
					{
						b[i][j] += G(i * h, j * h, p * h, q * h);
					}
				}
				b[i][j] *= h * h;
			}
		}
	}

	//печатаем время работы
	timeFinish = clock();
	float d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Time calculation of basic matrices " << d << endl;
	timeStart = clock();

	// записываем массивы в файлы
	ofstream f_a("matrix_a.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					f_a << fixed << setprecision(6) << a[i][j][p][q] << " ";
				}
			}
		}
	}
	f_a.close();

	ofstream f_overline_a("matrix_overline_a.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		for (size_t p = 0; p < N; ++p)
		{
			for (size_t q = 0; q < N; ++q)
			{
				f_overline_a << fixed << setprecision(6) << overline_a[j][p][q] << " ";
			}
		}
	}
	f_overline_a.close();

	ofstream f_b("matrix_b.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_b << fixed << setprecision(6) << b[i][j] << " ";
		}
	}
	f_b.close();

	//печатаем время работы
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Download time major arrays " << d << endl;
	timeStart = clock();

	// счет функции источника в R и X
	ofstream fileSource("Source.txt");
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource << fixed << setprecision(6) << source.Function(source.node[count], i * h, j * h) << " ";
			}
		}

		for (size_t j = 0; j <= N; ++j)
		{
			fileSource << fixed << setprecision(6) << source.Function(source.node[count], receiver, j * h) << " ";
		}
	}
	fileSource.close();

	//печатаем время работы
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "The computation time of the source function " << d << endl;
	timeStart = clock();

	// для нахождения u^(1) составляем СЛАУ основная матрица * u^(1) = правой части
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEequation[ii]

	const size_t N_squared = (N + 1) * (N + 1);
	vector<complex<float>> rightPartEquation(N_squared, (0.0f, 0.0f));
	vector<complex<float>> numbered_u(N_squared);
	vector<vector<complex<float>>> substantiveMatrix(N_squared, vector<complex<float>>(N_squared, (0.0f, 0.0f)));
	vector<complex<float>> overline_u(N + 1, (0.0f, 0.0f));

	//счет основной матрицы
	size_t ii, jj;
	complex<float> sumOfTheCoefficients;
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
						substantiveMatrix[ii][jj] += a[i][j][p][q] * xi[i][j];
						sumOfTheCoefficients += a[i][j][p][q];
					}
				}
			}
			substantiveMatrix[ii][ii] += (1.0f, 0.0f);
			substantiveMatrix[ii][ii] -= sumOfTheCoefficients * xi[i][j];
			substantiveMatrix[ii][ii] += b[i][j] * xi[i][j];
		}
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "The computation time of the matrix inside the squared " << d << endl;
	timeStart = clock();

	// для каждого источника
	ofstream file_overline_u("matrix_overline_u.txt");
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		// нахождение правой части
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				rightPartEquation[ii] = source.Function(source.node[count], i * h, j * h);
			}
		}
		// нахождение u^{(count)}
		SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);

		// Обратная перенумерация
		size_t coordinate_x;
		size_t coordinate_y;
		for (size_t i = 0; i < N_squared; ++i)
		{
			coordinate_x = i / (N + 1);
			coordinate_y = i % (N + 1);
			u[coordinate_x][coordinate_y] = numbered_u[i];
		}
		timeFinish = clock();
		d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
		cout << "Finding the acoustic pressure in R for " << count + 1 << " source " << d << endl;
		timeStart = clock();

		// находим overline_u_0
		for (size_t i = 0; i <= N; ++i)
		{
			overline_u[i] = source.Function(source.node[count], receiver, i * h);			
		}
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					overline_u[j] -= overline_a[j][p][q] * xi[p][q] * u[p][q];
				}
			}
		}
		// печать overline_u_0^(count) в файл в одну строчку
		for (size_t j = 0; j <= N; ++j)
		{
			file_overline_u << fixed << setprecision(6) << overline_u[j] << " ";
		}
		timeFinish = clock();
		d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
		cout << "Finding the acoustic pressure in X for " << count + 1 << " source " << d << endl;
		timeStart = clock();
	}
	file_overline_u.close();

	timeFinish = clock();
	d = (float)(timeFinish - timeBegin) / CLOCKS_PER_SEC;
	cout << "The total time of the program " << d << endl;

	return 0;
}
