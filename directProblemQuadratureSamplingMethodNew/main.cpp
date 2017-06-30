#include "stdafx.h"
#include "basicFunctions.h"
#include "Sources.h"
#include "taskData.h"
#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"

using namespace std;

int main()
{
	const size_t N = NUMBER_PARTITION_POINT;
	const double h = (double)DOMAIN_IN_HOMOGENEITY / N;

	const Source source;

	// выделение памяти для акустического поля u
	vector<vector<complex<double>>> u(N + 1, vector<complex<double>>(N + 1, complex<double>()));

	// задание точного решения \xi
	vector<vector<double>> xi(N + 1, vector<double>(N + 1, 0.0));
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			xi[i][j] = 0.8 * exp(-(i * h - 6.0) * (i * h - 6.0) - (j * h - 6.0) * (j * h - 6.0)) +
				0.2 * exp(-(i * h - 2.0) * (i * h - 2.0) - (j * h - 2.0) * (j * h - 2.0));
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
	vector<vector<vector<vector<complex<double>>>>> a(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	// выделение памяти под 3-х мерный "квадратный" комплексный массив
	vector<vector<vector<complex<double>>>> overline_a(N + 1, vector<vector<complex<double>>>(N + 1,
		vector<complex<double>>(N + 1, complex<double>())));

	// выделение памяти под 2-х мерный квадратный комплексный массив
	vector<vector<complex<double>>> b(N + 1, vector<complex<double>>(N + 1, complex<double>()));

	// счет индексов метода квадратур
	vector<double> index(N + 1);
	for (size_t i = 1; i < N; ++i)
	{
		if (i % 2 != 0)
		{
			index[i] = 4.0 / 3;
		}
		else
		{
			index[i] = 2.0 / 3;
		}
	}
	index[0] = 1.0 / 3;
	index[N] = 1.0 / 3;

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
	Lasting("Time calculation of basic matrices", timeStart, timeFinish);
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
	Lasting("Download time major arrays", timeStart, timeFinish);
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
	Lasting("The computation time of the source function", timeStart, timeFinish);
	timeStart = clock();

	// для нахождения u^(1) составляем СЛАУ основная матрица * u^(1) = правой части
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEequation[ii]

	const size_t N_squared = (N + 1) * (N + 1);
	vector<complex<double>> rightPartEquation(N_squared, complex<double>());
	vector<complex<double>> numbered_u(N_squared);
	vector<vector<complex<double>>> substantiveMatrix(N_squared, vector<complex<double>>(N_squared, complex<double>()));
	vector<complex<double>> overline_u(N + 1, complex<double>());

	//счет основной матрицы
	size_t ii, jj;
	complex<double> sumOfTheCoefficients;
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			sumOfTheCoefficients = complex<double>();
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					jj = p * (N + 1) + q;
					if ((i != p) || (q != j))
					{
						substantiveMatrix[ii][jj] += a[i][j][p][q] * xi[p][q];
						sumOfTheCoefficients += a[i][j][p][q];
					}
				}
			}
			substantiveMatrix[ii][ii] += 1.0;
			substantiveMatrix[ii][ii] -= sumOfTheCoefficients * xi[i][j];
			substantiveMatrix[ii][ii] += b[i][j] * xi[i][j];
		}
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
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
		Lasting("Finding the acoustic pressure in R", timeStart, timeFinish);
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
		Lasting("Finding the acoustic pressure in X", timeStart, timeFinish);
		timeStart = clock();
	}
	file_overline_u.close();

	timeFinish = clock();
	Lasting("The total time of the program", timeBegin, timeFinish);
}
