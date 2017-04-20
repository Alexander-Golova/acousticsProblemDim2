#include "stdafx.h"
#include "basicFunctions.h"
#include"Sources.h"
#include "taskData.h"
#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"

using namespace std;

int main()
{
	const size_t N = NUMBER_PARTITION_POINTS;
	const float h = (float)DOMAIN_IN_HOMOGENEITY / N;

	const Source source;

	// выделение пам€ти дл€ акустического пол€ u
	vector<vector<complex<float>>> u(N + 1, vector<complex<float>>(N + 1, complex<float>()));

	// задание точного решени€ \xi
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
	// Ќачало вычислений основных матриц
	//
	// начало счета времени
	clock_t timeStart, timeFinish, timeBegin;
	timeBegin = clock();
	timeStart = clock();


	// выделение пам€ти под 4-х мерный "квадратный" комплексный массив
	vector<vector<vector<vector<complex<float>>>>> a(N + 1,
		vector<vector<vector<complex<float>>>>(N + 1, vector<vector<complex<float>>>(N + 1,
			vector<complex<float>>(N + 1, complex<float>()))));

	// выделение пам€ти под 3-х мерный "квадратный" комплексный массив
	vector<vector<vector<complex<float>>>> overline_a(N + 1, vector<vector<complex<float>>>(N + 1,
		vector<complex<float>>(N + 1, complex<float>())));

	// выделение пам€ти под 2-х мерный квадратный комплексный массив
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

	//печатаем врем€ работы
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

	//печатаем врем€ работы
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Download time major arrays " << d << endl;
	timeStart = clock();

	// счет функции источника в R и X

	for (size_t count = 0; count < source.numberSource; ++count)
	{
		ofstream fileSource01("Source_01.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource01 << fixed << setprecision(6) << source.Function(source.node[count], i * h, j * h) << " ";
			}
		}
		fileSource01.close();
		ofstream fileSource02("Source_02.txt");
		for (size_t j = 0; j <= N; ++j)
		{
			fileSource02 << fixed << setprecision(6) << source.Function(source.node[count], receiver, j * h) << " ";
		}
		fileSource02.close();
	}

	return 0;
}

