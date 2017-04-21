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
	vector<vector<vector<complex<float>>>> Source_R(source.numberSource, vector<vector<complex<float>>>(N + 1, vector<complex<float>>(N + 1, complex<float>()));
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



	return 0;
}

