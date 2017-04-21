#include "stdafx.h"
#include "basicFunctions.h"
#include "sourceFunction.h"
#include "matrix_utils.h"

using namespace std;

int main()
{
	const size_t N = NUMBER_PARTITION_POINTS;
	const float h = (float)DOMAIN_IN_HOMOGENEITY / N;

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
			index[i] = 1.333333f;
		}
		else
		{
			index[i] = 0.666667f;
		}
	}
	index[0] = 0.333333f;
	index[N] = 0.333333f;

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
	// первый источник
	{
		ofstream fileSource01("Source_01.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource01 << fixed << setprecision(6) << SourceFunction(source01, i * h, j * h) << " ";
			}
		}
		fileSource01.close();
		ofstream fileSource02("Source_02.txt");
		for (size_t j = 0; j <= N; ++j)
		{
			fileSource02 << fixed << setprecision(6) << SourceFunction(source01, receiver, j * h) << " ";
		}
		fileSource02.close();
	}

	// второй источник
	{
		ofstream fileSource03("Source_03.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource03 << fixed << setprecision(6) << SourceFunction(source02, i * h, j * h) << " ";
			}
		}
		fileSource03.close();
		ofstream fileSource04("Source_04.txt");
		for (size_t j = 0; j <= N; ++j)
		{
			fileSource04 << fixed << setprecision(6) << SourceFunction(source02, receiver, j * h) << " ";
		}
		fileSource04.close();
	}

	// третий источник
	{
		ofstream fileSource05("Source_05.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource05 << fixed << setprecision(6) << SourceFunction(source03, i * h, j * h) << " ";
			}
		}
		fileSource05.close();
		ofstream fileSource06("Source_06.txt");
		for (size_t j = 0; j <= N; ++j)
		{
			fileSource06 << fixed << setprecision(6) << SourceFunction(source03, receiver, j * h) << " ";
		}
		fileSource06.close();
	}

	// четвертый источник
	{
		ofstream fileSource07("Source_07.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource07 << fixed << setprecision(6) << SourceFunction(source04, i * h, j * h) << " ";
			}
		}
		fileSource07.close();
		ofstream fileSource08("Source_08.txt");
		for (size_t j = 0; j <= N; ++j)
		{
			fileSource08 << fixed << setprecision(6) << SourceFunction(source04, receiver, j * h) << " ";
		}
		fileSource08.close();
	}

	// пятый источник
	{
		ofstream fileSource09("Source_09.txt");
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource09 << fixed << setprecision(6) << SourceFunction(source05, i * h, j * h) << " ";
			}
		}
		fileSource09.close();
		ofstream fileSource10("Source_10.txt");
		for (size_t j = 0; j <= N; ++j)
		{
			fileSource10 << fixed << setprecision(6) << SourceFunction(source05, receiver, j * h) << " ";
		}
		fileSource10.close();
	}

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

	//
	//счет основной матрицы
	//
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

	////////////////////////////////////////////////////////
	///Для первого источника
	////////////////////////////////////////////////////////
	// нахождение правой части
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			rightPartEquation[ii] = SourceFunction(source01, i * h, j * h);
		}
	}

	// нахождение u^{(1)}
	SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);
	//
	// Обратная перенумерация
	//
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
	cout << "Finding the acoustic pressure in R for 1 source " << d << endl;
	timeStart = clock();
	//
	// находим overline_u_0
	//
	for (size_t i = 0; i <= N; ++i)
	{
		overline_u[i] = SourceFunction(source01, receiver, i * h);
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
	// печать overline_u_0^(1) в файл в одну строчку
	ofstream file_overline_u_1("matrix_overline_u_1.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		file_overline_u_1 << fixed << setprecision(6) << overline_u[j] << " ";
	}
	file_overline_u_1.close();
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 1 source " << d << endl;
	timeStart = clock();
	
	////////////////////////////////////////////////////////
	///Для второго источника
	///////////////////////////////////////////////////////
	// нахождение правой части
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			rightPartEquation[ii] = SourceFunction(source02, i * h, j * h);
		}
	}

	// нахождение u^{(2)}
	SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);
	//
	// Обратная перенумерация
	//
	for (size_t i = 0; i < N_squared; ++i)
	{
		coordinate_x = i / (N + 1);
		coordinate_y = i % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[i];
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 2 source " << d << endl;
	timeStart = clock();
	//
	// находим overline_u_0
	//
	for (size_t i = 0; i <= N; ++i)
	{
		overline_u[i] = SourceFunction(source02, receiver, i * h);
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

	// печать overline_u_0^(2) в файл в одну строчку
	ofstream file_overline_u_2("matrix_overline_u_2.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		file_overline_u_2 << fixed << setprecision(6) << overline_u[j] << " ";
	}
	file_overline_u_2.close();
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 2 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///Для третьего источника
	///////////////////////////////////////////////////////
	// нахождение правой части
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			rightPartEquation[ii] = SourceFunction(source03, i * h, j * h);
		}
	}

	// нахождение u^{(3)}
	SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);
	//
	// Обратная перенумерация
	//
	for (size_t i = 0; i < N_squared; ++i)
	{
		coordinate_x = i / (N + 1);
		coordinate_y = i % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[i];
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 3 source " << d << endl;
	timeStart = clock();
	//
	// находим overline_u_0
	//
	for (size_t i = 0; i <= N; ++i)
	{
		overline_u[i] = SourceFunction(source03, receiver, i * h);
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

	// печать overline_u_0^(3) в файл в одну строчку
	ofstream file_overline_u_3("matrix_overline_u_3.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		file_overline_u_3 << fixed << setprecision(6) << overline_u[j] << " ";
	}
	file_overline_u_3.close();
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 3 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///Для четвёртого источника
	///////////////////////////////////////////////////////
	// нахождение правой части
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			rightPartEquation[ii] = SourceFunction(source04, i * h, j * h);
		}
	}

	// нахождение u^{(4)}
	SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);
	//
	// Обратная перенумерация
	//
	for (size_t i = 0; i < N_squared; ++i)
	{
		coordinate_x = i / (N + 1);
		coordinate_y = i % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[i];
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 4 source " << d << endl;
	timeStart = clock();
	//
	// находим overline_u_0
	//
	for (size_t i = 0; i <= N; ++i)
	{
		overline_u[i] = SourceFunction(source04, receiver, i * h);
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
	// печать overline_u_0^(4) в файл в одну строчку
	ofstream file_overline_u_4("matrix_overline_u_4.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		file_overline_u_4 << fixed << setprecision(6) << overline_u[j] << " ";
	}
	file_overline_u_4.close();
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 4 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///Для пятого источника
	///////////////////////////////////////////////////////
	// нахождение правой части
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			rightPartEquation[ii] = SourceFunction(source05, i * h, j * h);
		}
	}

	// нахождение u^{(5)}
	SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);
	//
	// Обратная перенумерация
	//
	for (size_t i = 0; i < N_squared; ++i)
	{
		coordinate_x = i / (N + 1);
		coordinate_y = i % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[i];
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 5 source " << d << endl;
	timeStart = clock();
	//
	// находим overline_u_0
	//
	for (size_t i = 0; i <= N; ++i)
	{
		overline_u[i] = SourceFunction(source05, receiver, i * h);
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
	// печать overline_u_0^(5) в файл в одну строчку
	ofstream file_overline_u_5("matrix_overline_u_5.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		file_overline_u_5 << fixed << setprecision(6) << overline_u[j] << " ";
	}
	file_overline_u_5.close();
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 5 source " << d << endl;

	d = (float)(timeFinish - timeBegin) / CLOCKS_PER_SEC;
	cout << "The total time of the program " << d << endl;
	
	return 0;
}
