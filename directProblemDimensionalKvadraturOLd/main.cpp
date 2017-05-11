#include "stdafx.h"
#include "functdirectProblem2D.h"

using namespace std;

int main()
{
	int N = numberPartitionPoints_N;
	double h = (double)domainInHomogeneity_R / numberPartitionPoints_N;

	// выделение пам€ти дл€ точного решени€ xi и акустического пол€ u
	double **xi;
	complex<double> ** u;
	xi = new double *[N + 1];
	u = new complex<double> *[N + 1];
	for (int i = 0; i <= N; i++)
	{
		xi[i] = new double[N + 1];
		u[i] = new complex<double>[N + 1];
	}
	// задание точного решени€ \xi
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			xi[i][j] = exp(-(i*h - 6.0)*(i*h - 6.0) - (j*h - 6.0)*(j*h - 6.0));
		}
	}
	//печатаем точное решение в файл
	ofstream f_xi("exact_xi.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_xi << fixed << setprecision(7) << xi[i][j] << " ";
		}
	}
	f_xi.close();
	// выдел€ем пам€ть под основные матрицы
	complex<double> ****a;
	a = new complex<double> ***[N + 1];
	for (int i = 0; i <= N; i++)
	{
		a[i] = new complex<double>**[N + 1];
		for (int j = 0; j <= N; j++)
		{
			a[i][j] = new complex<double>*[N + 1];
			{
				for (int p = 0; p <= N; p++)
				{
					a[i][j][p] = new complex<double>[N + 1];
				}
			}
		}
	}
	complex<double> **b;
	b = new complex<double> *[N + 1];
	for (int i = 0; i <= N; i++)
	{
		b[i] = new complex<double>[N + 1];
	}
	complex<double> ***overline_a;
	overline_a = new complex<double> **[N + 1];
	for (int i = 0; i <= N; i++)
	{
		overline_a[i] = new complex<double>*[N + 1];
		for (int j = 0; j <= N; j++)
		{
			overline_a[i][j] = new complex<double>[N + 1];
		}
	}
	// 
	// Ќачало вычислений основных матриц
	//
	// начало счета времени
	clock_t timeStart, timeFinish;
	timeStart = clock();
	// дл€ индексов метода квадратур
	double *index;
	index = new double[N + 1];
	for (int i = 1; i<N; i++)
	{
		if (i % 2 != 0)
		{
			index[i] = (double)4 / 3;
		}
		else
		{
			index[i] = (double)2 / 3;
		}
	}
	index[0] = (double)1 / 3;
	index[N] = (double)1 / 3;
	// нахождение массива a
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int q = 0; q<N; q++)
				{
					if ((i != p) || (q != j))
					{
						a[i][j][p][q] = index[p] * index[q];
						a[i][j][p][q] = a[i][j][p][q] * G(i*h, j*h, p*h, q*h);
						a[i][j][p][q] = a[i][j][p][q] * h*h;
					}
				}
			}
		}
	}
	// нахождение массива overline_a
	for (int j = 0; j <= N; j++)
	{
		for (int p = 0; p<N; p++)
		{
			for (int q = 0; q<N; q++)
			{
				overline_a[j][p][q] = index[p] * index[q];
				overline_a[j][p][q] = overline_a[j][p][q] * G(11, j*h, p*h, q*h);
				overline_a[j][p][q] = overline_a[j][p][q] * h*h;
			}
		}
	}
	// нахождение массива b
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			b[i][j] = 0.0;
			for (int p = 0; p<N; p++)
			{
				for (int q = 0; q<N; q++)
				{
					if (i != p)
					{
						b[i][j] += G(i*h, j*h, p*h, q*h);
					}
				}
				b[i][j] *= h*h;
			}
		}
	}

	//печатаем врем€ работы
	timeFinish = clock();
	double d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Time calculation of basic matrices " << d << endl;
	timeStart = clock();

	// записываем массивы в файлы
	ofstream f_a("matrix_a.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int q = 0; q<N; q++)
				{
					f_a << fixed << setprecision(12) << a[i][j][p][q] << " ";


				}
			}
		}
	}
	f_a.close();
	ofstream f_overline_a("matrix_overline_a.txt");
	for (int j = 0; j <= N; j++)
	{
		for (int p = 0; p<N; p++)
		{
			for (int q = 0; q<N; q++)
			{
				f_overline_a << fixed << setprecision(12) << overline_a[j][p][q] << " ";
			}
		}
	}
	f_overline_a.close();
	ofstream f_b("matrix_b.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_b << fixed << setprecision(12) << b[i][j] << " ";
		}
	}
	f_b.close();
	//печатаем врем€ работы
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Download time major arrays " << d << endl;
	timeStart = clock();
	// счет функции источника в R и X
	// первый источник
	ofstream f_Source_01("f_Source_01.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_Source_01 << fixed << setprecision(12) << f_01(i*h, j*h) << " ";
		}
	}
	f_Source_01.close();
	ofstream f_Source_02("f_Source_02.txt");
	for (int j = 0; j <= N; j++)
	{
		f_Source_02 << fixed << setprecision(12) << f_01(11, j*h) << " ";
	}
	f_Source_02.close();

	// второй источник
	ofstream f_Source_03("f_Source_03.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_Source_03 << fixed << setprecision(12) << f_02(i*h, j*h) << " ";
		}
	}
	f_Source_03.close();
	ofstream f_Source_04("f_Source_04.txt");
	for (int j = 0; j <= N; j++)
	{
		f_Source_04 << fixed << setprecision(12) << f_02(11, j*h) << " ";
	}
	f_Source_04.close();

	// третий источник
	ofstream f_Source_05("f_Source_05.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_Source_05 << fixed << setprecision(12) << f_03(i*h, j*h) << " ";
		}
	}
	f_Source_05.close();
	ofstream f_Source_06("f_Source_06.txt");
	for (int j = 0; j <= N; j++)
	{
		f_Source_06 << fixed << setprecision(12) << f_03(11, j*h) << " ";
	}
	f_Source_06.close();

	// четвертый источник
	ofstream f_Source_07("f_Source_07.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_Source_07 << fixed << setprecision(12) << f_04(i*h, j*h) << " ";
		}
	}
	f_Source_07.close();
	ofstream f_Source_08("f_Source_08.txt");
	for (int j = 0; j <= N; j++)
	{
		f_Source_08 << fixed << setprecision(12) << f_04(11, j*h) << " ";
	}
	f_Source_08.close();

	// п€тый источник
	ofstream f_Source_09("f_Source_09.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_Source_09 << fixed << setprecision(12) << f_05(i*h, j*h) << " ";
		}
	}
	f_Source_09.close();
	ofstream f_Source_10("f_Source_10.txt");
	for (int j = 0; j <= N; j++)
	{
		f_Source_10 << fixed << setprecision(12) << f_05(11, j*h) << " ";
	}
	f_Source_10.close();

	//печатаем врем€ работы
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "The computation time of the source function " << d << endl;
	timeStart = clock();
	// дл€ нахождени€ u^(1) составл€ем —Ћј” основна€ матрица * u^(1) = правой части
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEequation[ii]
	int N_squared;
	N_squared = (N + 1)*(N + 1);
	complex<double> *numbered_u, *rightPartEquation;
	rightPartEquation = new complex<double>[N_squared];
	numbered_u = new complex<double>[N_squared];
	//пам€ть под основную матрицу
	complex<double> **substantiveMatrix;
	substantiveMatrix = new complex<double> *[N_squared];
	for (int i = 0; i<N_squared; i++)
	{
		substantiveMatrix[i] = new complex<double>[N_squared];
	}
	for (int i = 0; i<N_squared; i++)
	{
		for (int j = 0; j<N_squared; j++)
		{
			substantiveMatrix[i][j] = 0.0;
		}
	}
	// выделение пам€ти дл€ overline_u_0
	complex<double> *overline_u;
	overline_u = new complex<double>[N + 1];
	//
	//счет основной матрицы
	//
	int ii, jj;
	//
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			complex<double> sumOfTheCoefficients = 0.0;
			for (int p = 0; p<N; p++)
			{
				for (int q = 0; q<N; q++)
				{
					jj = p*(N + 1) + q;
					if ((i != p) || (q != j))
					{
						substantiveMatrix[ii][jj] += a[i][j][p][q] * xi[p][q];
						sumOfTheCoefficients += a[i][j][p][q];
					}
				}
			}
			substantiveMatrix[ii][ii] += 1.0;
			substantiveMatrix[ii][ii] -= sumOfTheCoefficients*xi[i][j];
			substantiveMatrix[ii][ii] += b[i][j] * xi[i][j];
		}
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "The computation time of the matrix inside the squared " << d << endl;
	timeStart = clock();
	////////////////////////////////////////////////////////
	///ƒл€ первого источника
	///////////////////////////////////////////////////////
	// нахождение правой части
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			rightPartEquation[ii] = f_01(i*h, j*h);
		}
	}
	// нахождение u^{(1)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// ќбратна€ перенумераци€
	//
	for (int l = 0; l<N_squared; l++)
	{
		int coordinate_x = l / (N + 1);
		int coordinate_y = l % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[l];
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 1 source " << d << endl;
	timeStart = clock();
	//
	// находим overline_u_0
	//
	for (int i = 0; i <= N; i++)
	{
		overline_u[i] = f_01(11, i*h);
	}
	for (int j = 0; j <= N; j++)
	{
		for (int p = 0; p<N; p++)
		{
			for (int q = 0; q<N; q++)
			{
				overline_u[j] -= overline_a[j][p][q] * xi[p][q] * u[p][q];
			}
		}
	}
	// печать overline_u_0^(1) в файл в одну строчку
	ofstream f_overline_u_1("matrix_overline_u_1.txt");
	for (int j = 0; j <= N; j++)
	{
		f_overline_u_1 << fixed << setprecision(12) << overline_u[j] << " ";
	}
	f_overline_u_1.close();
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 1 source " << d << endl;
	timeStart = clock();



	////////////////////////////////////////////////////////
	///ƒл€ второго источника
	///////////////////////////////////////////////////////
	// нахождение правой части
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			rightPartEquation[ii] = f_02(i*h, j*h);
		}
	}
	// нахождение u^{(2)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// ќбратна€ перенумераци€
	//
	for (int l = 0; l<N_squared; l++)
	{
		int coordinate_x = l / (N + 1);
		int coordinate_y = l % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[l];
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 2 source " << d << endl;
	timeStart = clock();
	//
	// находим overline_u_0
	//
	for (int i = 0; i <= N; i++)
	{
		overline_u[i] = f_02(11, i*h);
	}
	for (int j = 0; j <= N; j++)
	{
		for (int p = 0; p<N; p++)
		{
			for (int q = 0; q<N; q++)
			{
				overline_u[j] -= overline_a[j][p][q] * xi[p][q] * u[p][q];
			}
		}
	}
	// печать overline_u_0^(2) в файл в одну строчку
	ofstream f_overline_u_2("matrix_overline_u_2.txt");
	for (int j = 0; j <= N; j++)
	{
		f_overline_u_2 << fixed << setprecision(12) << overline_u[j] << " ";
	}
	f_overline_u_2.close();
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 2 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///ƒл€ третьего источника
	///////////////////////////////////////////////////////
	// нахождение правой части
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			rightPartEquation[ii] = f_03(i*h, j*h);
		}
	}
	// нахождение u^{(3)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// ќбратна€ перенумераци€
	//
	for (int l = 0; l<N_squared; l++)
	{
		int coordinate_x = l / (N + 1);
		int coordinate_y = l % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[l];
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 3 source " << d << endl;
	timeStart = clock();
	//
	// находим overline_u_0
	//
	for (int i = 0; i <= N; i++)
	{
		overline_u[i] = f_03(11, i*h);
	}
	for (int j = 0; j <= N; j++)
	{
		for (int p = 0; p<N; p++)
		{
			for (int q = 0; q<N; q++)
			{
				overline_u[j] -= overline_a[j][p][q] * xi[p][q] * u[p][q];
			}
		}
	}
	// печать overline_u_0^(3) в файл в одну строчку
	ofstream f_overline_u_3("matrix_overline_u_3.txt");
	for (int j = 0; j <= N; j++)
	{
		f_overline_u_3 << fixed << setprecision(12) << overline_u[j] << " ";
	}
	f_overline_u_3.close();
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 3 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///ƒл€ четвЄртого источника
	///////////////////////////////////////////////////////
	// нахождение правой части
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			rightPartEquation[ii] = f_04(i*h, j*h);
		}
	}
	// нахождение u^{(4)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// ќбратна€ перенумераци€
	//
	for (int l = 0; l<N_squared; l++)
	{
		int coordinate_x = l / (N + 1);
		int coordinate_y = l % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[l];
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 4 source " << d << endl;
	timeStart = clock();
	//
	// находим overline_u_0
	//
	for (int i = 0; i <= N; i++)
	{
		overline_u[i] = f_04(11, i*h);
	}
	for (int j = 0; j <= N; j++)
	{
		for (int p = 0; p<N; p++)
		{
			for (int q = 0; q<N; q++)
			{
				overline_u[j] -= overline_a[j][p][q] * xi[p][q] * u[p][q];
			}
		}
	}
	// печать overline_u_0^(4) в файл в одну строчку
	ofstream f_overline_u_4("matrix_overline_u_4.txt");
	for (int j = 0; j <= N; j++)
	{
		f_overline_u_4 << fixed << setprecision(12) << overline_u[j] << " ";
	}
	f_overline_u_4.close();
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 4 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///ƒл€ п€того источника
	///////////////////////////////////////////////////////
	// нахождение правой части
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			rightPartEquation[ii] = f_05(i*h, j*h);
		}
	}
	// нахождение u^{(5)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// ќбратна€ перенумераци€
	//
	for (int l = 0; l<N_squared; l++)
	{
		int coordinate_x = l / (N + 1);
		int coordinate_y = l % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[l];
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 5 source " << d << endl;
	timeStart = clock();
	//
	// находим overline_u_0
	//
	for (int i = 0; i <= N; i++)
	{
		overline_u[i] = f_05(11, i*h);
	}
	for (int j = 0; j <= N; j++)
	{
		for (int p = 0; p<N; p++)
		{
			for (int q = 0; q<N; q++)
			{
				overline_u[j] -= overline_a[j][p][q] * xi[p][q] * u[p][q];
			}
		}
	}
	// печать overline_u_0^(5) в файл в одну строчку
	ofstream f_overline_u_5("matrix_overline_u_5.txt");
	for (int j = 0; j <= N; j++)
	{
		f_overline_u_5 << fixed << setprecision(12) << overline_u[j] << " ";
	}
	f_overline_u_5.close();
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 5 source " << d << endl;
	timeStart = clock();



	// освобождение пам€ти	
	for (int i = 0; i <= N; i++)
	{
		delete[] xi[i];
		delete[] u[i];
	}
	delete[] xi;
	delete[] u;
	// освобождаем основные массивы
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p <= N; p++)
			{
				delete[] a[i][j][p];
			}
			delete[] a[i][j];
		}
		delete[] a[i];
	}
	delete[] a;

	for (int i = 0; i <= N; i++)
	{
		delete[] b[i];
	}
	delete[] b;
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			delete[] overline_a[i][j];
		}
		delete[] overline_a[i];
	}
	delete[] overline_a;

	delete[] numbered_u;
	delete[] rightPartEquation;
	for (int i = 0; i<N_squared; i++)
	{
		delete[] substantiveMatrix[i];
	}
	delete[] substantiveMatrix;

	delete[] overline_u;

	return 0;
}

