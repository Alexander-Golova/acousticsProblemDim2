#include "stdafx.h"
#include "functdirectProblem2D.h"

using namespace std;

int main()
{
	size_t N = numberPartitionPosize_ts_N;
	float h = (float)domainInHomogeneity_R / numberPartitionPosize_ts_N;

	// ��������� ������ ��� ������� ������� xi � ������������� ���� u
	float **xi;
	complex<float> ** u;
	xi = new float *[N + 1];
	u = new complex<float> *[N + 1];
	for (size_t i = 0; i <= N; ++i)
	{
		xi[i] = new float[N + 1];
		u[i] = new complex<float>[N + 1];
	}
	// ������� ������� ������� \xi
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			xi[i][j] = exp(-(i * h - 6.0f) * (i * h - 6.0f) - (j * h - 6.0f) * (j * h - 6.0f));
		}
	}

	//�������� ������ ������� � ����
	ofstream f_xi("exact_xi.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_xi << fixed << setprecision(7) << xi[i][j] << " ";
		}
	}
	f_xi.close();

	// �������� ������ ��� �������� �������
	complex<float> ****a;
	a = new complex<float> ***[N + 1];
	for (size_t i = 0; i <= N; ++i)
	{
		a[i] = new complex<float>**[N + 1];
		for (size_t j = 0; j <= N; ++j)
		{
			a[i][j] = new complex<float>*[N + 1];
			{
				for (size_t p = 0; p <= N; ++p)
				{
					a[i][j][p] = new complex<float>[N + 1];
				}
			}
		}
	}
	complex<float> **b;
	b = new complex<float> *[N + 1];
	for (size_t i = 0; i <= N; ++i)
	{
		b[i] = new complex<float>[N + 1];
	}
	complex<float> ***overline_a;
	overline_a = new complex<float> **[N + 1];
	for (size_t i = 0; i <= N; ++i)
	{
		overline_a[i] = new complex<float>*[N + 1];
		for (size_t j = 0; j <= N; ++j)
		{
			overline_a[i][j] = new complex<float>[N + 1];
		}
	}

	// ������ ���������� �������� ������
	//
	// ������ ����� �������
	clock_t timeStart, timeFinish;
	timeStart = clock();
	// ��� �������� ������ ���������
	float *index;
	index = new float[N + 1];
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
	// ���������� ������� a
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
	// ���������� ������� overline_a
	for (size_t j = 0; j <= N; ++j)
	{
		for (size_t p = 0; p < N; ++p)
		{
			for (size_t q = 0; q < N; ++q)
			{
				overline_a[j][p][q] = index[p] * index[q];
				overline_a[j][p][q] = overline_a[j][p][q] * G(11.0f, j * h, p * h, q * h);
				overline_a[j][p][q] = overline_a[j][p][q] * h * h;
			}
		}
	}
	// ���������� ������� b
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			b[i][j] = { 0.0f, 0.0f };
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

	//�������� ����� ������
	timeFinish = clock();
	float d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Time calculation of basic matrices " << d << endl;
	timeStart = clock();

	// ���������� ������� � �����
	ofstream f_a("matrix_a.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t q = 0; q < N; ++q)
				{
					f_a << fixed << setprecision(12) << a[i][j][p][q] << " ";


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
				f_overline_a << fixed << setprecision(12) << overline_a[j][p][q] << " ";
			}
		}
	}
	f_overline_a.close();
	ofstream f_b("matrix_b.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_b << fixed << setprecision(12) << b[i][j] << " ";
		}
	}
	f_b.close();
	//�������� ����� ������
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Download time major arrays " << d << endl;
	timeStart = clock();

	// ���� ������� ��������� � R � X
	// ������ ��������
	ofstream f_Source_01("f_Source_01.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_Source_01 << fixed << setprecision(12) << f_01(i * h, j * h) << " ";
		}
	}
	f_Source_01.close();
	ofstream f_Source_02("f_Source_02.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		f_Source_02 << fixed << setprecision(12) << f_01(11.0f, j * h) << " ";
	}
	f_Source_02.close();

	// ������ ��������
	ofstream f_Source_03("f_Source_03.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_Source_03 << fixed << setprecision(12) << f_02(i * h, j * h) << " ";
		}
	}
	f_Source_03.close();
	ofstream f_Source_04("f_Source_04.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		f_Source_04 << fixed << setprecision(12) << f_02(11.0f, j * h) << " ";
	}
	f_Source_04.close();

	// ������ ��������
	ofstream f_Source_05("f_Source_05.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_Source_05 << fixed << setprecision(12) << f_03(i * h, j * h) << " ";
		}
	}
	f_Source_05.close();
	ofstream f_Source_06("f_Source_06.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		f_Source_06 << fixed << setprecision(12) << f_03(11.0f, j * h) << " ";
	}
	f_Source_06.close();

	// ��������� ��������
	ofstream f_Source_07("f_Source_07.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_Source_07 << fixed << setprecision(12) << f_04(i * h, j * h) << " ";
		}
	}
	f_Source_07.close();
	ofstream f_Source_08("f_Source_08.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		f_Source_08 << fixed << setprecision(12) << f_04(11.0f, j * h) << " ";
	}
	f_Source_08.close();

	// ����� ��������
	ofstream f_Source_09("f_Source_09.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_Source_09 << fixed << setprecision(12) << f_05(i * h, j * h) << " ";
		}
	}
	f_Source_09.close();
	ofstream f_Source_10("f_Source_10.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		f_Source_10 << fixed << setprecision(12) << f_05(11.0f, j * h) << " ";
	}
	f_Source_10.close();

	//�������� ����� ������
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "The computation time of the source function " << d << endl;
	timeStart = clock();

	// ��� ���������� u^(1) ���������� ���� �������� ������� * u^(1) = ������ �����
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEequation[ii]
	size_t N_squared = (N + 1) * (N + 1);
	complex<float> *numbered_u, *rightPartEquation;
	rightPartEquation = new complex<float>[N_squared];
	numbered_u = new complex<float>[N_squared];
	//������ ��� �������� �������
	complex<float> **substantiveMatrix;
	substantiveMatrix = new complex<float> *[N_squared];
	for (size_t i = 0; i < N_squared; ++i)
	{
		substantiveMatrix[i] = new complex<float>[N_squared];
	}
	for (size_t i = 0; i < N_squared; ++i)
	{
		for (size_t j = 0; j < N_squared; ++j)
		{
			substantiveMatrix[i][j] = { 0.0f, 0.0f };
		}
	}

	// ��������� ������ ��� overline_u_0
	complex<float> *overline_u;
	overline_u = new complex<float>[N + 1];
	//
	//���� �������� �������
	//
	size_t ii, jj;
	//
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			complex<float> sumOfTheCoefficients = 0.0f;
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
			substantiveMatrix[ii][ii] += 1.0f;
			substantiveMatrix[ii][ii] -= sumOfTheCoefficients * xi[i][j];
			substantiveMatrix[ii][ii] += b[i][j] * xi[i][j];
		}
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "The computation time of the matrix inside the squared " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///��� ������� ���������
	///////////////////////////////////////////////////////
	// ���������� ������ �����
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			rightPartEquation[ii] = f_01(i * h, j * h);
		}
	}
	// ���������� u^{(1)}
	SolveSlauGaussa(substantiveMatrix, (int)N_squared, rightPartEquation, numbered_u);
	//
	// �������� �������������
	//
	
	for (size_t l = 0; l < N_squared; ++l)
	{
		size_t coordinate_x = l / (N + 1);
		size_t coordinate_y = l % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[l];
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 1 source " << d << endl;
	timeStart = clock();
	//
	// ������� overline_u_0
	//
	for (size_t i = 0; i <= N; ++i)
	{
		overline_u[i] = f_01(11.0f, i * h);
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
	// ������ overline_u_0^(1) � ���� � ���� �������
	ofstream f_overline_u_1("matrix_overline_u_1.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		f_overline_u_1 << fixed << setprecision(12) << overline_u[j] << " ";
	}
	f_overline_u_1.close();
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 1 source " << d << endl;
	timeStart = clock();


	////////////////////////////////////////////////////////
	///��� ������� ���������
	///////////////////////////////////////////////////////
	// ���������� ������ �����
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			rightPartEquation[ii] = f_02(i * h, j * h);
		}
	}
	// ���������� u^{(2)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// �������� �������������
	//
	for (size_t l = 0; l < N_squared; ++l)
	{
		size_t coordinate_x = l / (N + 1);
		size_t coordinate_y = l % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[l];
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 2 source " << d << endl;
	timeStart = clock();
	//
	// ������� overline_u_0
	//
	for (size_t i = 0; i <= N; ++i)
	{
		overline_u[i] = f_02(11.0f, i * h);
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
	// ������ overline_u_0^(2) � ���� � ���� �������
	ofstream f_overline_u_2("matrix_overline_u_2.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		f_overline_u_2 << fixed << setprecision(12) << overline_u[j] << " ";
	}
	f_overline_u_2.close();
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 2 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///��� �������� ���������
	///////////////////////////////////////////////////////
	// ���������� ������ �����
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			rightPartEquation[ii] = f_03(i * h, j * h);
		}
	}
	// ���������� u^{(3)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// �������� �������������
	//
	for (size_t l = 0; l < N_squared; ++l)
	{
		size_t coordinate_x = l / (N + 1);
		size_t coordinate_y = l % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[l];
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 3 source " << d << endl;
	timeStart = clock();
	//
	// ������� overline_u_0
	//
	for (size_t i = 0; i <= N; ++i)
	{
		overline_u[i] = f_03(11.0f, i * h);
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
	// ������ overline_u_0^(3) � ���� � ���� �������
	ofstream f_overline_u_3("matrix_overline_u_3.txt");
	for (size_t j = 0; j <= N; j++)
	{
		f_overline_u_3 << fixed << setprecision(12) << overline_u[j] << " ";
	}
	f_overline_u_3.close();
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 3 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///��� ��������� ���������
	///////////////////////////////////////////////////////
	// ���������� ������ �����
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			rightPartEquation[ii] = f_04(i * h, j * h);
		}
	}
	// ���������� u^{(4)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// �������� �������������
	//
	for (size_t l = 0; l < N_squared; ++l)
	{
		size_t coordinate_x = l / (N + 1);
		size_t coordinate_y = l % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[l];
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 4 source " << d << endl;
	timeStart = clock();
	//
	// ������� overline_u_0
	//
	for (size_t i = 0; i <= N; ++i)
	{
		overline_u[i] = f_04(11.0f, i * h);
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
	// ������ overline_u_0^(4) � ���� � ���� �������
	ofstream f_overline_u_4("matrix_overline_u_4.txt");
	for (size_t j = 0; j <= N; ++j)
	{
		f_overline_u_4 << fixed << setprecision(12) << overline_u[j] << " ";
	}
	f_overline_u_4.close();
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 4 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///��� ������ ���������
	///////////////////////////////////////////////////////
	// ���������� ������ �����
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;
			rightPartEquation[ii] = f_05(i * h, j * h);
		}
	}
	// ���������� u^{(5)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// �������� �������������
	//
	for (size_t l = 0; l < N_squared; ++l)
	{
		size_t coordinate_x = l / (N + 1);
		size_t coordinate_y = l % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[l];
	}
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 5 source " << d << endl;
	timeStart = clock();
	//
	// ������� overline_u_0
	//
	for (size_t i = 0; i <= N; ++i)
	{
		overline_u[i] = f_05(11.0f, i * h);
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
	// ������ overline_u_0^(5) � ���� � ���� �������
	ofstream f_overline_u_5("matrix_overline_u_5.txt");
	for (size_t j = 0; j <= N; j++)
	{
		f_overline_u_5 << fixed << setprecision(12) << overline_u[j] << " ";
	}
	f_overline_u_5.close();
	timeFinish = clock();
	d = (float)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 5 source " << d << endl;
	timeStart = clock();


	// ������������ ������	
	for (size_t i = 0; i <= N; ++i)
	{
		delete[] xi[i];
		delete[] u[i];
	}
	delete[] xi;
	delete[] u;
	// ����������� �������� �������
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

	delete[] numbered_u;
	delete[] rightPartEquation;
	for (size_t i = 0; i < N_squared; ++i)
	{
		delete[] substantiveMatrix[i];
	}
	delete[] substantiveMatrix;

	delete[] overline_u;

	return 0;
}
