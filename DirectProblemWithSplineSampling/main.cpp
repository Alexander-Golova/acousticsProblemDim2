#include "stdafx.h"
#include "basicFunctions.h"
#include "Sources.h"
#include "taskData.h"
#include "basicArrays.h"
#include "exact_solution.h"
#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"

using namespace std;

int main()
{
	const size_t N = NUMBER_PARTITION_POINT;
	const Source source;

	// ��������� ������ ��� ������������� ���� u
	vector<vector<complex<double>>> u(NUMBER_PARTITION_POINT + 1,
		vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()));

	// ������� ������� ������� xi
	vector<vector<double>> xi(NUMBER_PARTITION_POINT + 1, vector<double>(NUMBER_PARTITION_POINT + 1, 0.0));
	GetExactSolution(xi);
	WriteSolutionFile(xi);

	// �������� ������ ��� �������� �������
	vector<vector<vector<vector<complex<double>>>>> aa_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> ab_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> ac_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> bb_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> bc_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> cc_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> aa_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> ab_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> ac_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> bb_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> bc_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> cc_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_aa_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_ab_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_ac_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_bb_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_bc_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_cc_1(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_aa_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_ab_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_ac_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_bb_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_bc_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_cc_2(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	// ������ ����� �������
	clock_t time = clock();
	clock_t timeBegin = clock();

	// ���������� �������� ������
	GetBasicArrays();
	// ������������ ������������ ������� �������� �������
	double x1, x2;
	complex<double> temp;
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					//������ ������������
					x1 = h / 3.0 + p * h;
					x2 = h / 3.0 + s * h;
					temp = -27.0 * G(i * h, j * h, x1, x2) / 9.0;
					aa_1[i][j][p][s] = temp;
					ab_1[i][j][p][s] = temp;
					ac_1[i][j][p][s] = temp;
					bb_1[i][j][p][s] = temp;
					bc_1[i][j][p][s] = temp;
					cc_1[i][j][p][s] = temp;

					temp = -27.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2) / 9.0;
					overline_aa_1[i][j][p][s] = temp;
					overline_ab_1[i][j][p][s] = temp;
					overline_ac_1[i][j][p][s] = temp;
					overline_bb_1[i][j][p][s] = temp;
					overline_bc_1[i][j][p][s] = temp;
					overline_cc_1[i][j][p][s] = temp;

					x1 = h * 0.2 + p * h;
					x2 = h * 0.2 + s * h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					aa_1[i][j][p][s] += temp * 0.04;
					ab_1[i][j][p][s] += temp * 0.04;
					ac_1[i][j][p][s] += temp * 0.12;
					bb_1[i][j][p][s] += temp * 0.04;
					bc_1[i][j][p][s] += temp * 0.12;
					cc_1[i][j][p][s] += temp * 0.36;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					overline_aa_1[i][j][p][s] += temp * 0.04;
					overline_ab_1[i][j][p][s] += temp * 0.04;
					overline_ac_1[i][j][p][s] += temp * 0.12;
					overline_bb_1[i][j][p][s] += temp * 0.04;
					overline_bc_1[i][j][p][s] += temp * 0.12;
					overline_cc_1[i][j][p][s] += temp * 0.36;

					x1 = h * 0.6 + p * h;
					x2 = h * 0.2 + s * h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					aa_1[i][j][p][s] += temp * 0.36;
					ab_1[i][j][p][s] += temp * 0.12;
					ac_1[i][j][p][s] += temp * 0.12;
					bb_1[i][j][p][s] += temp * 0.04;
					bc_1[i][j][p][s] += temp * 0.04;
					cc_1[i][j][p][s] += temp * 0.04;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					overline_aa_1[i][j][p][s] += temp * 0.36;
					overline_ab_1[i][j][p][s] += temp * 0.12;
					overline_ac_1[i][j][p][s] += temp * 0.12;
					overline_bb_1[i][j][p][s] += temp * 0.04;
					overline_bc_1[i][j][p][s] += temp * 0.04;
					overline_cc_1[i][j][p][s] += temp * 0.04;

					x1 = h * 0.2 + p * h;
					x2 = h * 0.6 + s * h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					aa_1[i][j][p][s] += temp * 0.04;
					ab_1[i][j][p][s] += temp * 0.12;
					ac_1[i][j][p][s] += temp * 0.04;
					bb_1[i][j][p][s] += temp * 0.36;
					bc_1[i][j][p][s] += temp * 0.12;
					cc_1[i][j][p][s] += temp * 0.04;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					overline_aa_1[i][j][p][s] += temp * 0.04;
					overline_ab_1[i][j][p][s] += temp * 0.12;
					overline_ac_1[i][j][p][s] += temp * 0.04;
					overline_bb_1[i][j][p][s] += temp * 0.36;
					overline_bc_1[i][j][p][s] += temp * 0.12;
					overline_cc_1[i][j][p][s] += temp * 0.04;

					temp = h * h / 96.0;
					aa_1[i][j][p][s] *= temp;
					ab_1[i][j][p][s] *= temp;
					ac_1[i][j][p][s] *= temp;
					bb_1[i][j][p][s] *= temp;
					bc_1[i][j][p][s] *= temp;
					cc_1[i][j][p][s] *= temp;

					overline_aa_1[i][j][p][s] *= temp;
					overline_ab_1[i][j][p][s] *= temp;
					overline_ac_1[i][j][p][s] *= temp;
					overline_bb_1[i][j][p][s] *= temp;
					overline_bc_1[i][j][p][s] *= temp;
					overline_cc_1[i][j][p][s] *= temp;

					//������ ������������
					x1 = -h / 3.0 + p * h + h;
					x2 = -h / 3.0 + s * h + h;
					temp = -27.0 * G(i * h, j * h, x1, x2) / 9.0;
					aa_2[i][j][p][s] = temp;
					ab_2[i][j][p][s] = temp;
					ac_2[i][j][p][s] = temp;
					bb_2[i][j][p][s] = temp;
					bc_2[i][j][p][s] = temp;
					cc_2[i][j][p][s] = temp;

					temp = -27.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2) / 9.0;
					overline_aa_2[i][j][p][s] = temp;
					overline_ab_2[i][j][p][s] = temp;
					overline_ac_2[i][j][p][s] = temp;
					overline_bb_2[i][j][p][s] = temp;
					overline_bc_2[i][j][p][s] = temp;
					overline_cc_2[i][j][p][s] = temp;

					x1 = -h * 0.2 + p * h + h;
					x2 = -h * 0.2 + s * h + h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					aa_2[i][j][p][s] += temp * 0.04;
					ab_2[i][j][p][s] += temp * 0.04;
					ac_2[i][j][p][s] += temp * 0.12;
					bb_2[i][j][p][s] += temp * 0.04;
					bc_2[i][j][p][s] += temp * 0.12;
					cc_2[i][j][p][s] += temp * 0.36;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					overline_aa_2[i][j][p][s] += temp * 0.04;
					overline_ab_2[i][j][p][s] += temp * 0.04;
					overline_ac_2[i][j][p][s] += temp * 0.12;
					overline_bb_2[i][j][p][s] += temp * 0.04;
					overline_bc_2[i][j][p][s] += temp * 0.12;
					overline_cc_2[i][j][p][s] += temp * 0.36;

					x1 = -h * 0.2 + p * h + h;
					x2 = -h * 0.6 + s * h + h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					aa_2[i][j][p][s] += temp * 0.36;
					ab_2[i][j][p][s] += temp * 0.12;
					ac_2[i][j][p][s] += temp * 0.12;
					bb_2[i][j][p][s] += temp * 0.04;
					bc_2[i][j][p][s] += temp * 0.04;
					cc_2[i][j][p][s] += temp * 0.04;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					overline_aa_2[i][j][p][s] += temp * 0.36;
					overline_ab_2[i][j][p][s] += temp * 0.12;
					overline_ac_2[i][j][p][s] += temp * 0.12;
					overline_bb_2[i][j][p][s] += temp * 0.04;
					overline_bc_2[i][j][p][s] += temp * 0.04;
					overline_cc_2[i][j][p][s] += temp * 0.04;

					x1 = -h * 0.6 + p * h + h;
					x2 = -h * 0.2 + s * h + h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					aa_2[i][j][p][s] += temp * 0.04;
					ab_2[i][j][p][s] += temp * 0.12;
					ac_2[i][j][p][s] += temp * 0.04;
					bb_2[i][j][p][s] += temp * 0.36;
					bc_2[i][j][p][s] += temp * 0.12;
					cc_2[i][j][p][s] += temp * 0.04;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					overline_aa_2[i][j][p][s] += temp * 0.04;
					overline_ab_2[i][j][p][s] += temp * 0.12;
					overline_ac_2[i][j][p][s] += temp * 0.04;
					overline_bb_2[i][j][p][s] += temp * 0.36;
					overline_bc_2[i][j][p][s] += temp * 0.12;
					overline_cc_2[i][j][p][s] += temp * 0.04;

					temp = temp = h * h / 96.0;
					aa_2[i][j][p][s] *= temp;
					ab_2[i][j][p][s] *= temp;
					ac_2[i][j][p][s] *= temp;
					bb_2[i][j][p][s] *= temp;
					bc_2[i][j][p][s] *= temp;
					cc_2[i][j][p][s] *= temp;

					overline_aa_2[i][j][p][s] *= temp;
					overline_ab_2[i][j][p][s] *= temp;
					overline_ac_2[i][j][p][s] *= temp;
					overline_bb_2[i][j][p][s] *= temp;
					overline_bc_2[i][j][p][s] *= temp;
					overline_cc_2[i][j][p][s] *= temp;
				}
			}
		}
	}

	//�������� ����� ������
	Lasting("Time calculation of basic matrices", time);

	//�������� ������������ � �����
	ofstream f_matrix("matrix.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << aa_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << ab_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << ac_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << bb_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << bc_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << cc_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << aa_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << ab_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << ac_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << bb_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << bc_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << cc_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_aa_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_ab_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_ac_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_bb_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_bc_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_cc_1[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_aa_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_ab_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_ac_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_bb_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_bc_2[i][j][p][s] << " ";
				}
			}
		}
	}
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			for (size_t p = 0; p < N; ++p)
			{
				for (size_t s = 0; s < N; ++s)
				{
					f_matrix << fixed << setprecision(6) << overline_cc_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_matrix.close();
	//�������� ����� ������
	Lasting("Time recording the basic matrix file", time);

	// ���� ������� ��������� � R � X
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

		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				fileSource << fixed << setprecision(6) << source.Function(source.node[count], RECEIVER + i * stepReceiver, j * h) << " ";
			}
		}
	}
	fileSource.close();

	//�������� ����� ������
	Lasting("The computation time of the source function", time);

	// ��� ���������� u^(1) ���������� ���� �������� ������� * u^(1) = ������ �����
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEequation[ii]

	const size_t N_squared = (N + 1) * (N + 1);
	vector<complex<double>> rightPartEquation(N_squared, complex<double>());
	vector<complex<double>> numbered_u(N_squared);
	vector<vector<complex<double>>> substantiveMatrix(N_squared, vector<complex<double>>(N_squared, complex<double>()));
	vector<vector<complex<double>>> overline_u(N + 1, vector<complex<double>>(N + 1, complex<double>()));

	//���� �������� �������
	size_t ii, jj;
	size_t auxInd;

	// ������� ��� ������ ��������
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;

			// ������� O(0; 0)
			//1 �����������
			substantiveMatrix[ii][0] = aa_1[i][j][0][0] * xi[0][0];
			substantiveMatrix[ii][0] += ab_1[i][j][0][0] * xi[0][1];
			substantiveMatrix[ii][0] += ac_1[i][j][0][0] * xi[1][0];

			//2 �����������
			// ������� A(N; 0)
			auxInd = N * (N + 1);
			//1 �����������
			substantiveMatrix[ii][auxInd] = aa_1[i][j][N][0] * xi[N][0];
			substantiveMatrix[ii][auxInd] += ab_1[i][j][N][0] * xi[N][1];
			substantiveMatrix[ii][auxInd] += ac_1[i][j][N - 1][0] * xi[N - 1][0];
			substantiveMatrix[ii][auxInd] += bc_1[i][j][N - 1][0] * xi[N - 1][1];
			substantiveMatrix[ii][auxInd] += cc_1[i][j][N - 1][0] * xi[N][0];
			//2 �����������
			substantiveMatrix[ii][auxInd] += ac_2[i][j][N - 1][0] * xi[N][1];
			substantiveMatrix[ii][auxInd] += bc_2[i][j][N - 1][0] * xi[N - 1][1];
			substantiveMatrix[ii][auxInd] += cc_2[i][j][N - 1][0] * xi[N][0];

			// ������� B(0; N)
			auxInd = N;
			//1 �����������
			substantiveMatrix[ii][auxInd] = aa_1[i][j][0][N] * xi[0][N];
			substantiveMatrix[ii][auxInd] += ab_1[i][j][0][N - 1] * xi[0][N - 1];
			substantiveMatrix[ii][auxInd] += ac_1[i][j][0][N] * xi[1][N];
			substantiveMatrix[ii][auxInd] += bb_1[i][j][0][N - 1] * xi[0][N];
			substantiveMatrix[ii][auxInd] += bc_1[i][j][0][N - 1] * xi[1][N - 1];
			//2 �����������
			substantiveMatrix[ii][auxInd] += ab_2[i][j][0][N - 1] * xi[1][N];
			substantiveMatrix[ii][auxInd] += bb_2[i][j][0][N - 1] * xi[0][N];
			substantiveMatrix[ii][auxInd] += bc_2[i][j][0][N - 1] * xi[1][N - 1];

			// ������� C(N; N)
			auxInd = N * (N + 1) + N;
			//1 �����������
			substantiveMatrix[ii][auxInd] = aa_1[i][j][N][N] * xi[N][N];
			substantiveMatrix[ii][auxInd] += ac_1[i][j][N - 1][N] * xi[N - 1][N];
			substantiveMatrix[ii][auxInd] += bb_1[i][j][N][N - 1] * xi[N][N];
			substantiveMatrix[ii][auxInd] += cc_1[i][j][N - 1][N] * xi[N][N];
			//2 �����������
			substantiveMatrix[ii][auxInd] += aa_2[i][j][N - 1][N - 1] * xi[N][N];
			substantiveMatrix[ii][auxInd] += ab_2[i][j][N - 1][N] * xi[N][N];
			substantiveMatrix[ii][auxInd] += ac_2[i][j][N - 1][N - 1] * xi[N][N - 1];
			substantiveMatrix[ii][auxInd] += bb_2[i][j][N][N - 1] * xi[N][N];
			substantiveMatrix[ii][auxInd] += cc_2[i][j][N - 1][N] * xi[N][N];
		}
	}

	// ������� ��� ������ ��������
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			ii = i * (N + 1) + j;

			// ������� OA(p; 0)
			for (size_t p = 1; p < N; ++p)
			{
				auxInd = p * (N + 1);
				//1 �����������
				substantiveMatrix[ii][auxInd] = aa_1[i][j][p][0] * xi[p][0];
				substantiveMatrix[ii][auxInd] += ab_1[i][j][p][0] * xi[p][1];
				substantiveMatrix[ii][auxInd] += ac_1[i][j][p - 1][0] * xi[p - 1][0];
				substantiveMatrix[ii][auxInd] += ac_1[i][j][p][0] * xi[p + 1][0];
				substantiveMatrix[ii][auxInd] += bc_1[i][j][p - 1][0] * xi[p - 1][1];
				substantiveMatrix[ii][auxInd] += cc_1[i][j][p - 1][0] * xi[p][0];
				//2 �����������
				substantiveMatrix[ii][auxInd] += ac_2[i][j][p - 1][0] * xi[p][1];
				substantiveMatrix[ii][auxInd] += bc_2[i][j][p - 1][0] * xi[p - 1][1];
				substantiveMatrix[ii][auxInd] += cc_2[i][j][p - 1][0] * xi[p][0];
			}

			// ������� OB(0; s)
			for (size_t s = 1; s < N; ++s)
			{
				auxInd = s;
				//1 �����������
				substantiveMatrix[ii][auxInd] = aa_1[i][j][0][s] * xi[0][s];
				substantiveMatrix[ii][auxInd] += ab_1[i][j][0][s - 1] * xi[0][s - 1];
				substantiveMatrix[ii][auxInd] += ab_1[i][j][0][s] * xi[0][s + 1];
				substantiveMatrix[ii][auxInd] += ac_1[i][j][0][s] * xi[1][s];
				substantiveMatrix[ii][auxInd] += bb_1[i][j][0][s - 1] * xi[0][s];
				substantiveMatrix[ii][auxInd] += bc_1[i][j][0][s - 1] * xi[1][s - 1];
				//2 �����������
				substantiveMatrix[ii][auxInd] += ab_2[i][j][0][s - 1] * xi[1][s];
				substantiveMatrix[ii][auxInd] += bb_2[i][j][0][s - 1] * xi[0][s];
				substantiveMatrix[ii][auxInd] += bc_2[i][j][0][s - 1] * xi[1][s - 1];
			}
			//
			// ������� BC(p; N)
			//
			for (size_t p = 1; p < N; ++p)
			{
				auxInd = p * (N + 1) + N;
				//1 �����������
				substantiveMatrix[ii][auxInd] = aa_1[i][j][p][N] * xi[p][N];
				substantiveMatrix[ii][auxInd] += ac_1[i][j][p - 1][N] * xi[p - 1][N];
				substantiveMatrix[ii][auxInd] += ac_1[i][j][p][N] * xi[p + 1][N];
				substantiveMatrix[ii][auxInd] += bb_1[i][j][p][N - 1] * xi[p][N];
				substantiveMatrix[ii][auxInd] += bc_1[i][j][p][N - 1] * xi[p + 1][N - 1];
				substantiveMatrix[ii][auxInd] += cc_1[i][j][p - 1][N] * xi[p][N];
				//2 �����������
				substantiveMatrix[ii][auxInd] += aa_2[i][j][p - 1][N - 1] * xi[p][N];
				substantiveMatrix[ii][auxInd] += ab_2[i][j][p][N - 1] * xi[p + 1][N];
				substantiveMatrix[ii][auxInd] += ab_2[i][j][p - 1][N - 1] * xi[p - 1][N];
				substantiveMatrix[ii][auxInd] += ac_2[i][j][p - 1][N - 1] * xi[p][N - 1];
				substantiveMatrix[ii][auxInd] += bb_2[i][j][p][N - 1] * xi[p][N];
				substantiveMatrix[ii][auxInd] += bc_2[i][j][p][N - 1] * xi[p + 1][N - 1];
				substantiveMatrix[ii][auxInd] += cc_2[i][j][p - 1][N] * xi[p][N];
			}

			// ������� AC(N; s)
			for (size_t s = 1; s < N; ++s)
			{
				auxInd = N * (N + 1) + s;
				//1 �����������
				substantiveMatrix[ii][auxInd] = aa_1[i][j][N][s] * xi[N][s];
				substantiveMatrix[ii][auxInd] += ab_1[i][j][N][s - 1] * xi[N][s - 1];
				substantiveMatrix[ii][auxInd] += ab_1[i][j][N][s] * xi[N][s + 1];
				substantiveMatrix[ii][auxInd] += ac_1[i][j][N - 1][s] * xi[N - 1][s];
				substantiveMatrix[ii][auxInd] += bb_1[i][j][N][s - 1] * xi[N][s];
				substantiveMatrix[ii][auxInd] += bc_1[i][j][N - 1][s] * xi[N - 1][s + 1];
				substantiveMatrix[ii][auxInd] += cc_1[i][j][N - 1][s] * xi[N][s];
				//2 �����������
				substantiveMatrix[ii][auxInd] += aa_2[i][j][N - 1][s - 1] * xi[N][s];
				substantiveMatrix[ii][auxInd] += ab_2[i][j][N - 1][s - 1] * xi[N - 1][s];
				substantiveMatrix[ii][auxInd] += ac_2[i][j][N - 1][s] * xi[N][s + 1];
				substantiveMatrix[ii][auxInd] += ac_2[i][j][N - 1][s - 1] * xi[N][s - 1];
				substantiveMatrix[ii][auxInd] += bb_2[i][j][N][s - 1] * xi[N][s];
				substantiveMatrix[ii][auxInd] += bc_2[i][j][N - 1][s] * xi[N - 1][s + 1];
				substantiveMatrix[ii][auxInd] += cc_2[i][j][N - 1][s] * xi[N][s];
			}
		}
	}

	// ������� ��� ���������� ����� ��������
	for (size_t i = 1; i < N; ++i)
	{
		for (size_t j = 1; j < N; ++j)
		{
			ii = i * (N + 1) + j;
			for (size_t p = 1; p < N; ++p)
			{
				for (size_t s = 1; s < N; ++s)
				{
					jj = p * (N + 1) + s;
					//1 �����������
					substantiveMatrix[ii][jj] = aa_1[i][j][p][s] * xi[p][s];
					substantiveMatrix[ii][jj] += ab_1[i][j][p][s - 1] * xi[p][s - 1];
					substantiveMatrix[ii][jj] += ab_1[i][j][p][s] * xi[p][s + 1];
					substantiveMatrix[ii][jj] += ac_1[i][j][p - 1][s] * xi[p - 1][s];
					substantiveMatrix[ii][jj] += ac_1[i][j][p][s] * xi[p + 1][s];
					substantiveMatrix[ii][jj] += bb_1[i][j][p][s - 1] * xi[p][s];
					substantiveMatrix[ii][jj] += bc_1[i][j][p - 1][s] * xi[p - 1][s + 1];
					substantiveMatrix[ii][jj] += bc_1[i][j][p][s - 1] * xi[p + 1][s - 1];
					substantiveMatrix[ii][jj] += cc_1[i][j][p - 1][s] * xi[p][s];
					//2 �����������
					substantiveMatrix[ii][jj] += aa_2[i][j][p - 1][s - 1] * xi[p][s];
					substantiveMatrix[ii][jj] += ab_2[i][j][p][s - 1] * xi[p + 1][s];
					substantiveMatrix[ii][jj] += ab_2[i][j][p - 1][s - 1] * xi[p - 1][s];
					substantiveMatrix[ii][jj] += ac_2[i][j][p - 1][s] * xi[p][s + 1];
					substantiveMatrix[ii][jj] += ac_2[i][j][p - 1][s - 1] * xi[p][s - 1];
					substantiveMatrix[ii][jj] += bb_2[i][j][p][s - 1] * xi[p][s];
					substantiveMatrix[ii][jj] += bc_2[i][j][p - 1][s] * xi[p - 1][s + 1];
					substantiveMatrix[ii][jj] += bc_2[i][j][p][s - 1] * xi[p + 1][s - 1];
					substantiveMatrix[ii][jj] += cc_2[i][j][p - 1][s] * xi[p][s];
				}
			}
		}
	}

	// ��������� ������� � ������� ���������
	for (size_t i = 0; i < N_squared; ++i)
	{
		substantiveMatrix[i][i] += 1.0;
	}

	Lasting("The computation time of the matrix inside the cube", time);

	// ������� ������������ ���� � ����������
	// ��� ������� ���������
	ofstream file_overline_u("matrix_overline_u.txt");
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		// ���������� ������ �����
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				ii = i * (N + 1) + j;
				rightPartEquation[ii] = source.Function(source.node[count], i * h, j * h);
			}
		}
		// ���������� u^{(count)}
		SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);

		// �������� �������������
		size_t coordinate_x;
		size_t coordinate_y;
		for (size_t i = 0; i < N_squared; ++i)
		{
			coordinate_x = i / (N + 1);
			coordinate_y = i % (N + 1);
			u[coordinate_x][coordinate_y] = numbered_u[i];
		}

		Lasting("Finding the acoustic pressure in R", time);

		// ������� overline_u_0
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				overline_u[i][j] = source.Function(source.node[count], RECEIVER + i * stepReceiver, j * h);
			}
		}
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				for (size_t p = 0; p < N; ++p)
				{
					for (size_t s = 0; s < N; ++s)
					{
						//1 �����������
						overline_u[i][j] -= overline_aa_1[i][j][p][s] * xi[p][s] * u[p][s];
						overline_u[i][j] -= overline_ab_1[i][j][p][s] * xi[p][s] * u[p][s + 1];
						overline_u[i][j] -= overline_ab_1[i][j][p][s] * xi[p][s + 1] * u[p][s];
						overline_u[i][j] -= overline_ac_1[i][j][p][s] * xi[p][s] * u[p + 1][s];
						overline_u[i][j] -= overline_ac_1[i][j][p][s] * xi[p + 1][s] * u[p][s];
						overline_u[i][j] -= overline_bb_1[i][j][p][s] * xi[p][s + 1] * u[p][s + 1];
						overline_u[i][j] -= overline_bc_1[i][j][p][s] * xi[p][s + 1] * u[p + 1][s];
						overline_u[i][j] -= overline_bc_1[i][j][p][s] * xi[p + 1][s] * u[p][s + 1];
						overline_u[i][j] -= overline_cc_1[i][j][p][s] * xi[p + 1][s] * u[p + 1][s];
						//2 �����������
						overline_u[i][j] -= overline_aa_2[i][j][p][s] * xi[p + 1][s + 1] * u[p + 1][s + 1];
						overline_u[i][j] -= overline_ab_2[i][j][p][s] * xi[p + 1][s + 1] * u[p][s + 1];
						overline_u[i][j] -= overline_ab_2[i][j][p][s] * xi[p][s + 1] * u[p + 1][s + 1];
						overline_u[i][j] -= overline_ac_2[i][j][p][s] * xi[p + 1][s + 1] * u[p + 1][s];
						overline_u[i][j] -= overline_ac_2[i][j][p][s] * xi[p + 1][s] * u[p + 1][s + 1];
						overline_u[i][j] -= overline_bb_2[i][j][p][s] * xi[p][s + 1] * u[p][s + 1];
						overline_u[i][j] -= overline_bc_2[i][j][p][s] * xi[p][s + 1] * u[p + 1][s];
						overline_u[i][j] -= overline_bc_2[i][j][p][s] * xi[p + 1][s] * u[p][s + 1];
						overline_u[i][j] -= overline_cc_2[i][j][p][s] * xi[p + 1][s] * u[p + 1][s];
					}
				}
			}
		}
		// ������ overline_u_0^(count) � ���� � ���� �������
		for (size_t i = 0; i <= N; ++i)
		{
			for (size_t j = 0; j <= N; ++j)
			{
				file_overline_u << fixed << setprecision(6) << overline_u[i][j] << " ";
			}
		}

		Lasting("Finding the acoustic pressure in X", time);
	}
	file_overline_u.close();

	Lasting("FThe total time of the program", timeBegin);
}
