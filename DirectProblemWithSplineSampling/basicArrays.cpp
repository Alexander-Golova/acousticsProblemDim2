#include "stdafx.h"
#include "taskData.h"
#include "basicArrays.h"
#include "basicFunctions.h"

using namespace std;

void GetBasicArrays(vector<vector<vector<vector<complex<double>>>>> & aa_1,
	vector<vector<vector<vector<complex<double>>>>> & ab_1,
	vector<vector<vector<vector<complex<double>>>>> & ac_1,
	vector<vector<vector<vector<complex<double>>>>> & bb_1,
	vector<vector<vector<vector<complex<double>>>>> & bc_1,
	vector<vector<vector<vector<complex<double>>>>> & cc_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_aa_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_ab_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_ac_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_bb_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_bc_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_cc_1,
	vector<vector<vector<vector<complex<double>>>>> & aa_2,
	vector<vector<vector<vector<complex<double>>>>> & ab_2,
	vector<vector<vector<vector<complex<double>>>>> & ac_2,
	vector<vector<vector<vector<complex<double>>>>> & bb_2,
	vector<vector<vector<vector<complex<double>>>>> & bc_2,
	vector<vector<vector<vector<complex<double>>>>> & cc_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_aa_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_ab_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_ac_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_bb_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_bc_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_cc_2)
{
	// Используются квадратурные формулы третьего порядка
	double x1, x2;
	complex<double> temp;
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t s = 0; s < NUMBER_PARTITION_POINT; ++s)
				{
					//первые треугольники
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

					//вторые треугольники
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
}

void WriteBasicArraysFile(vector<vector<vector<vector<complex<double>>>>> & aa_1,
	vector<vector<vector<vector<complex<double>>>>> & ab_1,
	vector<vector<vector<vector<complex<double>>>>> & ac_1,
	vector<vector<vector<vector<complex<double>>>>> & bb_1,
	vector<vector<vector<vector<complex<double>>>>> & bc_1,
	vector<vector<vector<vector<complex<double>>>>> & cc_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_aa_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_ab_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_ac_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_bb_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_bc_1,
	vector<vector<vector<vector<complex<double>>>>> & overline_cc_1,
	vector<vector<vector<vector<complex<double>>>>> & aa_2,
	vector<vector<vector<vector<complex<double>>>>> & ab_2,
	vector<vector<vector<vector<complex<double>>>>> & ac_2,
	vector<vector<vector<vector<complex<double>>>>> & bb_2,
	vector<vector<vector<vector<complex<double>>>>> & bc_2,
	vector<vector<vector<vector<complex<double>>>>> & cc_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_aa_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_ab_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_ac_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_bb_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_bc_2,
	vector<vector<vector<vector<complex<double>>>>> & overline_cc_2)
{
	//печатаем коэффициенты в файлы
	ofstream f_matrix("matrix.txt");
	f_matrix << fixed << setprecision(6);
	for (size_t i = 0; i <= NUMBER_PARTITION_POINT; ++i)
	{
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			for (size_t p = 0; p < NUMBER_PARTITION_POINT; ++p)
			{
				for (size_t s = 0; s < NUMBER_PARTITION_POINT; ++s)
				{
					f_matrix << aa_1[i][j][p][s] << " ";
					f_matrix << ab_1[i][j][p][s] << " ";
					f_matrix << ac_1[i][j][p][s] << " ";
					f_matrix << bb_1[i][j][p][s] << " ";
					f_matrix << bc_1[i][j][p][s] << " ";
					f_matrix << cc_1[i][j][p][s] << " ";
					f_matrix << aa_2[i][j][p][s] << " ";
					f_matrix << ab_2[i][j][p][s] << " ";
					f_matrix << ac_2[i][j][p][s] << " ";
					f_matrix << bb_2[i][j][p][s] << " ";
					f_matrix << bc_2[i][j][p][s] << " ";
					f_matrix << cc_2[i][j][p][s] << " ";
					f_matrix << overline_aa_1[i][j][p][s] << " ";
					f_matrix << overline_ab_1[i][j][p][s] << " ";
					f_matrix << overline_ac_1[i][j][p][s] << " ";
					f_matrix << overline_bb_1[i][j][p][s] << " ";
					f_matrix << overline_bc_1[i][j][p][s] << " ";
					f_matrix << overline_cc_1[i][j][p][s] << " ";
					f_matrix << overline_aa_2[i][j][p][s] << " ";
					f_matrix << overline_ab_2[i][j][p][s] << " ";
					f_matrix << overline_ac_2[i][j][p][s] << " ";
					f_matrix << overline_bb_2[i][j][p][s] << " ";
					f_matrix << overline_bc_2[i][j][p][s] << " ";
					f_matrix << overline_cc_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_matrix.close();
}
