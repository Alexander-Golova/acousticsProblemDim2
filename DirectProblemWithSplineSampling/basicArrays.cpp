#include "stdafx.h"
#include "taskData.h"
#include "basicArrays.h"
#include "basicFunctions.h"

using namespace std;

void GetBasicArrays(BasicArrays & basicArrays)
{
	// Используются квадратурные формулы третьего порядка
	double x1, x2;
	complex<double> temp;
	for (size_t i = 0; i <= SPLITTING; ++i)
	{
		for (size_t j = 0; j <= SPLITTING; ++j)
		{
			for (size_t p = 0; p < SPLITTING; ++p)
			{
				for (size_t s = 0; s < SPLITTING; ++s)
				{
					//первые треугольники
					x1 = h / 3.0 + p * h;
					x2 = h / 3.0 + s * h;
					temp = -27.0 * G(i * h, j * h, x1, x2) / 9.0;
					basicArrays.aa_1[i][j][p][s] = temp;
					basicArrays.ab_1[i][j][p][s] = temp;
					basicArrays.ac_1[i][j][p][s] = temp;
					basicArrays.bb_1[i][j][p][s] = temp;
					basicArrays.bc_1[i][j][p][s] = temp;
					basicArrays.cc_1[i][j][p][s] = temp;

					temp = -27.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2) / 9.0;
					basicArrays.overline_aa_1[i][j][p][s] = temp;
					basicArrays.overline_ab_1[i][j][p][s] = temp;
					basicArrays.overline_ac_1[i][j][p][s] = temp;
					basicArrays.overline_bb_1[i][j][p][s] = temp;
					basicArrays.overline_bc_1[i][j][p][s] = temp;
					basicArrays.overline_cc_1[i][j][p][s] = temp;

					x1 = h * 0.2 + p * h;
					x2 = h * 0.2 + s * h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					basicArrays.aa_1[i][j][p][s] += temp * 0.04;
					basicArrays.ab_1[i][j][p][s] += temp * 0.04;
					basicArrays.ac_1[i][j][p][s] += temp * 0.12;
					basicArrays.bb_1[i][j][p][s] += temp * 0.04;
					basicArrays.bc_1[i][j][p][s] += temp * 0.12;
					basicArrays.cc_1[i][j][p][s] += temp * 0.36;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					basicArrays.overline_aa_1[i][j][p][s] += temp * 0.04;
					basicArrays.overline_ab_1[i][j][p][s] += temp * 0.04;
					basicArrays.overline_ac_1[i][j][p][s] += temp * 0.12;
					basicArrays.overline_bb_1[i][j][p][s] += temp * 0.04;
					basicArrays.overline_bc_1[i][j][p][s] += temp * 0.12;
					basicArrays.overline_cc_1[i][j][p][s] += temp * 0.36;

					x1 = h * 0.6 + p * h;
					x2 = h * 0.2 + s * h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					basicArrays.aa_1[i][j][p][s] += temp * 0.36;
					basicArrays.ab_1[i][j][p][s] += temp * 0.12;
					basicArrays.ac_1[i][j][p][s] += temp * 0.12;
					basicArrays.bb_1[i][j][p][s] += temp * 0.04;
					basicArrays.bc_1[i][j][p][s] += temp * 0.04;
					basicArrays.cc_1[i][j][p][s] += temp * 0.04;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					basicArrays.overline_aa_1[i][j][p][s] += temp * 0.36;
					basicArrays.overline_ab_1[i][j][p][s] += temp * 0.12;
					basicArrays.overline_ac_1[i][j][p][s] += temp * 0.12;
					basicArrays.overline_bb_1[i][j][p][s] += temp * 0.04;
					basicArrays.overline_bc_1[i][j][p][s] += temp * 0.04;
					basicArrays.overline_cc_1[i][j][p][s] += temp * 0.04;

					x1 = h * 0.2 + p * h;
					x2 = h * 0.6 + s * h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					basicArrays.aa_1[i][j][p][s] += temp * 0.04;
					basicArrays.ab_1[i][j][p][s] += temp * 0.12;
					basicArrays.ac_1[i][j][p][s] += temp * 0.04;
					basicArrays.bb_1[i][j][p][s] += temp * 0.36;
					basicArrays.bc_1[i][j][p][s] += temp * 0.12;
					basicArrays.cc_1[i][j][p][s] += temp * 0.04;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					basicArrays.overline_aa_1[i][j][p][s] += temp * 0.04;
					basicArrays.overline_ab_1[i][j][p][s] += temp * 0.12;
					basicArrays.overline_ac_1[i][j][p][s] += temp * 0.04;
					basicArrays.overline_bb_1[i][j][p][s] += temp * 0.36;
					basicArrays.overline_bc_1[i][j][p][s] += temp * 0.12;
					basicArrays.overline_cc_1[i][j][p][s] += temp * 0.04;

					temp = h * h / 96.0;
					basicArrays.aa_1[i][j][p][s] *= temp;
					basicArrays.ab_1[i][j][p][s] *= temp;
					basicArrays.ac_1[i][j][p][s] *= temp;
					basicArrays.bb_1[i][j][p][s] *= temp;
					basicArrays.bc_1[i][j][p][s] *= temp;
					basicArrays.cc_1[i][j][p][s] *= temp;

					basicArrays.overline_aa_1[i][j][p][s] *= temp;
					basicArrays.overline_ab_1[i][j][p][s] *= temp;
					basicArrays.overline_ac_1[i][j][p][s] *= temp;
					basicArrays.overline_bb_1[i][j][p][s] *= temp;
					basicArrays.overline_bc_1[i][j][p][s] *= temp;
					basicArrays.overline_cc_1[i][j][p][s] *= temp;

					//вторые треугольники
					x1 = -h / 3.0 + p * h + h;
					x2 = -h / 3.0 + s * h + h;
					temp = -27.0 * G(i * h, j * h, x1, x2) / 9.0;
					basicArrays.aa_2[i][j][p][s] = temp;
					basicArrays.ab_2[i][j][p][s] = temp;
					basicArrays.ac_2[i][j][p][s] = temp;
					basicArrays.bb_2[i][j][p][s] = temp;
					basicArrays.bc_2[i][j][p][s] = temp;
					basicArrays.cc_2[i][j][p][s] = temp;

					temp = -27.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2) / 9.0;
					basicArrays.overline_aa_2[i][j][p][s] = temp;
					basicArrays.overline_ab_2[i][j][p][s] = temp;
					basicArrays.overline_ac_2[i][j][p][s] = temp;
					basicArrays.overline_bb_2[i][j][p][s] = temp;
					basicArrays.overline_bc_2[i][j][p][s] = temp;
					basicArrays.overline_cc_2[i][j][p][s] = temp;

					x1 = -h * 0.2 + p * h + h;
					x2 = -h * 0.2 + s * h + h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					basicArrays.aa_2[i][j][p][s] += temp * 0.04;
					basicArrays.ab_2[i][j][p][s] += temp * 0.04;
					basicArrays.ac_2[i][j][p][s] += temp * 0.12;
					basicArrays.bb_2[i][j][p][s] += temp * 0.04;
					basicArrays.bc_2[i][j][p][s] += temp * 0.12;
					basicArrays.cc_2[i][j][p][s] += temp * 0.36;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					basicArrays.overline_aa_2[i][j][p][s] += temp * 0.04;
					basicArrays.overline_ab_2[i][j][p][s] += temp * 0.04;
					basicArrays.overline_ac_2[i][j][p][s] += temp * 0.12;
					basicArrays.overline_bb_2[i][j][p][s] += temp * 0.04;
					basicArrays.overline_bc_2[i][j][p][s] += temp * 0.12;
					basicArrays.overline_cc_2[i][j][p][s] += temp * 0.36;

					x1 = -h * 0.2 + p * h + h;
					x2 = -h * 0.6 + s * h + h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					basicArrays.aa_2[i][j][p][s] += temp * 0.36;
					basicArrays.ab_2[i][j][p][s] += temp * 0.12;
					basicArrays.ac_2[i][j][p][s] += temp * 0.12;
					basicArrays.bb_2[i][j][p][s] += temp * 0.04;
					basicArrays.bc_2[i][j][p][s] += temp * 0.04;
					basicArrays.cc_2[i][j][p][s] += temp * 0.04;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					basicArrays.overline_aa_2[i][j][p][s] += temp * 0.36;
					basicArrays.overline_ab_2[i][j][p][s] += temp * 0.12;
					basicArrays.overline_ac_2[i][j][p][s] += temp * 0.12;
					basicArrays.overline_bb_2[i][j][p][s] += temp * 0.04;
					basicArrays.overline_bc_2[i][j][p][s] += temp * 0.04;
					basicArrays.overline_cc_2[i][j][p][s] += temp * 0.04;

					x1 = -h * 0.6 + p * h + h;
					x2 = -h * 0.2 + s * h + h;
					temp = 25.0 * G(i * h, j * h, x1, x2);
					basicArrays.aa_2[i][j][p][s] += temp * 0.04;
					basicArrays.ab_2[i][j][p][s] += temp * 0.12;
					basicArrays.ac_2[i][j][p][s] += temp * 0.04;
					basicArrays.bb_2[i][j][p][s] += temp * 0.36;
					basicArrays.bc_2[i][j][p][s] += temp * 0.12;
					basicArrays.cc_2[i][j][p][s] += temp * 0.04;

					temp = 25.0 * G(RECEIVER + i * stepReceiver, j * h, x1, x2);
					basicArrays.overline_aa_2[i][j][p][s] += temp * 0.04;
					basicArrays.overline_ab_2[i][j][p][s] += temp * 0.12;
					basicArrays.overline_ac_2[i][j][p][s] += temp * 0.04;
					basicArrays.overline_bb_2[i][j][p][s] += temp * 0.36;
					basicArrays.overline_bc_2[i][j][p][s] += temp * 0.12;
					basicArrays.overline_cc_2[i][j][p][s] += temp * 0.04;

					temp = temp = h * h / 96.0;
					basicArrays.aa_2[i][j][p][s] *= temp;
					basicArrays.ab_2[i][j][p][s] *= temp;
					basicArrays.ac_2[i][j][p][s] *= temp;
					basicArrays.bb_2[i][j][p][s] *= temp;
					basicArrays.bc_2[i][j][p][s] *= temp;
					basicArrays.cc_2[i][j][p][s] *= temp;

					basicArrays.overline_aa_2[i][j][p][s] *= temp;
					basicArrays.overline_ab_2[i][j][p][s] *= temp;
					basicArrays.overline_ac_2[i][j][p][s] *= temp;
					basicArrays.overline_bb_2[i][j][p][s] *= temp;
					basicArrays.overline_bc_2[i][j][p][s] *= temp;
					basicArrays.overline_cc_2[i][j][p][s] *= temp;
				}
			}
		}
	}
}

void WriteBasicArraysFile(BasicArrays & basicArrays)
{
	//печатаем коэффициенты в файлы
	ofstream f_matrix("matrix.txt");
	f_matrix << fixed << setprecision(6);
	for (size_t i = 0; i <= SPLITTING; ++i)
	{
		for (size_t j = 0; j <= SPLITTING; ++j)
		{
			for (size_t p = 0; p < SPLITTING; ++p)
			{
				for (size_t s = 0; s < SPLITTING; ++s)
				{
					f_matrix << basicArrays.aa_1[i][j][p][s] << " ";
					f_matrix << basicArrays.ab_1[i][j][p][s] << " ";
					f_matrix << basicArrays.ac_1[i][j][p][s] << " ";
					f_matrix << basicArrays.bb_1[i][j][p][s] << " ";
					f_matrix << basicArrays.bc_1[i][j][p][s] << " ";
					f_matrix << basicArrays.cc_1[i][j][p][s] << " ";
					f_matrix << basicArrays.aa_2[i][j][p][s] << " ";
					f_matrix << basicArrays.ab_2[i][j][p][s] << " ";
					f_matrix << basicArrays.ac_2[i][j][p][s] << " ";
					f_matrix << basicArrays.bb_2[i][j][p][s] << " ";
					f_matrix << basicArrays.bc_2[i][j][p][s] << " ";
					f_matrix << basicArrays.cc_2[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_aa_1[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_ab_1[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_ac_1[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_bb_1[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_bc_1[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_cc_1[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_aa_2[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_ab_2[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_ac_2[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_bb_2[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_bc_2[i][j][p][s] << " ";
					f_matrix << basicArrays.overline_cc_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_matrix.close();
}
