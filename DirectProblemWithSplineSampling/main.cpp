#include "stdafx.h"
#include "basicFunctions.h"
#include "Sources.h"
#include "taskData.h"
#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"

using namespace std;

int main()
{
	const size_t N = NUMBER_PARTITION_POINT;
	const double h = DOMAIN_IN_HOMOGENEITY / N;
	const double stepReceiver = 0.1 / N;

	const Source source;

	// ��������� ������ ��� ������������� ���� u
	vector<vector<complex<double>>> u(N + 1, vector<complex<double>>(N + 1, complex<double>()));

	// ������� ������� ������� \xi
	vector<vector<double>> xi(N + 1, vector<double>(N + 1, 0.0));
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			xi[i][j] = exp(-64.0 * (i * h - 0.6) * (i * h - 0.6) - 64.0 * (j * h - 0.6) * (j * h - 0.6));
		}
	}
	//�������� ������ ������� � ����
	ofstream f_xi("exact_xi.txt");
	for (size_t i = 0; i <= N; ++i)
	{
		for (size_t j = 0; j <= N; ++j)
		{
			f_xi << fixed << setprecision(6) << xi[i][j] << " ";
		}
	}
	f_xi.close();

	// �������� ������ ��� �������� �������
	vector<vector<vector<vector<complex<double>>>>> aa_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> ab_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> ac_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> bb_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> bc_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> cc_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> aa_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> ab_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> ac_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> bb_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> bc_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> cc_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_aa_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_ab_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_ac_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_bb_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_bc_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_cc_1(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_aa_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_ab_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_ac_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_bb_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_bc_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	vector<vector<vector<vector<complex<double>>>>> overline_cc_2(N + 1,
		vector<vector<vector<complex<double>>>>(N + 1, vector<vector<complex<double>>>(N + 1,
			vector<complex<double>>(N + 1, complex<double>()))));

	
	// ������ ���������� �������� ������

	// ������ ����� �������
	clock_t timeStart, timeFinish, timeBegin;
	timeBegin = clock();
	timeStart = clock();

	// ���������� �������� aa1, ab1,...,cc1
	// ������������ ������� �������� �������
	double x1, x2;
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
					aa_1[i][j][p][s] = -27.0 * G(i * h, j * h, x1, x2) / 9.0;
					ab_1[i][j][p][s] = -27.0 * G(i * h, j * h, x1, x2) / 9.0;
					ac_1[i][j][p][s] = -27.0 * G(i * h, j * h, x1, x2) / 9.0;
					bb_1[i][j][p][s] = -27.0 * G(i * h, j * h, x1, x2) / 9.0;
					bc_1[i][j][p][s] = -27.0 * G(i * h, j * h, x1, x2) / 9.0;
					cc_1[i][j][p][s] = -27.0 * G(i * h, j * h, x1, x2) / 9.0;

					overline_aa_1[i][j][p][s] = -27.0 * G(1.1 + i * stepReceiver, j * h, x1, x2) / 9.0;
					overline_ab_1[i][j][p][s] = -27.0 * G(1.1 + i * stepReceiver, j * h, x1, x2) / 9.0;
					overline_ac_1[i][j][p][s] = -27.0 * G(1.1 + i * stepReceiver, j * h, x1, x2) / 9.0;
					overline_bb_1[i][j][p][s] = -27.0 * G(1.1 + i * stepReceiver, j * h, x1, x2) / 9.0;
					overline_bc_1[i][j][p][s] = -27.0 * G(1.1 + i * stepReceiver, j * h, x1, x2) / 9.0;
					overline_cc_1[i][j][p][s] = -27.0 * G(1.1 + i * stepReceiver, j * h, x1, x2) / 9.0;

					x1 = h * 0.2 + p * h;
					x2 = h * 0.2 + s * h;
					aa_1[i][j][p][s] += 25.0 * G(i * h, j * h, x1, x2) * 0.04;
					ab_1[i][j][p][s] += 25.0 * G(i * h, j * h, x1, x2) * 0.04;
					ac_1[i][j][p][s] += 25.0 * G(i * h, j * h, x1, x2) * 0.12;
					bb_1[i][j][p][s] += 25.0 * G(i*h, j*h, x1, x2)*(0.2*0.2);
					bc_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.6);
					cc_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.6);

					overline_aa_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_ab_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_ac_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.6);
					overline_bb_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_bc_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.6);
					overline_cc_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.6);

					x1 = h*0.6 + p*h;
					x2 = h*0.2 + s*h;
					aa_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.6);
					ab_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.2);
					ac_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.2);
					bb_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					bc_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					cc_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);

					overline_aa_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.6);
					overline_ab_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.2);
					overline_ac_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.2);
					overline_bb_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_bc_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_cc_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);

					x1 = h*0.2 + p*h;
					x2 = h*0.6 + s*h;
					aa_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					ab_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.6);
					ac_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					bb_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.6);
					bc_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.2);
					cc_1[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);

					overline_aa_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_ab_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.6);
					overline_ac_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_bb_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.6);
					overline_bc_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.2);
					overline_cc_1[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);

					aa_1[i][j][p][s] *= 0.01041667*h*h;
					ab_1[i][j][p][s] *= 0.01041667*h*h;
					ac_1[i][j][p][s] *= 0.01041667*h*h;
					bb_1[i][j][p][s] *= 0.01041667*h*h;
					bc_1[i][j][p][s] *= 0.01041667*h*h;
					cc_1[i][j][p][s] *= 0.01041667*h*h;

					overline_aa_1[i][j][p][s] *= 0.01041667*h*h;
					overline_ab_1[i][j][p][s] *= 0.01041667*h*h;
					overline_ac_1[i][j][p][s] *= 0.01041667*h*h;
					overline_bb_1[i][j][p][s] *= 0.01041667*h*h;
					overline_bc_1[i][j][p][s] *= 0.01041667*h*h;
					overline_cc_1[i][j][p][s] *= 0.01041667*h*h;

					//������ ������������
					x1 = -h*0.33333333 + p*h + h;
					x2 = -h*0.33333333 + s*h + h;
					aa_2[i][j][p][s] = -27.0*G(i*h, j*h, x1, x2)*(0.33333333*0.33333333);
					ab_2[i][j][p][s] = -27.0*G(i*h, j*h, x1, x2)*(0.33333333*0.33333333);
					ac_2[i][j][p][s] = -27.0*G(i*h, j*h, x1, x2)*(0.33333333*0.33333333);
					bb_2[i][j][p][s] = -27.0*G(i*h, j*h, x1, x2)*(0.33333333*0.33333333);
					bc_2[i][j][p][s] = -27.0*G(i*h, j*h, x1, x2)*(0.33333333*0.33333333);
					cc_2[i][j][p][s] = -27.0*G(i*h, j*h, x1, x2)*(0.33333333*0.33333333);

					overline_aa_2[i][j][p][s] = -27.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.33333333*0.33333333);
					overline_ab_2[i][j][p][s] = -27.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.33333333*0.33333333);
					overline_ac_2[i][j][p][s] = -27.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.33333333*0.33333333);
					overline_bb_2[i][j][p][s] = -27.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.33333333*0.33333333);
					overline_bc_2[i][j][p][s] = -27.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.33333333*0.33333333);
					overline_cc_2[i][j][p][s] = -27.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.33333333*0.33333333);

					x1 = -h*0.2 + p*h + h;
					x2 = -h*0.2 + s*h + h;
					aa_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					ab_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					ac_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.6);
					bb_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					bc_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.6);
					cc_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.6);

					overline_aa_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_ab_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_ac_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.6);
					overline_bb_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_bc_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.6);
					overline_cc_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.6);

					x1 = -h*0.2 + p*h + h;
					x2 = -h*0.6 + s*h + h;
					aa_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.6);
					ab_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.2);
					ac_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.2);
					bb_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					bc_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					cc_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);

					overline_aa_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.6);
					overline_ab_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.2);
					overline_ac_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.2);
					overline_bb_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_bc_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_cc_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);

					x1 = -h*0.6 + p*h + h;
					x2 = -h*0.2 + s*h + h;
					aa_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					ab_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.6);
					ac_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);
					bb_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.6);
					bc_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.6*0.2);
					cc_2[i][j][p][s] += 25.0*G(i*h, j*h, x1, x2)*(0.2*0.2);

					overline_aa_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_ab_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.6);
					overline_ac_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);
					overline_bb_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.6);
					overline_bc_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.6*0.2);
					overline_cc_2[i][j][p][s] += 25.0*G(1.1 + i*stepReceiver, j*h, x1, x2)*(0.2*0.2);


					aa_2[i][j][p][s] *= 0.01041667*h*h;
					ab_2[i][j][p][s] *= 0.01041667*h*h;
					ac_2[i][j][p][s] *= 0.01041667*h*h;
					bb_2[i][j][p][s] *= 0.01041667*h*h;
					bc_2[i][j][p][s] *= 0.01041667*h*h;
					cc_2[i][j][p][s] *= 0.01041667*h*h;

					overline_aa_2[i][j][p][s] *= 0.01041667*h*h;
					overline_ab_2[i][j][p][s] *= 0.01041667*h*h;
					overline_ac_2[i][j][p][s] *= 0.01041667*h*h;
					overline_bb_2[i][j][p][s] *= 0.01041667*h*h;
					overline_bc_2[i][j][p][s] *= 0.01041667*h*h;
					overline_cc_2[i][j][p][s] *= 0.01041667*h*h;
				}
			}
		}
	}

	//�������� ����� ������
	timeFinish = clock();
	double d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Time calculation of basic matrices " << d << endl;
	timeStart = clock();
	//
	//�������� ������������ � �����
	//
	ofstream f_aa_1("matrix_aa_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_aa_1 << fixed << setprecision(12) << aa_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_aa_1.close();

	ofstream f_ab_1("matrix_ab_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_ab_1 << fixed << setprecision(12) << ab_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_ab_1.close();

	ofstream f_ac_1("matrix_ac_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_ac_1 << fixed << setprecision(12) << ac_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_ac_1.close();

	ofstream f_bb_1("matrix_bb_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_bb_1 << fixed << setprecision(12) << bb_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_bb_1.close();

	ofstream f_bc_1("matrix_bc_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_bc_1 << fixed << setprecision(12) << bc_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_bc_1.close();

	ofstream f_cc_1("matrix_cc_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_cc_1 << fixed << setprecision(12) << cc_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_cc_1.close();

	ofstream f_aa_2("matrix_aa_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_aa_2 << fixed << setprecision(12) << aa_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_aa_2.close();

	ofstream f_ab_2("matrix_ab_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_ab_2 << fixed << setprecision(12) << ab_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_ab_2.close();

	ofstream f_ac_2("matrix_ac_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_ac_2 << fixed << setprecision(12) << ac_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_ac_2.close();

	ofstream f_bb_2("matrix_bb_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_bb_2 << fixed << setprecision(12) << bb_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_bb_2.close();

	ofstream f_bc_2("matrix_bc_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_bc_2 << fixed << setprecision(12) << bc_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_bc_2.close();

	ofstream f_cc_2("matrix_cc_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_cc_2 << fixed << setprecision(12) << cc_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_cc_2.close();

	ofstream f_overline_aa_1("matrix_overline_aa_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_aa_1 << fixed << setprecision(12) << overline_aa_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_aa_1.close();

	ofstream f_overline_ab_1("matrix_overline_ab_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_ab_1 << fixed << setprecision(12) << overline_ab_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_ab_1.close();

	ofstream f_overline_ac_1("matrix_overline_ac_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_ac_1 << fixed << setprecision(12) << overline_ac_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_ac_1.close();

	ofstream f_overline_bb_1("matrix_overline_bb_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_bb_1 << fixed << setprecision(12) << overline_bb_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_bb_1.close();

	ofstream f_overline_bc_1("matrix_overline_bc_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_bc_1 << fixed << setprecision(12) << overline_bc_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_bc_1.close();

	ofstream f_overline_cc_1("matrix_overline_cc_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_cc_1 << fixed << setprecision(12) << overline_cc_1[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_cc_1.close();

	ofstream f_overline_aa_2("matrix_overline_aa_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_aa_2 << fixed << setprecision(12) << overline_aa_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_aa_2.close();

	ofstream f_overline_ab_2("matrix_overline_ab_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_ab_2 << fixed << setprecision(12) << overline_ab_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_ab_2.close();

	ofstream f_overline_ac_2("matrix_overline_ac_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_ac_2 << fixed << setprecision(12) << overline_ac_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_ac_2.close();

	ofstream f_overline_bb_2("matrix_overline_bb_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_bb_2 << fixed << setprecision(12) << overline_bb_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_bb_2.close();

	ofstream f_overline_bc_2("matrix_overline_bc_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_bc_2 << fixed << setprecision(12) << overline_bc_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_bc_2.close();

	ofstream f_overline_cc_2("matrix_overline_cc_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
				{
					f_overline_cc_2 << fixed << setprecision(12) << overline_cc_2[i][j][p][s] << " ";
				}
			}
		}
	}
	f_overline_cc_2.close();
	//�������� ����� ������
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Time recording the basic matrix file " << d << endl;
	timeStart = clock();
	//
	// ���� ������� ��������� � R � X
	//
	// ������ ��������
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
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_Source_02 << fixed << setprecision(12) << f_01(1.1 + i*stepReceiver, j*h) << " ";
		}
	}
	f_Source_02.close();

	// ������ ��������
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
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_Source_04 << fixed << setprecision(12) << f_02(1.1 + i*stepReceiver, j*h) << " ";
		}
	}
	f_Source_04.close();

	// ������ ��������
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
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_Source_06 << fixed << setprecision(12) << f_03(1.1 + i*stepReceiver, j*h) << " ";
		}
	}
	f_Source_06.close();

	// ��������� ��������
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
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_Source_08 << fixed << setprecision(12) << f_04(1.1 + i*stepReceiver, j*h) << " ";
		}
	}
	f_Source_08.close();

	// ����� ��������
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
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_Source_10 << fixed << setprecision(12) << f_05(1.1 + i*stepReceiver, j*h) << " ";
		}
	}
	f_Source_10.close();

	//�������� ����� ������
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "The computation time of the source function " << d << endl;
	timeStart = clock();
	//
	// ��� ���������� u^(1) ���������� ���� �������� ������� * u^(1) = ������ �����
	//
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEequation[ii]
	//
	int N_squared;
	N_squared = (N + 1)*(N + 1);
	//
	// ������ ��� ���������������� ������ � ������ �����
	complex<double> *numbered_u, *rightPartEquation;
	rightPartEquation = new complex<double>[N_squared];
	numbered_u = new complex<double>[N_squared];
	//������ ��� �������� �������
	complex<double> **substantiveMatrix;
	substantiveMatrix = new complex<double> *[N_squared];
	for (int i = 0; i<N_squared; i++)
	{
		substantiveMatrix[i] = new complex<double>[N_squared];
	}
	// ��������� ������ ��� overline_u_0
	complex<double> **overline_u;
	overline_u = new complex<double> *[N + 1];
	for (int i = 0; i <= N; i++)
	{
		overline_u[i] = new complex<double>[N + 1];
	}
	//
	//���� �������� �������
	//
	int ii, jj;
	//
	// ������� ��� ������ ��������
	//
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			//
			// ������� O(0; 0)
			//
			//1 �����������
			substantiveMatrix[ii][0] = aa_1[i][j][0][0] * xi[0][0];
			substantiveMatrix[ii][0] += ab_1[i][j][0][0] * xi[0][1];
			substantiveMatrix[ii][0] += ac_1[i][j][0][0] * xi[1][0];
			//2 �����������
			//
			// ������� A(N; 0)
			//
			int auxInd = N*(N + 1);
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
			//
			// ������� B(0; N)
			//
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
			//
			// ������� C(N; N)
			//
			auxInd = N*(N + 1) + N;
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
	//
	// ������� ��� ������ ��������
	//
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			//
			// ������� OA(p; 0)
			//
			for (int p = 1; p<N; p++)
			{
				int auxInd = p*(N + 1);
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
			//
			// ������� OB(0; s)
			//
			for (int s = 1; s<N; s++)
			{
				int auxInd = s;
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
			for (int p = 1; p<N; p++)
			{
				int auxInd = p*(N + 1) + N;
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
			//
			// ������� AC(N; s)
			//
			for (int s = 1; s<N; s++)
			{
				int auxInd = N*(N + 1) + s;
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
	//
	// ������� ��� ���������� ����� ��������
	//
	for (int i = 1; i<N; i++)
	{
		for (int j = 1; j<N; j++)
		{
			ii = i*(N + 1) + j;
			for (int p = 1; p<N; p++)
			{
				for (int s = 1; s<N; s++)
				{
					jj = p*(N + 1) + s;
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
	//
	// ��������� ������� � ������� ���������
	//
	for (int ii = 0; ii<N_squared; ii++)
	{
		substantiveMatrix[ii][ii] += 1;
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "The computation time of the matrix inside the cube " << d << endl;
	timeStart = clock();
	//
	// ������� ������������ ���� � ����������
	//
	////////////////////////////////////////////////////////
	///��� ������� ���������
	///////////////////////////////////////////////////////
	// ���������� ������ �����
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			rightPartEquation[ii] = f_01(i*h, j*h);
		}
	}
	// ���������� u^{(1)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// �������� �������������
	//
	for (int ii = 0; ii<N_squared; ii++)
	{
		int coordinate_x = ii / (N + 1);
		int coordinate_y = ii % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[ii];
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 1 source " << d << endl;
	timeStart = clock();
	//
	// ������� overline_u_0
	//
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			overline_u[i][j] = f_01(1.1 + i*stepReceiver, j*h);
		}
	}

	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
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
	// ������ overline_u_0^(1) � ���� � ���� �������
	ofstream f_overline_u_1("matrix_overline_u_1.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_overline_u_1 << fixed << setprecision(12) << overline_u[i][j] << " ";
		}
	}
	f_overline_u_1.close();
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 1 source " << d << endl;
	timeStart = clock();
	////////////////////////////////////////////////////////
	///��� ������� ���������
	///////////////////////////////////////////////////////
	// ���������� ������ �����
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			rightPartEquation[ii] = f_02(i*h, j*h);
		}
	}
	// ���������� u^{(2)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// �������� �������������
	//
	for (int ii = 0; ii<N_squared; ii++)
	{
		int coordinate_x = ii / (N + 1);
		int coordinate_y = ii % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[ii];
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 2 source " << d << endl;
	timeStart = clock();
	//
	// ������� overline_u_0
	//
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			overline_u[i][j] = f_02(1.1 + i*stepReceiver, j*h);
		}
	}

	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
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
	// ������ overline_u_0^(2) � ���� � ���� �������
	ofstream f_overline_u_2("matrix_overline_u_2.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_overline_u_2 << fixed << setprecision(12) << overline_u[i][j] << " ";
		}
	}
	f_overline_u_2.close();
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 2 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///��� �������� ���������
	///////////////////////////////////////////////////////
	// ���������� ������ �����
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			rightPartEquation[ii] = f_03(i*h, j*h);
		}
	}
	// ���������� u^{(3)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// �������� �������������
	//
	for (int ii = 0; ii<N_squared; ii++)
	{
		int coordinate_x = ii / (N + 1);
		int coordinate_y = ii % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[ii];
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 3 source " << d << endl;
	timeStart = clock();
	//
	// ������� overline_u_0
	//
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			overline_u[i][j] = f_03(1.1 + i*stepReceiver, j*h);
		}
	}

	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
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
	// ������ overline_u_0^(3) � ���� � ���� �������
	ofstream f_overline_u_3("matrix_overline_u_3.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_overline_u_3 << fixed << setprecision(12) << overline_u[i][j] << " ";
		}
	}
	f_overline_u_3.close();
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 3 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///��� ��������� ���������
	///////////////////////////////////////////////////////
	// ���������� ������ �����
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			rightPartEquation[ii] = f_04(i*h, j*h);
		}
	}
	// ���������� u^{(4)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// �������� �������������
	//
	for (int ii = 0; ii<N_squared; ii++)
	{
		int coordinate_x = ii / (N + 1);
		int coordinate_y = ii % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[ii];
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 4 source " << d << endl;
	timeStart = clock();
	//
	// ������� overline_u_0
	//
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			overline_u[i][j] = f_04(1.1 + i*stepReceiver, j*h);
		}
	}

	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
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
	// ������ overline_u_0^(4) � ���� � ���� �������
	ofstream f_overline_u_4("matrix_overline_u_4.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_overline_u_4 << fixed << setprecision(12) << overline_u[i][j] << " ";
		}
	}
	f_overline_u_4.close();
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 4 source " << d << endl;
	timeStart = clock();

	////////////////////////////////////////////////////////
	///��� ������ ���������
	///////////////////////////////////////////////////////
	// ���������� ������ �����
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			ii = i*(N + 1) + j;
			rightPartEquation[ii] = f_05(i*h, j*h);
		}
	}
	// ���������� u^{(5)}
	SolveSlauGaussa(substantiveMatrix, N_squared, rightPartEquation, numbered_u);
	//
	// �������� �������������
	//
	for (int ii = 0; ii<N_squared; ii++)
	{
		int coordinate_x = ii / (N + 1);
		int coordinate_y = ii % (N + 1);
		u[coordinate_x][coordinate_y] = numbered_u[ii];
	}
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in R for 5 source " << d << endl;
	timeStart = clock();
	//
	// ������� overline_u_0
	//
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			overline_u[i][j] = f_05(1.1 + i*stepReceiver, j*h);
		}
	}

	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int p = 0; p<N; p++)
			{
				for (int s = 0; s<N; s++)
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
	// ������ overline_u_0^(5) � ���� � ���� �������
	ofstream f_overline_u_5("matrix_overline_u_5.txt");
	for (int i = 0; i <= N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			f_overline_u_5 << fixed << setprecision(12) << overline_u[i][j] << " ";
		}
	}
	f_overline_u_5.close();
	timeFinish = clock();
	d = (double)(timeFinish - timeStart) / CLOCKS_PER_SEC;
	cout << "Finding the acoustic pressure in X for 5 source " << d << endl;
	timeStart = clock();

}
