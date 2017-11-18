#include "stdafx.h"
#include "basicFunctions.h"
#include "Sources.h"
#include "taskData.h"
#include "basicArrays.h"
#include "exact_solution.h"
#include "array_utils.h"
#include "directProblem_utils.h"

using namespace std;

int main()
{
	const size_t N = NUMBER_PARTITION_POINT;
	const Source source;

	// ��������� ������ ��� ������������� ���� u
	array<array<complex<double>, SPLITTING + 1>, SPLITTING + 1> u;

	// ������� ������� ������� xi
	array<array<complex<double>, SPLITTING + 1>, SPLITTING + 1> xi;
	GetExactSolution(xi);
	WriteSolutionFile(xi);
	
	// ������ ����� �������
	clock_t time = clock();
	clock_t timeBegin = clock();

	BasicArrays basicArrays;
	GetBasicArrays(basicArrays);
	Lasting("Time calculation of basic matrices", time);

	WriteBasicArraysFile(basicArrays);
	Lasting("Download time major arrays", time);

	WriteSourceValues(source);
	Lasting("The computation time of the source function", time);

	//�������� ����� ������
	Lasting("The computation time of the source function", time);

	// ��� ���������� u^(1) ���������� ���� �������� ������� * u^(1) = ������ �����
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEequation[ii]
	array<complex<double>, N_SQUARED> rightPartEquation{};
	array<complex<double>, N_SQUARED> numbered_u{};
	array<array<complex<double>, N_SQUARED>, N_SQUARED> substantiveMatrix{};
	array<array<complex<double>, N_SQUARED>, N_SQUARED> overline_u{};

	//���� �������� �������
	GetSubstantiveMatrix(basicArrays, xi, substantiveMatrix);
	Lasting("The computation time of the matrix inside the squared", time);

	// ������� ������������ ���� � ����������
	ofstream file_overline_u("matrix_overline_u.txt");
	file_overline_u << fixed << setprecision(6);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		// ���������� ������ �����
		GetRightPartEquation(source, count, rightPartEquation);
		
		// ���������� u^{(count)}
		SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);

		// �������� �������������
		InverseRenumbering(numbered_u, u);
		Lasting("Finding the acoustic pressure in R", time);

		GetOverlineU(source, count, basicArrays, xi, u, overline_u);
		for (size_t i = 0; i <= SPLITTING; ++i)
		{
			for (size_t j = 0; j <= SPLITTING; ++j)
			{
				file_overline_u << overline_u[i][j] << " ";
			}
		}

		Lasting("Finding the acoustic pressure in X", time);
	}
	file_overline_u.close();

	Lasting("FThe total time of the program", timeBegin);
}
