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

	// выделение пам€ти дл€ акустического пол€ u
	array<array<complex<double>, SPLITTING + 1>, SPLITTING + 1> u;

	// задание точного решени€ xi
	array<array<complex<double>, SPLITTING + 1>, SPLITTING + 1> xi;
	GetExactSolution(xi);
	WriteSolutionFile(xi);
	
	// начало счета времени
	clock_t time = clock();
	clock_t timeBegin = clock();

	BasicArrays basicArrays;
	GetBasicArrays(basicArrays);
	Lasting("Time calculation of basic matrices", time);

	WriteBasicArraysFile(basicArrays);
	Lasting("Download time major arrays", time);

	WriteSourceValues(source);
	Lasting("The computation time of the source function", time);

	//печатаем врем€ работы
	Lasting("The computation time of the source function", time);

	// дл€ нахождени€ u^(1) составл€ем —Ћј” основна€ матрица * u^(1) = правой части
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEequation[ii]
	array<complex<double>, N_SQUARED> rightPartEquation{};
	array<complex<double>, N_SQUARED> numbered_u{};
	array<array<complex<double>, N_SQUARED>, N_SQUARED> substantiveMatrix{};
	array<array<complex<double>, N_SQUARED>, N_SQUARED> overline_u{};

	//счет основной матрицы
	GetSubstantiveMatrix(basicArrays, xi, substantiveMatrix);
	Lasting("The computation time of the matrix inside the squared", time);

	// Ќаходим акустическое поле в приемниках
	ofstream file_overline_u("matrix_overline_u.txt");
	file_overline_u << fixed << setprecision(6);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		// нахождение правой части
		GetRightPartEquation(source, count, rightPartEquation);
		
		// нахождение u^{(count)}
		SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);

		// ќбратна€ перенумераци€
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
