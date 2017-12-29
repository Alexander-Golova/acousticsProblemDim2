#include "stdafx.h"
#include "basicFunctions.h"
#include "basicArrays.h"
#include "array_utils.h"
#include "directProblem_utils.h"
#include "Sources.h"
#include "taskData.h"
#include "exact_solution.h"
#include "matrix_utils.h"

using namespace std;

int main()
{
	// задаём множество источников
	const Source source;

	vector<vector<complex<float>>> u(NUMBER_PARTITION_POINT + 1,
		vector<complex<float>>(NUMBER_PARTITION_POINT + 1, complex<float>()));

	vector<vector<float>> xi(NUMBER_PARTITION_POINT + 1, vector<float>(NUMBER_PARTITION_POINT + 1, 0.0f));

	GetExactSolution(xi);

	WriteSolutionFile(xi);

	clock_t time = clock();
	clock_t timeBegin = clock();

	vector<vector<vector<vector<float>>>> a(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<float>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<float>>(NUMBER_PARTITION_POINT + 1, vector<float>(NUMBER_PARTITION_POINT + 1, 0.0f))));

	vector<vector<vector<vector<float>>>> b(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<float>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<float>>(NUMBER_PARTITION_POINT + 1, vector<float>(NUMBER_PARTITION_POINT + 1, 0.0f))));

	vector<vector<float>> c(NUMBER_PARTITION_POINT + 1, vector<float>(NUMBER_PARTITION_POINT + 1, 0.0f));

	vector<vector<vector<float>>> overline_a(NUMBER_PARTITION_POINT + 1,
		vector<vector<float>>(NUMBER_PARTITION_POINT + 1, vector<float>(NUMBER_PARTITION_POINT + 1, 0.0f)));

	vector<vector<vector<float>>> overline_b(NUMBER_PARTITION_POINT + 1,
		vector<vector<float>>(NUMBER_PARTITION_POINT + 1, vector<float>(NUMBER_PARTITION_POINT + 1, 0.0f)));

	GetBasicArrays(a, b, c, overline_a, overline_b);
	Lasting("Time calculation of basic matrices", time);

	WriteBasicArraysFile(a, b, c, overline_a, overline_b);
	Lasting("Download time major arrays", time);
	return 0;

	WriteSourceValues(source);
	Lasting("The computation time of the source function", time);

	// для нахождения u^(1) составляем СЛАУ основная матрица * u^(1) = правой части
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEquation[ii]

	vector<complex<float>> rightPartEquation(N_SQUARED, complex<float>());
	vector<complex<float>> numbered_u(N_SQUARED);
	vector<vector<complex<float>>> substantiveMatrix(N_SQUARED, vector<complex<float>>(N_SQUARED, complex<float>()));
	vector<complex<float>> overline_u(NUMBER_PARTITION_POINT + 1, complex<float>());

	GetSubstantiveMatrix(a, b, c, xi, substantiveMatrix);
	Lasting("The computation time of the matrix inside the squared", time);

	ofstream file_overline_u("matrix_overline_u.txt");
	file_overline_u << fixed << setprecision(6);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		GetRightPartEquation(source, count, rightPartEquation);

		SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);

		InverseRenumbering(numbered_u, u);
		Lasting("Finding the acoustic pressure in R", time);

		GetOverlineU(source, count, overline_a, overline_b, xi, u, overline_u);
		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			file_overline_u << overline_u[j] << " ";
		}

		Lasting("Finding the acoustic pressure in X", time);
	}
	file_overline_u.close();

	Lasting("The total time of the program", timeBegin);
}
