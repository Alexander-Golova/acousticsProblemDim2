#include "stdafx.h"
#include "basicFunctions.h"
#include "basicArrays.h"
#include "array_utils.h"
#include "directProblem_utils.h"
#include "Sources.h"
#include "taskData.h"
#include"exact_solution.h"
#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"

using namespace std;

int main()
{
	const Source source;

	vector<vector<complex<double>>> u(NUMBER_PARTITION_POINT + 1,
		vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()));

	vector<vector<double>> xi(NUMBER_PARTITION_POINT + 1,
		vector<double>(NUMBER_PARTITION_POINT + 1, 0.0));

	GetExactSolution(xi);

	WriteSolutionFile(xi);

	clock_t timeStart, timeFinish, timeBegin;
	timeBegin = clock();
	timeStart = clock();

	vector<vector<vector<vector<complex<double>>>>> a(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<complex<double>>>> overline_a(NUMBER_PARTITION_POINT + 1,
		vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
		vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>())));

	vector<vector<complex<double>>> b(NUMBER_PARTITION_POINT + 1,
		vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()));

	GetBasicArrays(a, overline_a, b);

	timeFinish = clock();
	Lasting("Time calculation of basic matrices", timeStart, timeFinish);
	timeStart = clock();

	WriteBasicArraysFile(a, overline_a, b);

	timeFinish = clock();
	Lasting("Download time major arrays", timeStart, timeFinish);
	timeStart = clock();

	WriteSourceValues(source);

	timeFinish = clock();
	Lasting("The computation time of the source function", timeStart, timeFinish);
	timeStart = clock();

	// для нахождения u^(1) составляем СЛАУ основная матрица * u^(1) = правой части
	// substantiveMatrix[ii][jj] * numbered_u[jj] = rightPartEquation[ii]

	const size_t N_squared = (NUMBER_PARTITION_POINT + 1) * (NUMBER_PARTITION_POINT + 1);
	vector<complex<double>> rightPartEquation(N_squared, complex<double>());
	vector<complex<double>> numbered_u(N_squared);
	vector<vector<complex<double>>> substantiveMatrix(N_squared, vector<complex<double>>(N_squared, complex<double>()));
	vector<complex<double>> overline_u(NUMBER_PARTITION_POINT + 1, complex<double>());

	GetSubstantiveMatrix(a, b, xi, substantiveMatrix);

	timeFinish = clock();
	Lasting("The computation time of the matrix inside the squared", timeStart, timeFinish);
	timeStart = clock();

	ofstream file_overline_u("matrix_overline_u.txt");
	file_overline_u << fixed << setprecision(6);
	for (size_t count = 0; count < source.numberSource; ++count)
	{
		GetRightPartEquation(source, count, rightPartEquation);

		SolveSlauGaussa(substantiveMatrix, rightPartEquation, numbered_u);

		InverseRenumbering(numbered_u, u);

		timeFinish = clock();
		Lasting("Finding the acoustic pressure in R", timeStart, timeFinish);
		timeStart = clock();

		GetOverlineU(source, count, overline_a, xi, u, overline_u);

		for (size_t j = 0; j <= NUMBER_PARTITION_POINT; ++j)
		{
			file_overline_u << overline_u[j] << " ";
		}

		timeFinish = clock();
		Lasting("Finding the acoustic pressure in X", timeStart, timeFinish);
		timeStart = clock();
	}
	file_overline_u.close();

	timeFinish = clock();
	Lasting("The total time of the program", timeBegin, timeFinish);
}
