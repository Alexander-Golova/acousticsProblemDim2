#include "stdafx.h"

#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"
#include <math.h>

using namespace std;


BOOST_AUTO_TEST_SUITE(Test)

// проверка функции, находящей решение СЛАУ методом Гаусса
BOOST_AUTO_TEST_CASE(Solve_Slau_Gaussa_)
{
	// для единичной матрицы
	{
		const size_t dim = 3;
		vector<vector<complex<float>>> matrix1(dim, vector<complex<float>>(dim, complex<float>()));
		matrix1 = { { (1.0f, 0.0f), (0.0f, 0.0f), (0.0f, 0.0f) },
		{ (0.0f, 0.0f), (1.0f, 0.0f), (0.0f, 0.0f) },
		{ (0.0f, 0.0f), (0.0f, 0.0f), (1.0f, 0.0f) } };

		vector<complex<float>> rhs(dim, complex<float>());
		rhs = { (5.0f, 0.0f), (0.0f, 3.0f), (4.0f, 0.0f) };

		vector<complex<float>> correctSolution(dim, complex<float>());
		correctSolution = { (5.0f, 0.0f), (0.0f, 3.0f), (4.0f, 0.0f) };

		vector<complex<float>> solution(dim, complex<float>());
		SolveSlauGaussa(matrix1, rhs, solution);

		BOOST_CHECK(solution == correctSolution);
	}

}

BOOST_AUTO_TEST_SUITE_END()
