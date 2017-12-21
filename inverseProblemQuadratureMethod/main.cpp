#include "stdafx.h"
#include "../directProblemQuadratureMethod/basicFunctions.h"
#include "../directProblemQuadratureMethod/Sources.h"
#include "../directProblemQuadratureMethod/taskData.h"
#include "../directProblemQuadratureMethod/matrix_utils.h"
#include "exact_solution.h"
#include "initialValue.h"
#include "arrayLoading.h"
#include "regularization.h"

using namespace std;

int main()
{
	size_t numberOfIterations;
	cout << "Enter the number of iterations ";
	cin >> numberOfIterations;

	float alpha;
	cout << "Enter alpha ";
	cin >> alpha;

	float multiplier;
	cout << "Enter q ";
	cin >> multiplier;

	const Source source;

	// начало счета времени
	clock_t time, timeBegin;
	timeBegin = clock();
	time = clock();

	// выделение памяти
	// выделяем память под основные матрицы
	vector<vector<vector<vector<complex<float>>>>> a(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<float>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<float>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<float>>(NUMBER_PARTITION_POINT + 1, complex<float>()))));

	vector<vector<vector<complex<float>>>> overline_a(NUMBER_PARTITION_POINT + 1,
		vector<vector<complex<float>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<float>>(NUMBER_PARTITION_POINT + 1, complex<float>())));

	vector<vector<complex<float>>> b(NUMBER_PARTITION_POINT + 1,
		vector<complex<float>>(NUMBER_PARTITION_POINT + 1, complex<float>()));

	// выделяем память для значений источников
	vector<vector<vector<complex<float>>>> Source_R(source.numberSource,
		vector<vector<complex<float>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<float>>(NUMBER_PARTITION_POINT + 1, complex<float>())));
	vector<vector<complex<float>>> Source_X(source.numberSource,
		vector<complex<float>>(NUMBER_PARTITION_POINT + 1, complex<float>()));

	// Выделяем память для поля в приемниках
	vector<vector<complex<float>>> overline_u(source.numberSource,
		vector<complex<float>>(NUMBER_PARTITION_POINT + 1, complex<float>()));

	//выделение памяти под массивы производных  F_1, F_2, ...
	vector<vector<vector<complex<float>>>> F_odd(source.numberSource + 1,
		vector<vector<complex<float>>>(N_SQUARED, vector<complex<float>>(N_SQUARED, complex<float>())));
	vector<vector<vector<complex<float>>>> F_even(source.numberSource + 1,
		vector<vector<complex<float>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<float>>(N_SQUARED, complex<float>())));

	//выделение памяти под массивы A и B
	vector<vector<vector<complex<float>>>> A(source.numberSource + 1,
		vector<vector<complex<float>>>(N_SQUARED, vector<complex<float>>(N_SQUARED, complex<float>())));
	vector<vector<complex<float>>> B(N_SQUARED, vector<complex<float>>(N_SQUARED, complex<float>()));
	vector<vector<complex<float>>> inverseMatrixB(N_SQUARED, vector<complex<float>>(N_SQUARED, complex<float>()));
	
	// память для хранения значений основного оператора
	vector<vector<complex<float>>> F_part_odd(source.numberSource,
		vector<complex<float>>(N_SQUARED, complex<float>()));
	vector<vector<complex<float>>> F_part_even(source.numberSource,
		vector<complex<float>>(NUMBER_PARTITION_POINT + 1, complex<float>()));

	// память для b_0, b_1,...
	vector<vector<complex<float>>> b_right(source.numberSource + 1,
		vector<complex<float>>(N_SQUARED, complex<float>()));

	// память для u^(1), u^(2), u^(3)
	vector<vector<vector<complex<float>>>> u(source.numberSource + 1,
		vector<vector<complex<float>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<float>>(NUMBER_PARTITION_POINT + 1, complex<float>())));

	// память для xi
	vector<vector<complex<float>>> xi(NUMBER_PARTITION_POINT + 1,
		vector<complex<float>>(NUMBER_PARTITION_POINT + 1, complex<float>()));

	// память для перенумерованных переменных
	vector<vector<complex<float>>> numbered_u(source.numberSource,
		vector<complex<float>>(N_SQUARED, complex<float>()));
	vector<complex<float>> numbered_xi(N_SQUARED, complex<float>());

	Lasting("Time allocation", time);

	LoadData(source.numberSource, a, overline_a, b, Source_R, Source_X, overline_u);
	Lasting("Download time", time);

	// Начало вычислительной части
	InitialValueU(source.numberSource, u, Source_R);
	InitialValueXi(xi);

	// начало основных итераций
	for (size_t iteration = 0; iteration < numberOfIterations; ++iteration)
	{
		cout << endl;
		cout << "Iteration number " << (iteration + 1) << endl;
		cout << "alpha= " << alpha << endl;

		//строим левую часть СЛАУ основного метода Ньютона
		GetJacobian(source.numberSource, a, overline_a, b, xi, u, F_odd, F_even);
		Lasting("The counting time of the Jacobian matrices", time);

		GetMatrixA(source.numberSource, F_odd, F_even, A, alpha);
		Lasting("Counting time of matrices A", time);

		GetMatrixB(F_odd, F_even, B, alpha);
		Lasting("Counting time of matrices B", time);

		cout << "Calculate the left side" << endl;

		//строим правую часть СЛАУ основного метода Ньютона
		GetOperatorF(source.numberSource, a, overline_a, b, xi, u, overline_u, Source_R, Source_X, F_part_odd, F_part_even);
		Lasting("Counting time of the main matrix", time);

		Renumbering(xi, numbered_xi);
		for (size_t count = 0; count < source.numberSource; ++count)
		{
			Renumbering(u[count], numbered_u[count]);
		}

		GetValueDerivedFunction(source.numberSource, numbered_xi, numbered_u, F_odd, F_even, F_part_odd, F_part_even);

		Getb(source.numberSource, F_odd, F_even, F_part_odd, F_part_even, b_right);
		Lasting("Calculation time right side", time);

		InvertMatrix(B, inverseMatrixB);
		GetXi(source.numberSource, A, inverseMatrixB, b_right, numbered_xi);
		Lasting("Time of xi", time);

		GetU(source.numberSource, A, inverseMatrixB, b_right, numbered_xi, numbered_u);
		Lasting("Time of u", time);

		alpha = alpha * multiplier;

		InverseRenumbering(numbered_xi, xi);
		for (size_t count = 0; count < source.numberSource; ++count)
		{
			InverseRenumbering(numbered_u[count], u[count]);
		}

		ProjectionXi(xi);

		PrintXi(xi, iteration);
		Lasting("Calculation time solutions", time);
	}
	Lasting("The total time of the program", timeBegin);
}
