#include "stdafx.h"
#include "../directProblemQuadratureSamplingMethodNew/basicFunctions.h"
#include "../directProblemQuadratureSamplingMethodNew/Sources.h"
#include "../directProblemQuadratureSamplingMethodNew/taskData.h"
#include "../directProblemQuadratureSamplingMethod/matrix_utils.h"
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

	double alpha;
	cout << "Enter alpha ";
	cin >> alpha;

	double multiplier;
	cout << "Enter q ";
	cin >> multiplier;

	const Source source;

	// начало счета времени
	clock_t time, timeBegin;
	timeBegin = clock();
	time = clock();

	// выделение памяти
	// выделяем память под основные матрицы
	vector<vector<vector<vector<complex<double>>>>> a(NUMBER_PARTITION_POINT + 1,
		vector<vector<vector<complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
				vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<complex<double>>>> overline_a(NUMBER_PARTITION_POINT + 1,
		vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>())));

	vector<vector<complex<double>>> b(NUMBER_PARTITION_POINT + 1,
		vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()));

	// выделяем память для значений источников
	vector<vector<vector<complex<double>>>> Source_R(source.numberSource,
		vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>())));
	vector<vector<complex<double>>> Source_X(source.numberSource,
		vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()));

	// Выделяем память для поля в приемниках
	vector<vector<complex<double>>> overline_u(source.numberSource,
		vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()));

	//выделение памяти под массивы производных  F_1, F_2, ...
	vector<vector<vector<complex<double>>>> F_odd(source.numberSource + 1,
		vector<vector<complex<double>>>(N_SQUARED, vector<complex<double>>(N_SQUARED, complex<double>())));
	vector<vector<vector<complex<double>>>> F_even(source.numberSource + 1,
		vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(N_SQUARED, complex<double>())));

	//выделение памяти под массивы A и B
	vector<vector<vector<complex<double>>>> A(source.numberSource + 1,
		vector<vector<complex<double>>>(N_SQUARED, vector<complex<double>>(N_SQUARED, complex<double>())));
	vector<vector<complex<double>>> B(N_SQUARED, vector<complex<double>>(N_SQUARED, complex<double>()));
	vector<vector<complex<double>>> inverseMatrixB(N_SQUARED, vector<complex<double>>(N_SQUARED, complex<double>()));
	
	// память для хранения значений основного оператора
	vector<vector<complex<double>>> F_part_odd(source.numberSource,
		vector<complex<double>>(N_SQUARED, complex<double>()));
	vector<vector<complex<double>>> F_part_even(source.numberSource,
		vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()));

	// память для b_0, b_1,...
	vector<vector<complex<double>>> b_right(source.numberSource + 1,
		vector<complex<double>>(N_SQUARED, complex<double>()));

	// память для u^(1), u^(2), u^(3)
	vector<vector<vector<complex<double>>>> u(source.numberSource + 1,
		vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>())));

	// память для xi
	vector<vector<complex<double>>> xi(NUMBER_PARTITION_POINT + 1,
		vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()));

	// память для перенумерованных переменных
	vector<vector<complex<double>>> numbered_u(source.numberSource,
		vector<complex<double>>(N_SQUARED, complex<double>()));
	vector<complex<double>> numbered_xi(N_SQUARED, complex<double>());

	Lasting("Time allocation", time);

	// Загрузка данных
	ArrayLoadingA(a);
	ArrayLoadingOverlineA(overline_a);
	ArrayLoadingB(b);
	ArrayLoadingSource(source.numberSource, Source_R, Source_X);
	ArrayLoadingOverlineU(source.numberSource, overline_u);

	//печатаем время работы
	Lasting("Download time", time);

	// Начало вычислительной части
	// начальные значения
	InitialValueU(source.numberSource, u, Source_R);
	InitialValueXi(xi);

	// начало основных итераций
	for (size_t iteration = 0; iteration < numberOfIterations; ++iteration)
	{
		cout << endl;
		cout << "Iteration number " << (iteration + 1) << endl;
		cout << "alpha= " << alpha << endl;

		//строим левую часть СЛАУ основного метода Ньютона
		// получаем матрицы Якобиана
		GetJacobian(source.numberSource, a, overline_a, b, xi, u, F_odd, F_even);
		Lasting("The counting time of the Jacobian matrices", time);

		// находим матрицы А
		GetMatrixA(source.numberSource, F_odd, F_even, A, alpha);
		Lasting("Counting time of matrices A", time);

		// находим матрицу В
		GetMatrixB(F_odd, F_even, B, alpha);
		Lasting("Counting time of matrices B", time);
		cout << "Calculate the left side" << endl;

		//строим правую часть СЛАУ основного метода Ньютона
		// находим значения основного оператора F
		GetOperatorF(source.numberSource, a, overline_a, b, xi, u, overline_u, Source_R, Source_X, F_part_odd, F_part_even);
		Lasting("Counting time of the main matrix", time);

		// перенумерация xi, u
		Renumbering(xi, numbered_xi);
		for (size_t count = 0; count < source.numberSource; ++count)
		{
			Renumbering(u[count], numbered_u[count]);
		}

		// находим F^{\prime}(x)*x и вычитаем F(x) 
		GetValueDerivedFunction(source.numberSource, numbered_xi, numbered_u, F_odd, F_even, F_part_odd, F_part_even);

		// находим окончательно правую часть b0, b1,...
		Getb(source.numberSource, F_odd, F_even, F_part_odd, F_part_even, b_right);
		Lasting("Calculation time right side", time);

		// Находим xi
		InvertMatrix(B, inverseMatrixB);
		GetXi(source.numberSource, A, inverseMatrixB, b_right, numbered_xi);
		Lasting("Time of xi", time);

		//находим u^{i}
		GetU(source.numberSource, A, inverseMatrixB, b_right, numbered_xi, numbered_u);
		Lasting("Time of u", time);

		// изменяем alpha для следующей итерации
		alpha = alpha * multiplier;

		// Обратная перенумерация u^{i}, xi
		InverseRenumbering(numbered_xi, xi);
		for (size_t count = 0; count < source.numberSource; ++count)
		{
			InverseRenumbering(numbered_u[count], u[count]);
		}

		// проекция Re(xi) >= 0
		ProjectionXi(xi);

		// печать результатов итераций в файл
		PrintXi(xi, iteration);
		Lasting("Calculation time solutions", time);
	}
	Lasting("The total time of the program", timeBegin);
}
