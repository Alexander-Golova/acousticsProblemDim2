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

	float alpha;
	cout << "Enter alpha ";
	cin >> alpha;

	float multiplier;
	cout << "Enter q ";
	cin >> multiplier;

	const Source source;

	const size_t N = NUMBER_PARTITION_POsize_tS;
	const size_t N_squared = (N + 1) * (N + 1);
	const float h = (float)DOMAIN_IN_HOMOGENEITY / N;

	// начало счета времени
	clock_t timeStart, timeFinish, timeBegin;
	timeBegin = clock();
	timeStart = clock();

	// выделение памяти
	// выделяем память под основные матрицы
	vector<vector<vector<vector<complex<float>>>>> a(N + 1,
		vector<vector<vector<complex<float>>>>(N + 1, vector<vector<complex<float>>>(N + 1,
			vector<complex<float>>(N + 1, complex<float>()))));

	vector<vector<vector<complex<float>>>> overline_a(N + 1, vector<vector<complex<float>>>(N + 1,
		vector<complex<float>>(N + 1, complex<float>())));

	vector<vector<complex<float>>> b(N + 1, vector<complex<float>>(N + 1, complex<float>()));

	// выделяем память для значений источников
	vector<vector<vector<complex<float>>>> Source_R(source.numberSource, vector<vector<complex<float>>>(N + 1, vector<complex<float>>(N + 1, complex<float>())));
	vector<vector<complex<float>>> Source_X(source.numberSource, vector<complex<float>>(N + 1, complex<float>()));

	// Выделяем память для поля в приемниках
	vector<vector<complex<float>>> overline_u(source.numberSource, vector<complex<float>>(N + 1, complex<float>()));

	//выделение памяти под массивы производных  F_1, F_2, ...
	vector<vector<vector<complex<float>>>> F_odd(source.numberSource + 1, vector<vector<complex<float>>>(N_squared, vector<complex<float>>(N_squared, complex<float>())));
	vector<vector<vector<complex<float>>>> F_even(source.numberSource + 1, vector<vector<complex<float>>>(N + 1, vector<complex<float>>(N_squared, complex<float>())));

	//выделение памяти под массивы A и B
	vector<vector<vector<complex<float>>>> A(source.numberSource + 1, vector<vector<complex<float>>>(N_squared, vector<complex<float>>(N_squared, complex<float>())));
	vector<vector<complex<float>>> B(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> inverseMatrixB(N_squared, vector<complex<float>>(N_squared, complex<float>()));
	
	// память для хранения значений основного оператора
	vector<vector<complex<float>>> F_part_odd(source.numberSource, vector<complex<float>>(N_squared, complex<float>()));
	vector<vector<complex<float>>> F_part_even(source.numberSource, vector<complex<float>>(N + 1, complex<float>()));

	// память для b_0, b_1,...
	vector<vector<complex<float>>> b_right(source.numberSource + 1, vector<complex<float>>(N_squared, complex<float>()));

	// память для u^(1), u^(2), u^(3)
	vector<vector<vector<complex<float>>>> u(source.numberSource + 1, vector<vector<complex<float>>>(N + 1, vector<complex<float>>(N + 1, complex<float>())));

	// память для xi
	vector<vector<complex<float>>> xi(N + 1, vector<complex<float>>(N + 1, complex<float>()));

	// память для перенумерованных переменных
	vector<vector<complex<float>>> numbered_u(source.numberSource, vector<complex<float>>(N_squared, complex<float>()));
	vector<complex<float>> numbered_xi(N_squared, complex<float>());

	timeFinish = clock();
	Lasting("Time allocation", timeStart, timeFinish);
	timeStart = clock();

	// Загрузка данных
	ArrayLoadingA(a);
	ArrayLoadingOverlineA(overline_a);
	ArrayLoadingB(b);
	ArrayLoadingSource(source.numberSource, Source_R, Source_X);
	ArrayLoadingOverlineU(source.numberSource, overline_u);

	//печатаем время работы
	timeFinish = clock();
	Lasting("Download time", timeStart, timeFinish);

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
		timeStart = clock();

		//строим левую часть СЛАУ основного метода Ньютона
		// получаем матрицы Якобиана
		GetJacobian(source.numberSource, a, overline_a, b, xi, u, F_odd, F_even);
		timeFinish = clock();
		Lasting("The counting time of the Jacobian matrices", timeStart, timeFinish);
		timeStart = clock();

		// находим матрицы А
		GetMatrixA(source.numberSource, F_odd, F_even, A, alpha);
		timeFinish = clock();
		Lasting("Counting time of matrices A", timeStart, timeFinish);
		timeStart = clock();

		// находим матрицу В
		GetMatrixB(F_odd, F_even, B, alpha);
		timeFinish = clock();
		Lasting("Counting time of matrices B", timeStart, timeFinish);
		timeStart = clock();

		cout << "Calculate the left side" << endl;

		//строим правую часть СЛАУ основного метода Ньютона
		// находим значения основного оператора F
		GetOperatorF(source.numberSource, a, overline_a, b, xi, u, overline_u, Source_R, Source_X, F_part_odd, F_part_even);
		timeFinish = clock();
		Lasting("Counting time of the main matrix", timeStart, timeFinish);
		timeStart = clock();

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
		timeFinish = clock();
		Lasting("Calculation time right side", timeStart, timeFinish);
		timeStart = clock();

		// Находим xi
		InvertMatrix(B, inverseMatrixB);
		GetXi(source.numberSource, A, inverseMatrixB, b_right, numbered_xi);
		timeFinish = clock();
		Lasting("Time of xi", timeStart, timeFinish);
		timeStart = clock();

		//находим u^{i}
		GetU(source.numberSource, A, inverseMatrixB, b_right, numbered_xi, numbered_u);
		timeFinish = clock();
		Lasting("Time of u", timeStart, timeFinish);
		timeStart = clock();

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
		Prsize_tXi(xi, iteration);
	}

	timeFinish = clock();
	Lasting("Calculation time solutions", timeBegin, timeFinish);

	return 0;
}
