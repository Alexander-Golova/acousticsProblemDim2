#pragma once
#include "taskData.h"
#include <vector>
#include <complex>

// Основные матрицы задачи
/*
struct BasicArrays
{
	//const size_t NUMBER_PARTITION_POINT = 50;

	/*std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> a(NUMBER_PARTITION_POINT + 1,
		std::vector<std::vector<std::vector<std::complex<double>>>>(NUMBER_PARTITION_POINT + 1,
			std::vector<std::vector<std::complex<double>>>(NUMBER_PARTITION_POINT + 1,
				std::vector<std::complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>()))));

	vector<vector<vector<complex<double>>>> overline_a(NUMBER_PARTITION_POINT + 1,
		vector<vector<complex<double>>>(NUMBER_PARTITION_POINT + 1,
			vector<complex<double>>(NUMBER_PARTITION_POINT + 1, complex<double>())));

	std::vector<std::vector<std::complex<double>>> b(NUMBER_PARTITION_POINT + 1,
		std::vector<std::complex<double>>(NUMBER_PARTITION_POINT + 1, std::complex<double>()));

};
*/

void GetBasicArrays(std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> & a,
     std::vector<std::vector<std::vector<std::complex<double>>>> & overline_a,
     std::vector<std::vector<std::complex<double>>> & b);

