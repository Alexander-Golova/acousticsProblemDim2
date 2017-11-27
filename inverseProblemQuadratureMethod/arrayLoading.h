#pragma once

#include <vector>
#include <complex>
#include "../directProblemQuadratureMethod/taskData.h"

void ArrayLoadingA(std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> & a);

void ArrayLoadingOverlineA(std::vector<std::vector<std::vector<std::complex<double>>>> & overline_a);

void ArrayLoadingB(std::vector<std::vector<std::complex<double>>> & b);

void ArrayLoadingSource(const size_t numberSource,
	std::vector<std::vector<std::vector<std::complex<double>>>> & Source_R,
	std::vector<std::vector<std::complex<double>>> & Source_X);

void ArrayLoadingOverlineU(const size_t numberSource, std::vector<std::vector<std::complex<double>>> & overline_u);

void LoadData(const size_t numberSource, std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> & a,
	std::vector<std::vector<std::vector<std::complex<double>>>> & overline_a,
	std::vector<std::vector<std::complex<double>>> & b,
	std::vector<std::vector<std::vector<std::complex<double>>>> & Source_R,
	std::vector<std::vector<std::complex<double>>> & Source_X,
	std::vector<std::vector<std::complex<double>>> & overline_u);
