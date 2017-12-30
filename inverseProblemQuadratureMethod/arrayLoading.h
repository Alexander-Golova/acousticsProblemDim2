#pragma once

#include <vector>
#include <complex>
#include "../directProblemQuadratureMethod/taskData.h"

void ArrayLoadingA(std::vector<std::vector<std::vector<std::vector<float>>>> & a);
void ArrayLoadingB(std::vector<std::vector<std::vector<std::vector<float>>>> & b);
void ArrayLoadingC(std::vector<std::vector<float>> & c);

void ArrayLoadingOverlineA(std::vector<std::vector<std::vector<float>>> & overline_a);
void ArrayLoadingOverlineB(std::vector<std::vector<std::vector<float>>> & overline_b);

void ArrayLoadingB(std::vector<std::vector<std::complex<float>>> & b);

void ArrayLoadingSource(const size_t numberSource,
	std::vector<std::vector<std::vector<std::complex<float>>>> & Source_R,
	std::vector<std::vector<std::complex<float>>> & Source_X);

void ArrayLoadingOverlineU(const size_t numberSource, std::vector<std::vector<std::complex<float>>> & overline_u);

void LoadData(const size_t numberSource,
	std::vector<std::vector<std::vector<std::vector<float>>>> & a,
	std::vector<std::vector<std::vector<std::vector<float>>>> & b,
	std::vector<std::vector<float>> & c,
	std::vector<std::vector<std::vector<float>>> & overline_a,
	std::vector<std::vector<std::vector<float>>> & overline_b,
	std::vector<std::vector<std::vector<std::complex<float>>>> & Source_R,
	std::vector<std::vector<std::complex<float>>> & Source_X,
	std::vector<std::vector<std::complex<float>>> & overline_u);
