#pragma once

#include <vector>
#include <complex>
#include "../directProblemQuadratureMethod/taskData.h"

void ArrayLoadingA(std::vector<std::vector<std::vector<std::vector<float>>>> & a) noexcept;
void ArrayLoadingB(std::vector<std::vector<std::vector<std::vector<float>>>> & b) noexcept;
void ArrayLoadingC(std::vector<std::vector<float>> & c) noexcept;

void ArrayLoadingOverlineA(std::vector<std::vector<std::vector<float>>> & overline_a) noexcept;
void ArrayLoadingOverlineB(std::vector<std::vector<std::vector<float>>> & overline_b) noexcept;

void ArrayLoadingB(std::vector<std::vector<std::complex<float>>> & b) noexcept;

void ArrayLoadingSource(const size_t numberSource,
	std::vector<std::vector<std::vector<std::complex<float>>>> & Source_R,
	std::vector<std::vector<std::complex<float>>> & Source_X) noexcept;

void ArrayLoadingOverlineU(const size_t numberSource, std::vector<std::vector<std::complex<float>>> & overline_u) noexcept;

void LoadData(const size_t numberSource,
	std::vector<std::vector<std::vector<std::vector<float>>>> & a,
	std::vector<std::vector<std::vector<std::vector<float>>>> & b,
	std::vector<std::vector<float>> & c,
	std::vector<std::vector<std::vector<float>>> & overline_a,
	std::vector<std::vector<std::vector<float>>> & overline_b,
	std::vector<std::vector<std::vector<std::complex<float>>>> & Source_R,
	std::vector<std::vector<std::complex<float>>> & Source_X,
	std::vector<std::vector<std::complex<float>>> & overline_u) noexcept;
