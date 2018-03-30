#pragma once

#include <vector>
#include <complex>
#include "../directProblemQuadratureMethod/taskData.h"

void ArrayLoadingA(std::vector<std::vector<std::vector<std::vector<std::complex<float>>>>> & a) noexcept;

void ArrayLoadingOverlineA(std::vector<std::vector<std::vector<std::complex<float>>>> & overline_a) noexcept;

void ArrayLoadingSource(const size_t numberSource,
	std::vector<std::vector<std::vector<std::complex<float>>>> & Source_R,
	std::vector<std::vector<std::complex<float>>> & Source_X) noexcept;

void ArrayLoadingOverlineU(const size_t numberSource, std::vector<std::vector<std::complex<float>>> & overline_u) noexcept;

void LoadData(const size_t numberSource,
	std::vector<std::vector<std::vector<std::vector<std::complex<float>>>>> & a,
	std::vector<std::vector<std::vector<std::complex<float>>>> & overline_a,
	std::vector<std::vector<std::vector<std::complex<float>>>> & Source_R,
	std::vector<std::vector<std::complex<float>>> & Source_X,
	std::vector<std::vector<std::complex<float>>> & overline_u) noexcept;
