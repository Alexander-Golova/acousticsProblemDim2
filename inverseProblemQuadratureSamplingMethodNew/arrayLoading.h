#pragma once

#include <vector>
#include <complex>
#include "../directProblemQuadratureSamplingMethodNew/taskData.h"

void ArrayLoadingA(std::vector<std::vector<std::vector<std::vector<std::complex<float>>>>> & a);

void ArrayLoadingOverlineA(std::vector<std::vector<std::vector<std::complex<float>>>> & overline_a);

void ArrayLoadingB(std::vector<std::vector<std::complex<float>>> & b);

void ArrayLoadingSource(const size_t numberSource, std::vector<std::vector<std::vector<std::complex<float>>>> & Source_R, std::vector<std::vector<std::complex<float>>>Source_X);

void ArrayLoadingOverlineU(const size_t numberSource, std::vector<std::vector<std::complex<float>>> & overline_u);
