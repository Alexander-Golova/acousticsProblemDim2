#pragma once
#include "../directProblemQuadratureSamplingMethodNew/taskData.h"
#include "../directProblemQuadratureSamplingMethodNew/Sources.h"
#include <vector>
#include <complex>

void InitialValueU(const size_t numberSource, std::vector<std::vector<std::vector<std::complex<double>>>> &u,
	std::vector<std::vector<std::vector<std::complex<double>>>> Source_R);

void InitialValueXi(std::vector<std::vector<std::complex<double>>> &xi);
